/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2021, Jo Devriendt
Copyright (c) 2020-2021, Stephan Gocht
Copyright (c) 2014-2021, Jakob Nordström

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/


#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#ifdef IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif
#include <boost/iostreams/filtering_stream.hpp>

#include "run.hpp"
#include "ConstrExp.hpp"
#include "Solver.hpp"
#include "parsing.hpp"

namespace rs {

Solver run::solver; 

run::LazyVar::LazyVar(Solver& slvr, const Ce32 cardCore, int cardUpperBound, Var startVar)
    : solver(slvr), coveredVars(cardCore->getDegree()), upperBound(cardUpperBound) {
  assert(cardCore->isCardinality());
  cardCore->toSimple()->copyTo(atLeast);
  atLeast.toNormalFormLit();
  assert(atLeast.rhs == cardCore->getDegree());
  atMost.rhs = -atLeast.rhs;
  atMost.terms.reserve(atLeast.terms.size());
  for (auto& t : atLeast.terms) {
    atMost.terms.emplace_back(-t.c, t.l);
  }
  currentVar = startVar;
  atLeast.terms.emplace_back(-1, startVar);
  atMost.terms.emplace_back(remainingVars(), startVar);
  ++coveredVars;
}

run::LazyVar::~LazyVar() {
  solver.dropExternal(atLeastID, false, false);
  solver.dropExternal(atMostID, false, false);
}

int run::LazyVar::remainingVars() const { return upperBound - coveredVars; }

void run::LazyVar::setUpperBound(int cardUpperBound) {
  assert(upperBound >= cardUpperBound);
  upperBound = cardUpperBound;
}

void run::LazyVar::addVar(Var v, bool reified) {
  currentVar = v;
  if (reified) {
    Term<int>& last = atLeast.terms.back();
    last = {last.c - 1, v};
    --atMost.rhs;
    Term<int>& last2 = atMost.terms.back();
    last2 = {remainingVars(), v};
  } else {
    atLeast.terms.emplace_back(-1, v);
    Term<int>& last = atMost.terms.back();
    last = {1, last.l};
    atMost.terms.emplace_back(remainingVars(), v);
  }
  ++coveredVars;
}

ID run::LazyVar::addAtLeastConstraint(bool reified) {
  assert(atLeast.terms.back().l == currentVar);
  solver.dropExternal(atLeastID, !reified, false);
  atLeastID = solver.addConstraint(atLeast, Origin::COREGUIDED).second;
  return atLeastID;
}

ID run::LazyVar::addAtMostConstraint(bool reified) {
  assert(atMost.terms.back().l == currentVar);
  solver.dropExternal(atMostID, !reified, false);
  atMostID = solver.addConstraint(atMost, Origin::COREGUIDED).second;
  return atMostID;
}

ID run::LazyVar::addSymBreakingConstraint(Var prevvar) const {
  assert(prevvar < currentVar);
  // y-- + ~y >= 1 (equivalent to y-- >= y)
  return solver.addConstraint(ConstrSimple32({{1, prevvar}, {1, -currentVar}}, 1), Origin::COREGUIDED).second;
}

ID run::LazyVar::addFinalAtMost(bool reified) {
  solver.dropExternal(atMostID, !reified, false);
  Term<int>& last = atMost.terms.back();
  last = {1, last.l};
  atMostID = solver.addConstraint(atMost, Origin::COREGUIDED).second;
  return atMostID;
}

std::ostream& run::operator<<(std::ostream& o, const std::shared_ptr<LazyVar> lv) {
  o << lv->atLeast << "\n" << lv->atMost;
  return o;
}

void run::decide() {

  while (true) {
    SolveState reply = aux::timeCall<SolveState>([&] { return solver.solve().state; }, stats.SOLVETIME);
    assert(reply != SolveState::INCONSISTENT);
    switch (reply) {
      case SolveState::SAT:
        std::cout << "SAT\n";
        break;
      case SolveState::UNSAT:
        std::cout << "UNSAT\n";
        break;
      case SolveState::INPROCESSED:
        std::cout << "INPROCESSED\n";
        break;
      default:
        std::cout << "UNKNOWN\n";
        break;
    }
    if (options.calculateBackbones && reply == SolveState::SAT) {
      solver.calculateBackbones();
      quit::exit_SAT(solver);
      return;
    }
    else if (reply == SolveState::SAT) {
      quit::exit_SAT(solver);
      return;
    } else if (reply == SolveState::UNSAT) {
      quit::exit_UNSAT(solver);
      return;
    }
    if (options.timeout.get() != -1.0 && stats.getTime() > options.timeout.get()) {
      quit::exit_INDETERMINATE(solver);
      return;
    }
  }
}

void run::run(CeArb objective) {
  stats.RUNSTARTTIME = aux::cpuTime();
  if (options.verbosity.get() > 0)
    std::cout << "c #variables " << solver.getNbOrigVars() << " #constraints " << solver.getNbConstraints()
              << std::endl;
  try {
    if (objective->vars.size() > 0) {
      if (options.alreadySolved) {
        std::cout << "c Already solved\n";
        std::cout << "c Optimum value: " << options.optimumValue.get() << std::endl;
        
        CeArb upperBoundCons = solver.cePools.takeArb();

        objective->copyTo(upperBoundCons);
        upperBoundCons->invert();
        upperBoundCons->addRhs(-options.optimumValue.get());
        std::pair<ID, ID> res = solver.addConstraint(upperBoundCons, Origin::FORMULA);
        assert(res.second != ID_Unsat);

        CeArb lowerBoundCons = solver.cePools.takeArb();
        objective->copyTo(lowerBoundCons);
        lowerBoundCons->addRhs(options.optimumValue.get());
        std::pair<ID, ID> res2 = solver.addConstraint(lowerBoundCons, Origin::FORMULA);
        assert(res2.second != ID_Unsat);
        objective->reset();

        decide();
        return;
      }
      objective->stopLogging();
      objective->removeUnitsAndZeroes(solver.getLevel(), solver.getPos(), false);
      bigint bestUpperBound;

      BigVal maxVal = objective->getCutoffVal();
      if (maxVal <= limit32) {  // TODO: try to internalize this check in ConstrExp
        Ce32 result = solver.cePools.take32();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
        bestUpperBound = optim.getUpperBound();
      } else if (maxVal <= limit64) {
        Ce64 result = solver.cePools.take64();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
        bestUpperBound = optim.getUpperBound();
      } else if (maxVal <= BigVal(limit96)) {
        Ce96 result = solver.cePools.take96();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
        bestUpperBound = optim.getUpperBound();
      } else if (maxVal <= BigVal(limit128)) {
        Ce128 result = solver.cePools.take128();
        objective->copyTo(result);
        Optimization optim(result);
        optim.optimize();
        bestUpperBound = optim.getUpperBound();
      } else {
        CeArb result = solver.cePools.takeArb();
        objective->copyTo(result);
        Optimization<bigint, bigint> optim(result);
        optim.optimize();
        bestUpperBound = optim.getUpperBound();
      }

      // if (rs::options.calculateBackbones) {
      //   std::cout << "c Calculating backbones\n";
      //   std::cout << "Best objective value: " << bestUpperBound << std::endl;

      //   // reinitialize solver for backbone calculation
      //   rs::run::solver.reset();                              // call reset
      //   rs::run::solver.~Solver();                    // call destructor
      //   new (&rs::run::solver) Solver(/* optional args */);  // placement new

      //   rs::run::solver.init();
      //   rs::CeArb objective = rs::run::solver.cePools.takeArb();
      //   bool infeasible_or_error = false;

      //   // read the formula again
      //   if (!rs::options.formulaName.empty()) {
      //     std::ifstream fin(rs::options.formulaName, std::ifstream::in);
      //     boost::iostreams::filtering_istream in;

      // #ifdef IOSTREAMS_WITH_ZLIB
      //     if (boost::algorithm::ends_with(rs::options.formulaName, ".gz"))
      //       in.push(boost::iostreams::gzip_decompressor());
      // #endif
      // #ifdef IOSTREAMS_WITH_BZIP2
      //     if (boost::algorithm::ends_with(rs::options.formulaName, ".bz2"))
      //       in.push(boost::iostreams::bzip2_decompressor());
      // #endif
      //     in.push(fin);
      //     infeasible_or_error = rs::parsing::file_read(in, rs::run::solver, objective);
      //   }

      //   rs::run::solver.initLP(objective);

      //   std::cout << "c Adding upper bound constraint\n";

      //   // add upper bound constraint
      //   CeArb upperBound = rs::run::solver.cePools.takeArb();
      //   objective->copyTo(upperBound);
      //   upperBound->invert();
      //   upperBound->addRhs(-bestUpperBound);
      //   std::pair<ID, ID> res = rs::run::solver.addConstraint(upperBound, Origin::UPPERBOUND);
      //   assert(res.second != ID_Unsat);

      //   std::cout << "c Starting backbone calculation\n";

      //   decide();
      // }
    } else {
      if (options.alreadySolved) {
        std::cout << "Warning: alreadySolved is set but no objective function is given. " <<
          "Solving the instance without taking into consideration the optimum-value passed as parameter\n" << std::endl;
      }
      decide();
    }
  } catch (const AsynchronousInterrupt& ai) {
    std::cout << "c " << ai.what() << std::endl;
    quit::exit_INDETERMINATE(solver);
  }
}

}  // namespace rs
