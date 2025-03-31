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

#pragma once

#include "Solver.hpp"
#include "typedefs.hpp"

namespace rs {

namespace run {

extern Solver solver;

struct LazyVar {
  Solver& solver;
  int coveredVars;
  int upperBound;
  Var currentVar;
  ID atLeastID = ID_Undef;
  ID atMostID = ID_Undef;
  ConstrSimple32 atLeast;  // X >= k + y1 + ... + yi
  ConstrSimple32 atMost;   // k + y1 + ... + yi-1 + (1+n-k-i)yi >= X

  LazyVar(Solver& slvr, const Ce32 cardCore, int cardUpperBound, Var startVar);
  ~LazyVar();

  void addVar(Var v, bool reified);
  ID addAtLeastConstraint(bool reified);
  ID addAtMostConstraint(bool reified);
  ID addSymBreakingConstraint(Var prevvar) const;
  ID addFinalAtMost(bool reified);
  int remainingVars() const;
  void setUpperBound(int cardUpperBound);
};

std::ostream& operator<<(std::ostream& o, const std::shared_ptr<LazyVar> lv);

template <typename SMALL>
struct LvM {
  std::unique_ptr<LazyVar> lv;
  SMALL m;
};

template <typename SMALL, typename LARGE>
class Optimization {
  const CePtr<ConstrExp<SMALL, LARGE>> origObj;
  CePtr<ConstrExp<SMALL, LARGE>> reformObj;

  LARGE lower_bound;
  LARGE upper_bound;
  ID lastUpperBound = ID_Undef;
  ID lastUpperBoundUnprocessed = ID_Undef;
  ID lastLowerBound = ID_Undef;
  ID lastLowerBoundUnprocessed = ID_Undef;

  std::vector<LvM<SMALL>> lazyVars;
  IntSet assumps;

 public:
  Optimization(CePtr<ConstrExp<SMALL, LARGE>> obj) : origObj(obj) {
    assert(origObj->vars.size() > 0);
    // NOTE: -origObj->getDegree() keeps track of the offset of the reformulated objective (or after removing unit lits)
    lower_bound = -origObj->getDegree();
    upper_bound = origObj->absCoeffSum() - origObj->getRhs() + 1;

    reformObj = solver.cePools.take<SMALL, LARGE>();
    reformObj->stopLogging();
    if (!options.optMode.is("linear")) {
      origObj->copyTo(reformObj);
    }
  };

  LARGE normalizedLowerBound() { return lower_bound + origObj->getDegree(); }
  LARGE normalizedUpperBound() { return upper_bound + origObj->getDegree(); }

  LARGE getUpperBound() const { return upper_bound; }
  LARGE getLowerBound() const { return lower_bound; }

  void printObjBounds() {
    if (options.verbosity.get() == 0) return;
    std::cout << "c bounds ";
    if (solver.foundSolution()) {
      std::cout << bigint(upper_bound);  // TODO: remove bigint(...) hack
    } else {
      std::cout << "-";
    }
    std::cout << " >= " << bigint(lower_bound) << " @ " << stats.getTime() << "\n";  // TODO: remove bigint(...) hack
  }

  bool checkLazyVariables() {
    bool reified = options.cgEncoding.is("reified");
    for (int i = 0; i < (int)lazyVars.size(); ++i) {
      LazyVar& lv = *lazyVars[i].lv;
      if (reformObj->getLit(lv.currentVar) == 0) {
        int cardCoreUpper = lv.upperBound;
        if (options.cgCoreUpper) {
          cardCoreUpper = static_cast<int>(std::min<LARGE>(cardCoreUpper, normalizedUpperBound() / lazyVars[i].m));
          // NOTE: integer division rounds down
        }
        assert(cardCoreUpper >= 0);
        lv.setUpperBound(cardCoreUpper);
        if (lv.remainingVars() == 0 ||
            isUnit(solver.getLevel(), -lv.currentVar)) {  // binary constraints make all new auxiliary variables unit
          if (lv.addFinalAtMost(reified) == ID_Unsat) {
            quit::exit_UNSAT(solver, upper_bound);
            return true;
          }
          aux::swapErase(lazyVars, i--);  // fully expanded, no need to keep in memory
        } else {                          // add auxiliary variable
          long long newN = solver.getNbVars() + 1;
          solver.setNbVars(newN);
          Var oldvar = lv.currentVar;
          lv.addVar(newN, reified);
          // reformulate the objective
          reformObj->addLhs(lazyVars[i].m, newN);
          // add necessary lazy constraints
          if (lv.addAtLeastConstraint(reified) == ID_Unsat || lv.addAtMostConstraint(reified) == ID_Unsat ||
              lv.addSymBreakingConstraint(oldvar) == ID_Unsat) {
            quit::exit_UNSAT(solver, upper_bound);
            return true;
          } else if (lv.remainingVars() == 0) {
            aux::swapErase(lazyVars, i--);  // fully expanded, no need to keep in memory
          }
        }
      }
    }
    return false;
  }

  bool addLowerBound() {
    CePtr<ConstrExp<SMALL, LARGE>> aux = solver.cePools.take<SMALL, LARGE>();
    origObj->copyTo(aux);
    aux->addRhs(lower_bound);
    solver.dropExternal(lastLowerBound, true, true);
    std::pair<ID, ID> res = solver.addConstraint(aux, Origin::LOWERBOUND);
    lastLowerBoundUnprocessed = res.first;
    lastLowerBound = res.second;
    if (lastLowerBound == ID_Unsat) {
      quit::exit_UNSAT(solver, upper_bound);
      return true;
    }
    return false;
  }

  std::pair<Ce32, int> reduceToCardinality(const CeSuper& core) {
    int bestNbBlocksRemoved = 0;
    CeSuper card = core->clone(solver.cePools);
    if (options.cgReduction.is("clause")) {
      card->sortInDecreasingCoefOrder(
          [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
      card->simplifyToClause();
    } else if (options.cgReduction.is("minslack")) {
      card->sortInDecreasingCoefOrder(
          [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
      card->simplifyToMinLengthCardinality();
    } else {
      assert(options.cgReduction.is("bestbound"));
      CeSuper cloneCoefOrder = card->clone(solver.cePools);
      cloneCoefOrder->sortInDecreasingCoefOrder();
      std::reverse(cloneCoefOrder->vars.begin(), cloneCoefOrder->vars.end());  // *IN*creasing coef order
      card->sort([&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) > aux::abs(reformObj->coefs[v2]); });
      CeSuper clone = card->clone(solver.cePools);
      assert(clone->vars.size() > 0);
      LARGE bestLowerBound = 0;
      int bestNbVars = clone->vars.size();

      // find the optimum number of variables to weaken to
      int nbBlocksRemoved = 0;
      while (!clone->isTautology()) {
        int carddegree = cloneCoefOrder->getCardinalityDegreeWithZeroes();
        if (bestLowerBound < aux::abs(reformObj->coefs[clone->vars.back()]) * carddegree) {
          bestLowerBound = aux::abs(reformObj->coefs[clone->vars.back()]) * carddegree;
          bestNbVars = clone->vars.size();
          bestNbBlocksRemoved = nbBlocksRemoved;
        }
        SMALL currentObjCoef = aux::abs(reformObj->coefs[clone->vars.back()]);
        // weaken all lowest objective coefficient literals
        while (clone->vars.size() > 0 && currentObjCoef == aux::abs(reformObj->coefs[clone->vars.back()])) {
          ++nbBlocksRemoved;
          Var v = clone->vars.back();
          clone->weakenLast();
          cloneCoefOrder->weaken(v);
        }
      }

      // weaken to the optimum number of variables and generate cardinality constraint
      while ((int)card->vars.size() > bestNbVars) {
        card->weakenLast();
      }
      card->saturate();
      // sort in decreasing coef order to minimize number of auxiliary variables, but break ties so that *large*
      // objective coefficient literals are removed first, as the smallest objective coefficients are the ones that
      // will be eliminated from the objective function. Thanks Stephan!
      card->sortInDecreasingCoefOrder(
          [&](Var v1, Var v2) { return aux::abs(reformObj->coefs[v1]) < aux::abs(reformObj->coefs[v2]); });
      card->simplifyToCardinality(false);
    }

    Ce32 result = solver.cePools.take32();
    card->copyTo(result);
    return {result, bestNbBlocksRemoved};
  }

  bool handleInconsistency(std::vector<CeSuper>& cores) {
    // take care of derived unit lits
    for (Var v : reformObj->vars) {
      if (isUnit(solver.getLevel(), v) || isUnit(solver.getLevel(), -v)) {
        assumps.remove(v);
        assumps.remove(-v);
      }
    }
    reformObj->removeUnitsAndZeroes(solver.getLevel(), solver.getPos(), false);
    [[maybe_unused]] LARGE prev_lower = lower_bound;
    lower_bound = -reformObj->getDegree();

    std::vector<std::pair<Ce32, int>> cardCores;
    for (CeSuper& core : cores) {
      assert(core);
      if (core->isTautology()) {
        continue;
      }
      if (core->hasNegativeSlack(solver.getLevel())) {
        assert(solver.decisionLevel() == 0);
        if (solver.logger) core->logInconsistency(solver.getLevel(), solver.getPos(), stats);
        quit::exit_UNSAT(solver, upper_bound);
        return true;
      }
      // figure out an appropriate core
      cardCores.push_back(reduceToCardinality(core));
    }

    if (cardCores.size() == 0) {
      // only violated unit assumptions were derived
      ++stats.UNITCORES;
      assert(lower_bound > prev_lower);
      if (checkLazyVariables()) return true;
      if (addLowerBound()) return true;
      if (!options.cgIndCores) assumps.clear();
      return false;
    }

    stats.SINGLECORES += cardCores.size() == 1;

    LARGE bestLowerBound = -1;
    Ce32& bestCardCore = cardCores[0].first;
    int bestBlocksRemoved = 0;
    for (int i = 0; i < (int)cardCores.size(); ++i) {
      Ce32 cardCore = cardCores[i].first;
      assert(cardCore->hasNoZeroes());
      assert(cardCore->vars.size() > 0);
      SMALL lowestCoef = aux::abs(reformObj->coefs[cardCore->vars[0]]);
      for (Var v : cardCore->vars) {
        if (aux::abs(reformObj->coefs[v]) < lowestCoef) lowestCoef = aux::abs(reformObj->coefs[v]);
      }
      LARGE lowerBound = lowestCoef * cardCore->degree;
      if (i == 1) {
        stats.NOCOREBEST += lowerBound == bestLowerBound;
        stats.FIRSTCOREBEST += lowerBound < bestLowerBound;
        stats.DECCOREBEST += lowerBound > bestLowerBound;
      }
      if (lowerBound > bestLowerBound) {
        bestLowerBound = lowerBound;
        bestCardCore = cardCore;
        bestBlocksRemoved = cardCores[i].second;
      }
    }

    stats.REMOVEDBLOCKS += bestBlocksRemoved;
    stats.NCORECARDINALITIES += !bestCardCore->isClause();
    stats.COREDEGSUM += bestCardCore->getDegree();
    stats.CORESLACKSUM += bestCardCore->vars.size() - bestCardCore->getDegree();

    for (Var v : bestCardCore->vars) {
      assert(assumps.has(-bestCardCore->getLit(v)));
      assumps.remove(-bestCardCore->getLit(v));  // independent cores
    }

    // adjust the lower bound
    SMALL mult = 0;
    for (Var v : bestCardCore->vars) {
      assert(reformObj->getLit(v) != 0);
      if (mult == 0) {
        mult = aux::abs(reformObj->coefs[v]);
      } else if (mult == 1) {
        break;
      } else {
        mult = std::min(mult, aux::abs(reformObj->coefs[v]));
      }
    }
    assert(mult > 0);
    lower_bound += bestCardCore->getDegree() * mult;

    int cardCoreUpper = bestCardCore->vars.size();
    if (options.cgCoreUpper) {
      cardCoreUpper = static_cast<int>(std::min<LARGE>(cardCoreUpper, normalizedUpperBound() / mult));
      // NOTE: integer division rounds down
    }
    assert(cardCoreUpper >= 0);

    if (options.cgEncoding.is("sum") || bestCardCore->vars.size() - bestCardCore->getDegree() <= 1) {
      // add auxiliary variables
      long long oldN = solver.getNbVars();
      long long newN = oldN - static_cast<int>(bestCardCore->getDegree()) + cardCoreUpper;
      solver.setNbVars(newN);
      // reformulate the objective
      for (Var v = oldN + 1; v <= newN; ++v) bestCardCore->addLhs(-1, v);
      bestCardCore->invert();
      reformObj->addUp(bestCardCore, mult);
      assert(lower_bound == -reformObj->getDegree());
      // add channeling constraints
      if (solver.addConstraint(bestCardCore, Origin::COREGUIDED).second == ID_Unsat) {
        quit::exit_UNSAT(solver, upper_bound);
        return true;
      }
      bestCardCore->invert();
      if (solver.addConstraint(bestCardCore, Origin::COREGUIDED).second == ID_Unsat) {
        quit::exit_UNSAT(solver, upper_bound);
        return true;
      }
      for (Var v = oldN + 1; v < newN; ++v) {  // add symmetry breaking constraints
        if (solver.addConstraint(ConstrSimple32({{1, v}, {1, -v - 1}}, 1), Origin::COREGUIDED).second == ID_Unsat) {
          quit::exit_UNSAT(solver, upper_bound);
          return true;
        }
      }
    } else {
      bool reified = options.cgEncoding.is("reified");
      // add auxiliary variable
      long long newN = solver.getNbVars() + 1;
      solver.setNbVars(newN);
      // reformulate the objective
      bestCardCore->invert();
      reformObj->addUp(bestCardCore, mult);
      bestCardCore->invert();
      reformObj->addLhs(mult, newN);  // add only one variable for now
      assert(lower_bound == -reformObj->getDegree());
      // add first lazy constraint
      lazyVars.push_back({std::make_unique<LazyVar>(solver, bestCardCore, cardCoreUpper, newN), mult});
      lazyVars.back().lv->addAtLeastConstraint(reified);
      lazyVars.back().lv->addAtMostConstraint(reified);
    }
    if (checkLazyVariables()) return true;
    if (addLowerBound()) return true;
    if (!options.cgIndCores) assumps.clear();
    return false;
  }

  bool handleNewSolution(const std::vector<Lit>& sol) {
    [[maybe_unused]] LARGE prev_val = upper_bound;
    upper_bound = -origObj->getRhs();
    for (Var v : origObj->vars) upper_bound += origObj->coefs[v] * (int)(sol[v] > 0);
    assert(upper_bound < prev_val);

    CePtr<ConstrExp<SMALL, LARGE>> aux = solver.cePools.take<SMALL, LARGE>();
    origObj->copyTo(aux);
    aux->invert();
    aux->addRhs(-upper_bound + 1);
    solver.dropExternal(lastUpperBound, true, true);
    std::pair<ID, ID> res = solver.addConstraint(aux, Origin::UPPERBOUND);
    lastUpperBoundUnprocessed = res.first;
    lastUpperBound = res.second;
    if (lastUpperBound == ID_Unsat) {
      quit::exit_UNSAT(solver, upper_bound);
      return true;
    }
    return false;
  }

  void logProof() {
    if (!solver.logger) return;
    assert(lastUpperBound != ID_Undef);
    assert(lastUpperBound != ID_Unsat);
    assert(lastLowerBound != ID_Undef);
    assert(lastLowerBound != ID_Unsat);
    CePtr<ConstrExp<SMALL, LARGE>> coreAggregate = solver.cePools.take<SMALL, LARGE>();
    CePtr<ConstrExp<SMALL, LARGE>> aux = solver.cePools.take<SMALL, LARGE>();
    origObj->copyTo(aux);
    aux->invert();
    aux->addRhs(1 - upper_bound);
    aux->resetBuffer(lastUpperBoundUnprocessed);
    coreAggregate->addUp(aux);
    aux->reset();
    origObj->copyTo(aux);
    aux->addRhs(lower_bound);
    aux->resetBuffer(lastLowerBoundUnprocessed);
    coreAggregate->addUp(aux);
    assert(coreAggregate->hasNegativeSlack(solver.getLevel()));
    assert(solver.decisionLevel() == 0);
    coreAggregate->logInconsistency(solver.getLevel(), solver.getPos(), stats);
  }

  bool harden() {
    LARGE diff = upper_bound - lower_bound;
    for (Var v : reformObj->vars) {
      if (aux::abs(reformObj->coefs[v]) > diff &&
          solver.addUnitConstraint(-reformObj->getLit(v), Origin::HARDENEDBOUND).second == ID_Unsat) {
        quit::exit_UNSAT(solver, upper_bound);
        return true;
      }
    }
    return false;
  }

  void optimize() {
    size_t upper_time = 0, lower_time = 0;
    SolveState reply = SolveState::SAT;
    SMALL coeflim = options.cgStrat ? reformObj->getLargestCoef() : 0;

    enum class CoefLimStatus { START, ENCOUNTEREDSAT, REFINE };
    CoefLimStatus coefLimFlag = CoefLimStatus::REFINE;
    while (true) {
      size_t current_time = stats.getDetTime();
      if (options.timeout.get() != -1.0 && stats.getTime() > options.timeout.get()) {
        std::cout << "time limit reached\n";
        return;
      }
      if (reply != SolveState::INPROCESSED && reply != SolveState::RESTARTED) printObjBounds();

      // NOTE: only if assumptions are empty will they be refilled
      if (assumps.isEmpty() &&
          (options.optMode.is("coreguided") ||
           (options.optMode.is("coreboosted") && stats.getRunTime() < options.cgBoosted.get()) ||
           (options.optMode.is("hybrid") &&
            lower_time <
                options.cgHybrid.get() * (upper_time + lower_time)))) {  // use core-guided step by setting assumptions
        reformObj->removeZeroes();
        if (coeflim > 0 && coefLimFlag == CoefLimStatus::REFINE) {
          // find a new coeflim
          reformObj->sortInDecreasingCoefOrder();
          int numLargerCoefs = 0;
          int numLargerUniqueCoefs = 0;
          SMALL prevCoef = -1;
          bool changed = false;
          for (Var v : reformObj->vars) {
            SMALL coef = aux::abs(reformObj->coefs[v]);
            ++numLargerCoefs;
            numLargerUniqueCoefs += prevCoef != coef;
            prevCoef = coef;
            if (coeflim > coef && numLargerCoefs / numLargerUniqueCoefs > 1.25) {
              changed = true;
              coeflim = coef;
              break;
            }
          }
          if (!changed) {
            coeflim = 0;
          }
        }
        std::vector<Term<double>> litcoefs;  // using double will lead to faster sort than bigint
        litcoefs.reserve(reformObj->vars.size());
        for (Var v : reformObj->vars) {
          assert(reformObj->getLit(v) != 0);
          SMALL cf = aux::abs(reformObj->coefs[v]);
          if (cf >= coeflim) litcoefs.emplace_back(static_cast<double>(cf), v);
        }
        std::sort(litcoefs.begin(), litcoefs.end(), [](const Term<double>& t1, const Term<double>& t2) {
          return t1.c > t2.c || (t1.l < t2.l && t1.c == t2.c);
        });
        for (const Term<double>& t : litcoefs) assumps.add(-reformObj->getLit(t.l));
        coefLimFlag = CoefLimStatus::ENCOUNTEREDSAT;
      }
      assert(upper_bound > lower_bound);
      // std::cout << "c nAssumptions for solve: " << assumps.size() << " @ " << stats.getTime() << " s\n";
      SolveAnswer out = aux::timeCall<SolveAnswer>(
          [&] {
            solver.setAssumptions(assumps);
            return solver.solve();
          },
          assumps.isEmpty() ? stats.SOLVETIME : stats.SOLVETIMECG);
      reply = out.state;
      if (reply == SolveState::RESTARTED) continue;
      if (reply == SolveState::UNSAT) {
        printObjBounds();
        quit::exit_UNSAT(solver, upper_bound);
        return;
      }
      if (assumps.isEmpty()) {
        upper_time += stats.getDetTime() - current_time;
      } else {
        lower_time += stats.getDetTime() - current_time;
      }
      if (reply == SolveState::SAT) {
        assumps.clear();
        assert(solver.foundSolution());
        ++stats.NSOLS;
        if (handleNewSolution(out.solution)) return;
        if (harden()) return;
        if (coefLimFlag == CoefLimStatus::ENCOUNTEREDSAT) coefLimFlag = CoefLimStatus::REFINE;
      } else if (reply == SolveState::INCONSISTENT) {
        assert(!options.optMode.is("linear"));
        ++stats.NCORES;
        if (handleInconsistency(out.cores)) return;
        if (harden()) return;
        coefLimFlag = CoefLimStatus::START;
      } else {
        assert(reply == SolveState::INPROCESSED);  // keep looping
        // coefLimFlag = CoefLimStatus::REFINE;
        // NOTE: above avoids long non-terminating solve calls due to insufficient assumptions
        assumps.clear();
      }
      if (lower_bound >= upper_bound) {
        printObjBounds();
        logProof();
        quit::exit_UNSAT(solver, upper_bound);
        return;
      }
    }
  }
};

void decide();
void run(CeArb objective);

}  // namespace run

}  // namespace rs
