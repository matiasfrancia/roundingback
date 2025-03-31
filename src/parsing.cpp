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

#include "parsing.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <sstream>
#include "Solver.hpp"

namespace rs {

bigint parsing::read_number(const std::string& s) {
  bigint answer = 0;
  bool negate = false;
  for (char c : s) {
    if ('0' <= c && c <= '9') {
      answer *= 10;
      answer += c - '0';
    }
    negate = (negate || (c == '-'));
  }
  return negate ? -answer : answer;
}

bool parsing::opb_read(std::istream& in, Solver& solver, CeArb objective) {
  assert(objective->isReset());
  CeArb input = solver.cePools.takeArb();
  [[maybe_unused]] bool first_constraint = true;
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == '*') continue;
    for (char& c : line)
      if (c == ';') c = ' ';
    bool opt_line = line.substr(0, 4) == "min:";
    std::string line0;
    if (opt_line)
      line = line.substr(4), assert(first_constraint);
    else {
      std::string symbol;
      if (line.find(">=") != std::string::npos)
        symbol = ">=";
      else
        symbol = "=";
      assert(line.find(symbol) != std::string::npos);
      line0 = line;
      line = line.substr(0, line.find(symbol));
    }
    first_constraint = false;
    std::istringstream is(line);
    input->reset();
    std::vector<std::string> tokens;
    std::string tmp;
    while (is >> tmp) tokens.push_back(tmp);
    if (tokens.size() % 2 != 0) quit::exit_ERROR({"No support for non-linear constraints."});
    for (int i = 0; i < (long long)tokens.size(); i += 2)
      if (find(tokens[i].begin(), tokens[i].end(), 'x') != tokens[i].end())
        quit::exit_ERROR({"No support for non-linear constraints."});
    for (int i = 0; i < (long long)tokens.size(); i += 2) {
      std::string scoef = tokens[i];
      std::string var = tokens[i + 1];
      BigCoef coef = read_number(scoef);
      bool negated = false;
      std::string origvar = var;
      if (!var.empty() && var[0] == '~') {
        negated = true;
        var = var.substr(1);
      }
      if (var.empty() || var[0] != 'x') quit::exit_ERROR({"Invalid literal token: ", origvar});
      var = var.substr(1);
      Lit l = atoi(var.c_str());
      if (l < 1) quit::exit_ERROR({"Variable token less than 1: ", origvar});
      if (negated) l = -l;
      solver.setNbVars(std::abs(l), true);
      input->addLhs(coef, l);
    }

    // initialize the metrics
    if (stats.backboneMetrics.empty()) {
      stats.backboneMetrics.push_back({0.0, 0, 0, 0, 0.0, 0.0});
    }
    
    if (opt_line)
      input->copyTo(objective);
    else {
      input->addRhs(read_number(line0.substr(line0.find("=") + 1)));
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
      if (line0.find(" = ") != std::string::npos) {  // Handle equality case with second constraint
        input->invert();
        if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
          quit::exit_UNSAT(solver);
          return true;
        }
      }
    }
  }
  return false;
}

bool parsing::wcnf_read(std::istream& in, BigCoef top, Solver& solver, CeArb objective) {
  assert(objective->isReset());
  CeArb input = solver.cePools.takeArb();
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c')
      continue;
    else {
      std::istringstream is(line);
      std::string sweight;
      is >> sweight;
      BigCoef weight = read_number(sweight);
      if (weight == 0) continue;
      input->reset();
      input->addRhs(1);
      Lit l;
      while (is >> l, l) input->addLhs(1, l);
      if (weight < top) {  // soft clause
        std::stringstream s;
        s << "Negative clause weight: " << weight;
        if (weight < 0) quit::exit_ERROR({s.str()});
        solver.setNbVars(solver.getNbVars() + 1, true);  // increases n to n+1
        objective->addLhs(weight, solver.getNbVars());
        input->addLhs(1, solver.getNbVars());
      }  // else hard clause
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
    }
  }
  return false;
}

bool parsing::cnf_read(std::istream& in, Solver& solver) {
  Ce32 input = solver.cePools.take32();
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c')
      continue;
    else {
      std::istringstream is(line);
      input->reset();
      input->addRhs(1);
      Lit l;
      while (is >> l, l) {
        solver.setNbVars(std::abs(l), true);
        input->addLhs(1, l);
      }
      if (solver.addConstraint(input, Origin::FORMULA).second == ID_Unsat) {
        quit::exit_UNSAT(solver);
        return true;
      }
    }
  }
  return false;
}

bool parsing::file_read(boost::iostreams::filtering_istream& in, Solver& solver, CeArb objective) {
  for (std::string line; getline(in, line);) {
    if (line.empty() || line[0] == 'c') continue;
    if (line[0] == 'p') {
      std::istringstream is(line);
      is >> line;  // skip 'p'
      std::string type;
      is >> type;
      if (type == "cnf") {
        return cnf_read(in, solver);
      } else if (type == "wcnf") {
        long long val;
        is >> val;  // variable count
        solver.setNbVars(val);
        is >> line;  // skip nbConstraints
        std::string stop;
        is >> stop;  // top
        BigCoef top = read_number(stop);
        return wcnf_read(in, top, solver, objective);
      }
    } else if (line[0] == '*' && line.substr(0, 13) == "* #variable= ") {
      return opb_read(in, solver, objective);
    } else {
      quit::exit_ERROR({"No supported format [opb, cnf, wcnf] detected."});
      return true;
    }
  }
  return true;
}

}  // namespace rs
