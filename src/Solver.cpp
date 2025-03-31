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

#include "Solver.hpp"
#include <cmath>
#include "Constr.hpp"
#include "auxiliary.hpp"
#include "globals.hpp"

namespace rs {

// ---------------------------------------------------------------------
// Initialization

Solver::Solver() : order_heap(activity) {
  ca.capacity(1024 * 1024);  // 4MB
}

void Solver::setNbVars(long long nvars, bool orig) {
  assert(nvars > 0);
  assert(nvars < INF);
  if (nvars <= n) return;
  aux::resizeIntMap(_adj, adj, nvars, resize_factor, {});
  aux::resizeIntMap(_Level, Level, nvars, resize_factor, INF);
  Pos.resize(nvars + 1, INF);
  Reason.resize(nvars + 1, CRef_Undef);
  activity.resize(nvars + 1, 1 / actLimitV);
  phase.resize(nvars + 1);
  cePools.resize(nvars + 1);
  order_heap.resize(nvars + 1);
  for (Var v = n + 1; v <= nvars; ++v) phase[v] = -v, order_heap.insert(v);
  // if (lpSolver) lpSolver->setNbVariables(nvars + 1); // Currently, LP solver only reasons on formula constraints
  n = nvars;
  if (orig) {
    orig_n = n;
    stats.NORIGVARS = n;
    stats.NAUXVARS = 0;
    state.resize(n + 1, 0);
  } else {
    stats.NAUXVARS = n - orig_n;
  }
}

void Solver::init() {
  if (!options.proofLog.get().empty()) logger = std::make_shared<Logger>(options.proofLog.get());
  cePools.initializeLogging(logger);
}

void Solver::initLP([[maybe_unused]] const CeArb objective) {
#if WITHSOPLEX
  if (options.lpPivotRatio.get() == 0) return;
  bool pureCNF = objective->vars.size() == 0;
  for (CRef cr : constraints) {
    if (!pureCNF) break;
    pureCNF = (ca[cr].degree() == 1);
  }
  if (pureCNF) return;
  lpSolver = std::make_shared<LpSolver>(*this, objective);
#else
  return;
#endif  // WITHSOPLEX
}

void Solver::reset() {
  ca.~ConstraintAllocator();
  new (&ca) ConstraintAllocator();
}

// ---------------------------------------------------------------------
// VSIDS

void Solver::restartVSIDS() {
  v_vsids_inc = 1.0;
  for (Var x = 1; x <= n; ++x) {
    activity[x] = 1 / actLimitV;
  }
}

void Solver::vDecayActivity() { v_vsids_inc *= (1 / options.varDecay.get()); }
void Solver::vBumpActivity(Var v) {
  assert(v > 0);
  if ((activity[v] += v_vsids_inc) > actLimitV) {  // Rescale
    for (Var x = 1; x <= n; ++x) {
      activity[x] /= actLimitV;
      activity[x] /= actLimitV;
    }
    v_vsids_inc /= actLimitV;
    v_vsids_inc /= actLimitV;
  }
  // Update heap with respect to new activity:
  if (order_heap.inHeap(v)) order_heap.percolateUp(v);
}

void Solver::cDecayActivity() { c_vsids_inc *= (1 / options.clauseDecay.get()); }
void Solver::cBumpActivity(Constr& c) {
  c.act += c_vsids_inc;
  if (c.act > actLimitC) {  // Rescale:
    for (size_t i = 0; i < constraints.size(); i++) {
      ca[constraints[i]].act /= actLimitC;
      ca[constraints[i]].act /= actLimitC;
    }
    c_vsids_inc /= actLimitC;
    c_vsids_inc /= actLimitC;
  }
}

// ---------------------------------------------------------------------
// Assignment manipulation

void Solver::uncheckedEnqueue(Lit p, CRef from) {
  assert(!isTrue(Level, p));
  assert(!isFalse(Level, p));
  assert(isUnknown(Pos, p));
  Var v = toVar(p);
  Reason[v] = from;
  if (decisionLevel() == 0) {
    Reason[v] = CRef_Undef;  // no need to keep track of reasons for unit literals
    if (logger) {
      Constr& C = ca[from];
      C.toExpanded(cePools)->logUnit(Level, Pos, v, stats);
      assert(logger->unitIDs.size() == trail.size() + 1);
    }
  }
  Level[p] = decisionLevel();
  Pos[v] = (int)trail.size();
  trail.push_back(p);
}

void Solver::removeLastAssignment() {
  assert(!trail.empty());
  ++stats.NTRAILPOPS;
  Lit l = trail.back();
  if (qhead == (int)trail.size()) {
    for (const Watch& w : adj[-l])
      if (w.idx >= INF) ca[w.cref].undoFalsified(w.idx);
    --qhead;
  }
  Var v = toVar(l);
  trail.pop_back();
  Level[l] = INF;
  Pos[v] = INF;
  phase[v] = l;
  if (!trail_lim.empty() && trail_lim.back() == (int)trail.size()) trail_lim.pop_back();
  order_heap.insert(v);
}

void Solver::backjumpTo(int level) {
  assert(level >= 0);
  while (decisionLevel() > level) removeLastAssignment();
}

void Solver::decide(Lit l) {
  ++stats.NDECIDE;
  trail_lim.push_back(trail.size());
  uncheckedEnqueue(l, CRef_Undef);
}

void Solver::propagate(Lit l, CRef reason) {
  assert(reason != CRef_Undef);
  ++stats.NPROP;
  uncheckedEnqueue(l, reason);
}

/**
 * Unit propagation with watched literals.
 * @post: all watches up to trail[qhead] have been propagated
 */
CeSuper Solver::runPropagation(bool onlyUnitPropagation) {
  CeSuper confl = processLearnedStack();
  if (confl) {
    return confl;
  }
  while (qhead < (int)trail.size()) {
    Lit p = trail[qhead++];
    std::vector<Watch>& ws = adj[-p];
    for (int it_ws = 0; it_ws < (int)ws.size(); ++it_ws) {
      int idx = ws[it_ws].idx;
      if (idx < 0 && isTrue(Level, idx + INF)) {
        assert(dynamic_cast<Clause*>(&(ca[ws[it_ws].cref])) != nullptr);
        continue;
      }  // blocked literal check
      CRef cr = ws[it_ws].cref;
      WatchStatus wstat = checkForPropagation(cr, ws[it_ws].idx, -p);
      if (wstat == WatchStatus::DROPWATCH)
        aux::swapErase(ws, it_ws--);
      else if (wstat == WatchStatus::CONFLICTING) {  // clean up current level and stop propagation
        ++stats.NTRAILPOPS;
        for (int i = 0; i <= it_ws; ++i) {
          const Watch& w = ws[i];
          if (w.idx >= INF) ca[w.cref].undoFalsified(w.idx);
        }
        --qhead;
        Constr& C = ca[cr];
        if (!C.isLocked()) {
          cBumpActivity(C);
          recomputeLBD(C);
        }

        stats.NENCFORMULA += C.getOrigin() == Origin::FORMULA;
        stats.NENCLEARNED += C.getOrigin() == Origin::LEARNED;
        stats.NENCBOUND += (C.getOrigin() == Origin::LOWERBOUND || C.getOrigin() == Origin::UPPERBOUND ||
                            C.getOrigin() == Origin::HARDENEDBOUND);
        stats.NENCCOREGUIDED += C.getOrigin() == Origin::COREGUIDED;
        stats.NLPENCGOMORY += C.getOrigin() == Origin::GOMORY;
        stats.NLPENCLEARNEDFARKAS += C.getOrigin() == Origin::LEARNEDFARKAS;
        stats.NLPENCFARKAS += C.getOrigin() == Origin::FARKAS;

        return C.toExpanded(cePools);
      }
    }
  }
  if (onlyUnitPropagation) return CeNull();
  if (lpSolver) {
    std::pair<LpStatus, CeSuper> lpResult =
        aux::timeCall<std::pair<LpStatus, CeSuper>>([&] { return lpSolver->checkFeasibility(); }, stats.LPTOTALTIME);
    assert((lpResult.first == LpStatus::INFEASIBLE) == (lpResult.second && lpResult.second->hasNegativeSlack(Level)));
    return lpResult.second;
  }
  return CeNull();
  ;
}

WatchStatus Solver::checkForPropagation(CRef cr, int& idx, Lit p) {
  assert(isFalse(Level, p));
  Constr& C = ca[cr];
  if (C.isMarkedForDelete()) return WatchStatus::DROPWATCH;
  ++stats.NWATCHLOOKUPS;

  return C.checkForPropagation(cr, idx, p, *this);
}

// ---------------------------------------------------------------------
// Conflict analysis

void Solver::recomputeLBD(Constr& C) {
  if (C.lbd() > 2) {  // constraints with LBD <= 2 won't have score recomputed
    assert(tmpSet.isEmpty());
    for (int i = 0; i < (int)C.size(); i++)
      if (isFalse(Level, C.lit(i))) tmpSet.add(Level[-C.lit(i)]);
    if (C.lbd() > tmpSet.size() + 1) C.setLBD(tmpSet.size());  // simulate Glucose
    tmpSet.clear();
  }
}

CeSuper getAnalysisCE(const CeSuper& conflict, int bitsOverflow, ConstrExpPools& cePools) {
  if (bitsOverflow == 0 || bitsOverflow > conflLimit128) {
    CeArb confl = cePools.takeArb();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit96) {
    Ce128 confl = cePools.take128();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit64) {
    Ce96 confl = cePools.take96();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit32) {
    Ce64 confl = cePools.take64();
    conflict->copyTo(confl);
    return confl;
  } else {
    Ce32 confl = cePools.take32();
    conflict->copyTo(confl);
    return confl;
  }
}

CeSuper Solver::prepareConflictConstraint(CeSuper conflict) {
  conflict->removeUnitsAndZeroes(Level, Pos);
  conflict->saturateAndFixOverflow(getLevel(), (bool)options.weakenFull, options.bitsOverflow.get(),
                                   options.bitsReduced.get(), 0);
  CeSuper confl = getAnalysisCE(conflict, options.bitsOverflow.get(), cePools);
  conflict->reset();
  return confl;
}

void Solver::assignActiveSet(CeSuper confl) {
  assert(actSet.isEmpty());  // will hold the literals that need their activity bumped
  for (Var v : confl->vars) {
    if (options.bumpLits)
      actSet.add(confl->getLit(v));
    else if (!options.bumpOnlyFalse || isFalse(Level, confl->getLit(v)))
      actSet.add(v);
  }
}

Constr& Solver::getReasonConstraint(Lit l) {
  assert(isPropagated(Reason, l));
  Constr& reasonC = ca[Reason[toVar(l)]];
  if (!reasonC.isLocked()) {
    cBumpActivity(reasonC);
    recomputeLBD(reasonC);
  }

  trackReasonConstraintStats(reasonC);

  return reasonC;
}

void Solver::trackReasonConstraintStats(Constr& reasonC) {
  stats.NENCFORMULA += reasonC.getOrigin() == Origin::FORMULA;
  stats.NENCLEARNED += reasonC.getOrigin() == Origin::LEARNED;
  stats.NENCBOUND += (reasonC.getOrigin() == Origin::LOWERBOUND || reasonC.getOrigin() == Origin::UPPERBOUND ||
                      reasonC.getOrigin() == Origin::HARDENEDBOUND);
  stats.NENCCOREGUIDED += reasonC.getOrigin() == Origin::COREGUIDED;
  stats.NLPENCGOMORY += reasonC.getOrigin() == Origin::GOMORY;
  stats.NLPENCLEARNEDFARKAS += reasonC.getOrigin() == Origin::LEARNEDFARKAS;
  stats.NLPENCFARKAS += reasonC.getOrigin() == Origin::FARKAS;
  ++stats.NRESOLVESTEPS;
}

void Solver::bumpLiteralActivity() {
  for (Lit l : actSet.getKeys())
    if (l != 0) vBumpActivity(toVar(l));
  actSet.clear();
}

CeSuper Solver::analyze(CeSuper conflict) {
  assert(conflict->hasNegativeSlack(Level));

  if (logger) logger->logComment("Analyze", stats);
  stats.NADDEDLITERALS += conflict->vars.size();

  CeSuper confl = prepareConflictConstraint(conflict);
  assignActiveSet(confl);

  while (decisionLevel() > 0) {
    if (asynch_interrupt) throw asynchInterrupt;
    Lit l = trail.back();
    if (confl->hasLit(-l)) {
      assert(confl->hasNegativeSlack(Level));

      AssertionStatus status = confl->isAssertingBefore(Level, decisionLevel());
      // Conflict constraint could now be asserting after removing some assignments.
      if (status == AssertionStatus::ASSERTING) break;
      // Constraint is already falsified by before last decision on trail.
      else if (status == AssertionStatus::FALSIFIED) {
        backjumpTo(decisionLevel() - 1);
        continue;
      }

      Constr& reasonC = getReasonConstraint(l);
      reasonC.resolveWith(confl, l, &actSet, *this);
    }
    removeLastAssignment();
  }
  bumpLiteralActivity();

  assert(confl->hasNegativeSlack(Level));
  return confl;
}

std::vector<CeSuper> Solver::extractCore(CeSuper conflict, const IntSet& assumptions, Lit l_assump) {
  if (l_assump != 0) {  // l_assump is an assumption propagated to the opposite value
    assert(assumptions.has(l_assump));
    assert(isFalse(Level, l_assump));
    int pos = Pos[toVar(l_assump)];
    while ((int)trail.size() > pos) removeLastAssignment();
    assert(isUnknown(Pos, l_assump));
    decide(l_assump);
  }

  // Set all assumptions in front of the trail, all propagations later. This makes it easy to do decision learning.
  // For this, we first copy the trail, then backjump to 0, then rebuild the trail.
  // Otherwise, reordering the trail messes up the slacks of the watched constraints (see removeLastAssignment()).
  std::vector<Lit> decisions;  // holds the decisions
  decisions.reserve(decisionLevel());
  std::vector<Lit> props;  // holds the propagations
  props.reserve(trail.size());
  assert(trail_lim.size() > 0);
  for (int i = trail_lim[0]; i < (int)trail.size(); ++i) {
    Lit l = trail[i];
    if (assumptions.has(l) && (isDecided(Reason, l) || !options.cgResolveProp)) {
      decisions.push_back(l);
    } else {
      props.push_back(l);
    }
  }
  backjumpTo(0);

  for (Lit l : decisions) decide(l);
  for (Lit l : props) propagate(l, Reason[toVar(l)]);

  assert(conflict->hasNegativeSlack(Level));
  stats.NADDEDLITERALS += conflict->vars.size();
  conflict->removeUnitsAndZeroes(Level, Pos);
  conflict->saturateAndFixOverflow(getLevel(), (bool)options.weakenFull, options.bitsOverflow.get(),
                                   options.bitsReduced.get(), 0);
  assert(conflict->hasNegativeSlack(Level));
  CeSuper core = getAnalysisCE(conflict, options.bitsOverflow.get(), cePools);
  conflict->reset();

  std::vector<CeSuper> result;
  int resolvesteps = l_assump == 0;  // if l==0, we already had some resolve steps in conflict analysis
  // analyze conflict
  if (core->hasNegativeSlack(assumptions)) {  // early termination core
    result.push_back(core->clone(cePools));
    if (resolvesteps > 0) learnConstraint(result.back(), Origin::LEARNED);
    resolvesteps = 0;
  }
  while (decisionLevel() > 0) {
    if (asynch_interrupt) throw asynchInterrupt;
    if (!options.cgDecisionCore && result.size() > 0) break;
    Lit l = trail.back();
    if (core->hasLit(-l)) {
      if (isDecided(Reason, l)) break;  // no more propagated literals
      Constr& reasonC = ca[Reason[toVar(l)]];
      // TODO: stats? activity?
      reasonC.resolveWith(core, l, nullptr, *this);
      ++resolvesteps;
      if (result.size() == 0 && core->hasNegativeSlack(assumptions)) {  // early termination core
        result.push_back(core->clone(cePools));
        if (resolvesteps > 0) learnConstraint(result.back(), Origin::LEARNED);
        resolvesteps = 0;
      }
    }
    removeLastAssignment();
  }
  if (options.cgDecisionCore && resolvesteps > 0) {  // decision core
    result.push_back(core->clone(cePools));
    learnConstraint(result.back(), Origin::LEARNED);
  }

  // weaken non-falsifieds
  for (CeSuper& cnfl : result) {
    assert(cnfl->hasNegativeSlack(assumptions));
    assert(!cnfl->isTautology());
    assert(cnfl->isSaturated());
    for (Var v : cnfl->vars)
      if (!assumptions.has(-cnfl->getLit(v))) cnfl->weaken(v);
    cnfl->postProcess(Level, Pos, true, stats);
    assert(cnfl->hasNegativeSlack(assumptions));
  }
  backjumpTo(0);

  return result;
}

// ---------------------------------------------------------------------
// Constraint management

CRef Solver::attachConstraint(CeSuper constraint, bool locked) {
  assert(constraint->isSortedInDecreasingCoefOrder());
  assert(constraint->isSaturated());
  assert(constraint->hasNoZeroes());
  assert(constraint->hasNoUnits(getLevel()));
  assert(!constraint->isTautology());
  assert(constraint->vars.size() > 0);
  assert(!constraint->hasNegativeSlack(getLevel()));
  assert(constraint->orig != Origin::UNKNOWN);

  CRef cr = constraint->toConstr(ca, locked, logger ? constraint->logProofLineWithInfo("Attach", stats) : ++crefID);
  Constr& C = ca[cr];
  C.initializeWatches(cr, *this);
  constraints.push_back(cr);

  bool learned = (C.getOrigin() == Origin::LEARNED || C.getOrigin() == Origin::LEARNEDFARKAS ||
                  C.getOrigin() == Origin::FARKAS || C.getOrigin() == Origin::GOMORY);
  if (learned) {
    stats.LEARNEDLENGTHSUM += C.size();
    stats.LEARNEDDEGREESUM += C.degree();
  } else {
    stats.EXTERNLENGTHSUM += C.size();
    stats.EXTERNDEGREESUM += C.degree();
  }
  if (C.degree() == 1) {
    stats.NCLAUSESLEARNED += learned;
    stats.NCLAUSESEXTERN += !learned;
  } else if (C.largestCoef() == 1) {
    stats.NCARDINALITIESLEARNED += learned;
    stats.NCARDINALITIESEXTERN += !learned;
  } else {
    stats.NGENERALSLEARNED += learned;
    stats.NGENERALSEXTERN += !learned;
  }

  stats.NCONSFORMULA += C.getOrigin() == Origin::FORMULA;
  stats.NCONSLEARNED += C.getOrigin() == Origin::LEARNED;
  stats.NCONSBOUND += (C.getOrigin() == Origin::LOWERBOUND || C.getOrigin() == Origin::UPPERBOUND ||
                       C.getOrigin() == Origin::HARDENEDBOUND);
  stats.NCONSCOREGUIDED += C.getOrigin() == Origin::COREGUIDED;
  stats.NLPGOMORYCUTS += C.getOrigin() == Origin::GOMORY;
  stats.NLPLEARNEDFARKAS += C.getOrigin() == Origin::LEARNEDFARKAS;
  stats.NLPFARKAS += C.getOrigin() == Origin::FARKAS;

  return cr;
}

void Solver::learnConstraint(const CeSuper c, Origin orig) {
  assert(orig == Origin::LEARNED || orig == Origin::FARKAS || orig == Origin::LEARNEDFARKAS || orig == Origin::GOMORY);
  CeSuper learned = c->clone(cePools);
  learned->orig = orig;
  learned->saturateAndFixOverflow(getLevel(), (bool)options.weakenFull, options.bitsLearned.get(),
                                  options.bitsLearned.get(), 0);
  learnedStack.push_back(learned->toSimple());
}

// NOTE: backjumps to place where the constraint is not conflicting, as otherwise we might miss propagations
CeSuper Solver::processLearnedStack() {
  // loop back to front as the last constraint in the queue is a result of conflict analysis
  // and we want to first check this constraint to backjump.
  while (learnedStack.size() > 0) {
    CeSuper learned = learnedStack.back()->toExpanded(cePools);
    learnedStack.pop_back();
    learned->removeUnitsAndZeroes(Level, Pos);
    learned->sortInDecreasingCoefOrder();
    int assertionLevel = learned->getAssertionLevel(Level, Pos);
    if (assertionLevel < 0) {
      backjumpTo(0);
      return learned;
    }
    backjumpTo(assertionLevel);
    assert(!learned->hasNegativeSlack(Level));
    learned->heuristicWeakening(Level, Pos, stats);  // TODO: don't always weaken heuristically?
    learned->postProcess(Level, Pos, false, stats);
    assert(learned->isSaturated());
    if (learned->isTautology()) {
      continue;
    }
    CRef cr = attachConstraint(learned, false);
    Constr& C = ca[cr];
    if (assertionLevel < INF)
      recomputeLBD(C);
    else
      C.setLBD(C.size());  // the LBD of non-asserting constraints is undefined, so we take a safe upper bound
  }
  return CeNull();
}

std::pair<ID, ID> Solver::addInputConstraint(CeSuper ce) {
  assert(ce->orig == Origin::FORMULA || ce->orig == Origin::UPPERBOUND || ce->orig == Origin::LOWERBOUND ||
         ce->orig == Origin::HARDENEDBOUND || ce->orig == Origin::COREGUIDED || ce->orig == Origin::BACKBONECALCULATION);
  assert(decisionLevel() == 0);
  ID input = ID_Undef;
  if (logger) {
    if (ce->orig == Origin::FORMULA) {
      input = ce->logInput();
    } else {
      input = ce->logAsAssumption();
    }
  }
  ce->postProcess(Level, Pos, true, stats);
  if (ce->isTautology()) {
    return {input, ID_Undef};  // already satisfied.
  }

  if (ce->hasNegativeSlack(Level)) {
    if (options.verbosity.get() > 0) puts("c Inconsistent input constraint");
    if (logger) ce->logInconsistency(Level, Pos, stats);
    assert(decisionLevel() == 0);
    return {input, ID_Unsat};
  }

  if (options.bitsInput.get() != 0 && !ce->largestCoefFitsIn(options.bitsInput.get())) {
    ce->toStreamAsOPB(std::cerr);
    quit::exit_ERROR({"Input constraint coefficient exceeds bit limit."});
  }

  CRef cr = attachConstraint(ce, true);
  CeSuper confl = aux::timeCall<CeSuper>([&] { return runPropagation(true); }, stats.PROPTIME);
  if (confl) {
    assert(confl->hasNegativeSlack(Level));
    if (options.verbosity.get() > 0) puts("c Input conflict");
    if (logger) confl->logInconsistency(Level, Pos, stats);
    assert(decisionLevel() == 0);
    return {input, ID_Unsat};
  }
  ID id = ca[cr].id;
  Origin orig = ca[cr].getOrigin();
  if (orig != Origin::FORMULA) {
    external[id] = cr;
  }
  if (lpSolver && (orig == Origin::FORMULA || orig == Origin::UPPERBOUND || orig == Origin::LOWERBOUND)) {
    lpSolver->addConstraint(cr, false, orig == Origin::UPPERBOUND, orig == Origin::LOWERBOUND);
  }
  return {input, id};
}

std::pair<ID, ID> Solver::addConstraint(const CeSuper c, Origin orig) {
  // NOTE: copy to temporary constraint guarantees original constraint is not changed and does not need logger
  CeSuper ce = c->clone(cePools);
  ce->orig = orig;
  std::pair<ID, ID> result = addInputConstraint(ce);
  return result;
}

std::pair<ID, ID> Solver::addConstraint(const ConstrSimpleSuper& c, Origin orig) {
  CeSuper ce = c.toExpanded(cePools);
  ce->orig = orig;
  std::pair<ID, ID> result = addInputConstraint(ce);
  return result;
}

std::pair<ID, ID> Solver::addUnitConstraint(Lit l, Origin orig) {
  return addConstraint(ConstrSimple32({{1, l}}, 1), orig);
}

void Solver::removeConstraint(Constr& C, [[maybe_unused]] bool override) {
  assert(override || !C.isLocked());
  assert(!C.isMarkedForDelete());
  assert(!external.count(C.id));
  C.markForDel();
  ca.wasted += C.getMemSize();
}

void Solver::dropExternal(ID id, bool erasable, bool forceDelete) {
  assert(erasable || !forceDelete);
  if (id == ID_Undef) return;
  auto old_it = external.find(id);
  assert(old_it != external.end());
  Constr& constr = ca[old_it->second];
  external.erase(old_it);
  constr.setLocked(!erasable);
  if (forceDelete) removeConstraint(constr);
}

// ---------------------------------------------------------------------
// Assumptions

void Solver::setAssumptions(const IntSet& assumps) {
  assumptions = assumps;
  backjumpTo(0);
}

// ---------------------------------------------------------------------
// Garbage collection

// We assume in the garbage collection method that reduceDB() is the
// only place where constraints are deleted.
void Solver::garbage_collect() {
  // if (options.verbosity.get() > 0) puts("c GARBAGE COLLECT");

  ca.wasted = 0;
  ca.at = 0;
  std::unordered_map<uint32_t, CRef> crefmap;
  for (int i = 1; i < (int)constraints.size(); ++i) assert(constraints[i - 1].ofs < constraints[i].ofs);
  for (CRef& cr : constraints) {
    uint32_t offset = cr.ofs;
    size_t memSize = ca[cr].getMemSize();
    memmove(ca.memory + ca.at, ca.memory + cr.ofs, sizeof(uint32_t) * memSize);
    cr.ofs = ca.at;
    ca.at += memSize;
    crefmap[offset] = cr;
  }
#define update_ptr(cr) cr = crefmap[cr.ofs];
  for (Lit l = -n; l <= n; ++l) {
    for (size_t i = 0; i < adj[l].size(); i++) update_ptr(adj[l][i].cref);
  }
  for (Var v = 1; v <= n; ++v) {
    if (Reason[v] != CRef_Undef) update_ptr(Reason[v]);
  }
  for (auto& ext : external) update_ptr(ext.second);
#undef update_ptr
}

// We assume in the garbage collection method that reduceDB() is the
// only place where constraints are removed from memory.
void Solver::reduceDB() {
  std::vector<CRef> learnts;
  learnts.reserve(constraints.size() / 2);

  size_t totalLearnts = 0;
  size_t promisingLearnts = 0;
  for (CRef& cr : constraints) {
    Constr& C = ca[cr];
    if (C.isMarkedForDelete() || external.count(C.id)) continue;
    bool isReason = false;
    for (unsigned int i = 0; i < C.size() && !isReason; ++i) {
      isReason = Reason[toVar(C.lit(i))] == cr;
    }
    if (isReason) continue;
    if (C.isSatisfiedAtRoot(Level))
      removeConstraint(C, true);
    else if (!options.keepAll && !C.isLocked()) {
      if (C.size() > 2 && C.lbd() > 2) learnts.push_back(cr);  // Keep all binary clauses and short LBDs
      if (C.size() <= 2 || C.lbd() <= 3) ++promisingLearnts;
      ++totalLearnts;
    }
  }

  if (promisingLearnts > totalLearnts / 2)
    nconfl_to_reduce += 10 * options.dbCleanInc.get();
  else
    nconfl_to_reduce += options.dbCleanInc.get();
  std::sort(learnts.begin(), learnts.end(), [&](CRef x, CRef y) {
    return ca[x].lbd() > ca[y].lbd() || (ca[x].lbd() == ca[y].lbd() && ca[x].act < ca[y].act);
  });
  for (size_t i = 0; i < std::min(totalLearnts / 2, learnts.size()); ++i) removeConstraint(ca[learnts[i]]);

  for (Lit l = -n; l <= n; ++l)
    for (int i = 0; i < (int)adj[l].size(); ++i) {
      if (ca[adj[l][i].cref].isMarkedForDelete()) aux::swapErase(adj[l], i--);
    }

  size_t i = 0;
  size_t j = 0;
  for (; i < constraints.size(); ++i) {
    Constr& c = ca[constraints[i]];
    if (c.isMarkedForDelete()) {
      c.freeUp();  // free up indirectly owned memory before implicitly deleting c during garbage collect
    } else {
      constraints[j++] = constraints[i];
    }
  }
  constraints.resize(j);
  if ((double)ca.wasted / (double)ca.at > 0.2) garbage_collect();
}

// ---------------------------------------------------------------------
// Solving

double Solver::luby(double y, int i) const {
  // Find the finite subsequence that contains index 'i', and the
  // size of that subsequence:
  int size, seq;
  for (size = 1, seq = 0; size < i + 1; seq++, size = 2 * size + 1) {
  }
  while (size != i + 1) {
    size = (size - 1) >> 1;
    --seq;
    assert(size != 0);
    i = i % size;
  }
  return std::pow(y, seq);
}

bool Solver::checkSAT() {
  for (CRef cr : constraints) {
    const Constr& C = ca[cr];
    if (C.getOrigin() == Origin::FORMULA && !C.toExpanded(cePools)->isSatisfied(getLevel())) return false;
  }
  return true;
}

Lit Solver::pickBranchLit(bool lastSolPhase) {
  Var next = 0;
  // Activity based decision:
  while (next == 0 || !isUnknown(Pos, next)) {
    if (order_heap.empty())
      return 0;
    else
      next = order_heap.removeMax();
  }
  assert(phase[0] == 0);
  assert(lastSol[0] == 0);
  return (lastSolPhase && (int)lastSol.size() > next) ? lastSol[next] : phase[next];
}

// function to add the constraint of the negated literals
std::pair<ID,ID> Solver::activateBackboneConstraint(std::vector<int>& chunk, IntSet& assumps) {
  
  // create a new constraint object
  CeArb input = cePools.takeArb();
  input->reset();

  // add the literals to the left side of the constraint
  for (Lit l : chunk) {
    input->addLhs(1, -l);
  }

  // we are checking if at least one literal is not backbone
  input->addRhs(1); 

  // add auxiliary variables
  long long oldN = getNbVars();
  long long newN = oldN + 1;
  setNbVars(newN);

  // add the left side of the cardinality constraint
  input->addLhs(1, newN);

  // assume the negation of the aux variable
  assumps.add(-newN);

  // add the constraint to the solver
  std::pair<ID, ID> cons = addConstraint(input, Origin::BACKBONECALCULATION);

  // check the status of the constraint (it should never happen to be unsat)
  assert(cons.second != ID_Unsat);

  return {cons.second, newN};
}

bool isStabilizing(double currentValue, double threshold, int& stableCount, int stabilityPeriod) {
    if (currentValue <= threshold) {
        stableCount++;
    } else {
        stableCount = 0; // Reset stable count if the difference exceeds the threshold
    }
    return stableCount >= stabilityPeriod;
}

void Solver::cleanVarState(const std::vector<Var> ordered, std::vector<Lit> polarities, SolveState reply, 
  int& stableCount, bool& flippableActive) {

  int solSize = (int)lastSol.size() - stats.NAUXVARS;
  double flippableTimeStart;
  double flippableTimeEnd;

  std::vector<Lit> flippableSoFar(solSize, true); // flippable means that is not a candidate
  int nModelDiff = 0; // number of variables that changed from the last model

  // metrics for checking if it's necessary to call the flippable code
  int threshold = int(0.0005 * solSize) < 5 ? 5 : int(0.0005 * solSize); // threshold for recognizing stability
  int stabilityPeriod = 10; // number of consecutive iterations needed for stability
  int flippableCounter = 0; // counter for the number of flippable variables recognized in the last iteration

  for (int i = 1; i < ordered.size(); i++) {
    Var var = ordered[i];
    Lit lit = polarities[var] * var;

    if (state[var] != 0) continue;

    // check if it's fixed
    if (isUnit(Level, lit)) {
      state[var] = 1;
      backboneCounter++;
      continue;
    }

    // if the polarity changed from last model, it's not a backbone
    if (lastSol[var] * polarities[var] < 0) {
      state[var] = -1;
      noBackboneCounter++;
      nModelDiff++;
      continue;
    }
  }
  
  if (reply != SolveState::SAT) return;

  // for example, if the flippable code doesn't change the state of the variables in many iterations, 
  // reduce the calls to it (adjust it based on the utility)
  if (!flippableActive) return;
  else if (isStabilizing(flippableCounter, threshold, stableCount, stabilityPeriod)) {

    // as it's stable there's no need to call the flippable code
    stableCount = 0;
    return;
  }

  flippableTimeStart = rs::aux::cpuTime();

  // iterate over constraints
  for (CRef cr : constraints) {
    const Constr& C = ca[cr];
    if(!(C.getOrigin() == Origin::FORMULA)) break;

    CeSuper ce = C.toExpanded(cePools);
    bigint slack = ce->calculateSlack(std::cout, lastSol);

    // iterate over literals in the constraint
    for (int i = 0; i < C.size(); i++) {

      Lit lit = C.lit(i);
      Var var = toVar(lit);
      
      // if the literal cannot be modified without making this constraint false, then it's fixed (not flippable)
      if (lit == lastSol[var] && state[var] == 0 && flippableSoFar[var] && C.coef(i) > slack) {
        flippableSoFar[var] = false;
      }
    }
  }

  for (int i = 1; i < ordered.size(); i++) {
    Var var = ordered[i];
    if (state[var] == 0 && flippableSoFar[var]) {
      state[var] = -1;
      noBackboneCounter++;
      flippableCounter++;
    }
  }

  flippableTimeEnd = rs::aux::cpuTime();
  stats.backboneMetrics.back().flippableTime = flippableTimeEnd - flippableTimeStart;
}


void Solver::calculateBackbones() {

  std::cout << "=================== calculateBackbones was called ===================" << std::endl;
  std::cout << "Instance: " << options.formulaName << std::endl;

  stats.isCalculatingBackbones = true;
  int solSize = (int)lastSol.size() - stats.NAUXVARS;

  // initialize stats
  stats.STARTTIMEBACKBONE = rs::aux::cpuTime();
  stats.startTimeBackbone = std::chrono::high_resolution_clock::now();
  stats.NITBACKBONES = 0;
  stats.backboneMetrics.reserve(int(solSize * 0.05)); // initialized on a heuristically decided amount, to reduce resizing cost
  
  SolveState reply = SolveState::UNSAT; // reply of each call to the solver
  IntSet assumps; // assumps for each iteration
  std::pair<ID,ID> cons; // constraint added for checking backbone property on multiple variables
  int k = 1; // chunk size
  int j = 0;

  // counters
  backboneCounter = 0;
  noBackboneCounter = 0;

  // metrics
  int stableCount = 0; // consecutive stable iterations
  bool flippableActive = true; // flag to check if the flippable code is active

  double startTimeIter; //
  double endTimeIter; // keep track of the time used in each iter, for profiling purposes
  double startTimeProp; // 
  double endTimeProp; // keep track of the propagation time for each time we discover a new unit

  std::vector<int> firstSol = lastSol; // copy of sol to compare with the new one
  std::vector<Lit> chunk = {}; // contains the positions of the literals in the chunk
  std::vector<Lit> polarities(solSize, 0); // contains the polarity of the literals

  // initialize the ordered vector
  stats.STARTTIMEORDER = rs::aux::cpuTime();
  std::vector<Var> ordered(lastSol.begin(), lastSol.begin() + solSize); // contains the literals ordered lexicographically
  stats.ENDTIMEORDER = rs::aux::cpuTime();
  
  // save initial metrics: ORDERING TIME, NUMBER OF VARIABLES, NUMBER OF UNITS
  double orderingTime = stats.ENDTIMEORDER - stats.STARTTIMEORDER;
  int nUnits = 0;

  // count units in first solution
  for (int i = 1; i < ordered.size(); i++) {
    Var var = ordered[i];
    if (isUnit(Level, var)) {
      nUnits++;
    }
  }

  // stats.backboneMetrics.push_back({orderingTime, solSize, nUnits, 0, 0.0, 0.0});
  stats.backboneMetrics[0].time = orderingTime;
  stats.backboneMetrics[0].backbones = nUnits;
  stats.backboneMetrics[0].noBackbones = solSize-1;

  std::cout << "Ordering time: " << orderingTime << ", Sol size: " << solSize-1 << ", Units 1st sol: " <<
   nUnits << ", Time 1st sol:" << stats.getTime() << std::endl;

  // initializate polarity
  for (int i = 1; i < polarities.size(); i++) {
    polarities[i] = firstSol[i] > 0 ? 1 : -1;
    ordered[i] = toVar(ordered[i]);
    assert(ordered[i] != 0);
  }

  // print the ordered vector and first solution to a file
  stats.printOrderVariables(ordered, options.formulaName);
  stats.printFirstSolution(options.formulaName, firstSol);

  // change the state of the found variables
  cleanVarState(ordered, polarities, SolveState::SAT, stableCount, flippableActive);

  // iterate over ordered vector
  for (int i = 1; i < ordered.size(); i++) {

    stats.NITBACKBONES++;
    startTimeIter = rs::aux::cpuTime();

    // check if the timeout was reached or if the user interrupted the execution
    if (asynch_interrupt or stats.TLE_backbone_computing()) break;
    if (backboneCounter + noBackboneCounter == solSize - 1) break;

    stats.backboneMetrics.push_back({0.0, 0, 0, 0, 0.0, 0.0});

    // get chunk
    chunk.clear();
    j = i;
    k = std::min(k, solSize - backboneCounter - noBackboneCounter - 1);
    assert(k > 0);

    // at this point the state vector is updated
    while (chunk.size() < k && j < ordered.size()) {

      assert(ordered[j] != 0); // check if there's a malformation in ordered vector (only 1st element should be 0)

      if (state[ordered[j]]) {
        j++;
        continue;
      }

      Var var = ordered[j];
      assert(var > 0 || polarities[var] == 0);
      Lit lit = var * polarities[var];
      chunk.push_back(lit);
      j++;
    }

    // reduce the size of j in one
    j--;
    
    assert(chunk.size() == k);

    // create constraint with the negated literals in the chunk
    if (chunk.size() == 1) assumps.add(-chunk[0]);
    else cons = activateBackboneConstraint(chunk, assumps);

    // solve the problem again with that assumptions
    reply = SolveState::UNSAT;

    while (reply != SolveState::SAT && reply != SolveState::INCONSISTENT) {
      setAssumptions(assumps);
      SolveAnswer out = solve();
      reply = out.state;
    }
    
    assumps.clear();

    // if the solver didn't learn the unit clause
    if (chunk.size() != 1 && !isUnit(Level, cons.second)) {

      // change decision level to 0 and add the aux var as unit
      backjumpTo(0);
      uncheckedEnqueue(cons.second, CRef_Undef);

      // at this point the aux variable must be a unit
      assert(isUnit(Level, cons.second));
    }

    assert(reply == SolveState::SAT || reply == SolveState::INCONSISTENT);

    if(reply == SolveState::SAT) {
      
      // TODO Matias: if it's SAT store in the metrics vector a 1 
      stats.backboneMetrics[stats.NITBACKBONES].sat = 1;

      // if there's only one literal in the chunk
      if (chunk.size() == 1) {
        state[toVar(chunk[0])] = -1;
        noBackboneCounter++;
      }

      else i--; // dirty hack, could be improved

      // set the k value to 1 and repeat this iteration
      k = 1;
    }

    else if(reply == SolveState::INCONSISTENT) {

      // add all the literals as unit constraints
      for (Lit chunkLit : chunk) {
        
        // create unit clause
        if (!isUnit(Level, chunkLit)) {
          
          // change decision level to 0 and add the literal to the trail
          backjumpTo(0);
          uncheckedEnqueue(chunkLit, CRef_Undef);

          // propagate the unit clause to prevent input conflicts
          // there should be no conflicts at this point
          startTimeProp = rs::aux::cpuTime();
          CeSuper confl = aux::timeCall<CeSuper>([&] { return runPropagation(false); }, stats.PROPTIME);
          endTimeProp = rs::aux::cpuTime();
        }

        assert(isUnit(Level, chunkLit));
      }

      // skip the literals that we know that are not unknown
      i = j;
      k *= 2;
    }
    
    cleanVarState(ordered, polarities, reply, stableCount, flippableActive);

    endTimeIter = rs::aux::cpuTime();
    stats.backboneMetrics[stats.NITBACKBONES].time = endTimeIter - startTimeIter;
    stats.backboneMetrics[stats.NITBACKBONES].propTime = endTimeProp - startTimeProp;
    stats.backboneMetrics[stats.NITBACKBONES].backbones = backboneCounter;
    stats.backboneMetrics[stats.NITBACKBONES].noBackbones = noBackboneCounter;
  
    std::cout << stats.NITBACKBONES << ": " << stats.backboneMetrics[stats.NITBACKBONES].time << ", " << 
      stats.backboneMetrics[stats.NITBACKBONES].noBackbones << ", " << stats.backboneMetrics[stats.NITBACKBONES].backbones << 
      ", " << stats.backboneMetrics[stats.NITBACKBONES].sat << ", " << stats.backboneMetrics[stats.NITBACKBONES].propTime << 
      ", " << stats.backboneMetrics[stats.NITBACKBONES].flippableTime << std::endl;
  }

  // print the metrics to a file
  double printMetricsTimeStart = rs::aux::cpuTime();
  stats.printMetricsBackbones(options.formulaName);
  // std::cout << "Time metrics print: " << rs::aux::cpuTime() - printMetricsTimeStart << "s" << std::endl;

  // print the backbone set to a file
  stats.printBackbones(state, firstSol, options.formulaName);

  std::cout << "s CALCULATED BACKBONE\n";
}

void Solver::presolve() {
  firstRun = false;
  if (lpSolver) aux::timeCall<void>([&] { lpSolver->inProcess(); }, stats.LPTOTALTIME);
}

SolveAnswer Solver::solve() {
  if (firstRun) presolve();
  std::vector<int> assumptions_lim = {0};
  assumptions_lim.reserve((int)assumptions.size() + 1);
  bool runLP = false;
  while (true) {
    if (asynch_interrupt || (options.timeout.get() != -1.0 && stats.getTime() > options.timeout.get())) {
      if (stats.isCalculatingBackbones) stats.printIncompleteBackbones(state, lastSol, options.formulaName);
      throw asynchInterrupt;
    }
    CeSuper confl = aux::timeCall<CeSuper>([&] { return runPropagation(runLP); }, stats.PROPTIME);
    runLP = !confl;
    if (confl) {
      assert(confl->hasNegativeSlack(Level));
      vDecayActivity();
      cDecayActivity();
      stats.NCONFL++;
      nconfl_to_restart--;
      if (stats.NCONFL % 1000 == 0 && options.verbosity.get() > 0) {
        // printf("c #Conflicts: %10lld | #Constraints: %10lld\n", stats.NCONFL, (long long)constraints.size());
        if (options.verbosity.get() > 2) {
          // memory usage
          std::cout << "c total constraint space: " << ca.cap * 4 / 1024. / 1024. << "MB" << std::endl;
          std::cout << "c total #watches: ";
          long long cnt = 0;
          for (Lit l = -n; l <= n; l++) cnt += (long long)adj[l].size();
          std::cout << cnt << std::endl;
        }
      }
      if (decisionLevel() == 0) {
        if (logger) {
          confl->logInconsistency(Level, Pos, stats);
        }
        return {SolveState::UNSAT, {}, lastSol};
      } else if (decisionLevel() >= (int)assumptions_lim.size()) {
        CeSuper analyzed = aux::timeCall<CeSuper>([&] { return analyze(confl); }, stats.CATIME);
        assert(analyzed);
        assert(analyzed->hasNegativeSlack(getLevel()));
        assert(analyzed->isSaturated());
        if (learnedStack.size() > 0 && learnedStack.back()->orig == Origin::FARKAS)
          learnConstraint(analyzed, Origin::LEARNEDFARKAS);  // TODO: ugly hack
        else
          learnConstraint(analyzed, Origin::LEARNED);
      } else {
        std::vector<CeSuper> result =
            aux::timeCall<std::vector<CeSuper>>([&] { return extractCore(confl, assumptions); }, stats.CATIME);
        for ([[maybe_unused]] const CeSuper& core : result) {
          assert(core);
          assert(core->hasNegativeSlack(assumptions));
        }
        return {SolveState::INCONSISTENT, result, lastSol};
      }
    } else {  // no conflict
      if (nconfl_to_restart <= 0) {
        backjumpTo(0);
        double rest_base = luby(options.lubyBase.get(), ++stats.NRESTARTS);
        nconfl_to_restart = (long long)rest_base * options.lubyMult.get();
        //        return {SolveState::RESTARTED, {}, lastSol}; // avoid this overhead for now
      }
      if (stats.NCONFL >= (stats.NCLEANUP + 1) * nconfl_to_reduce) {
        // if (options.verbosity.get() > 0) puts("c INPROCESSING");
        ++stats.NCLEANUP;
        reduceDB();
        while (stats.NCONFL >= stats.NCLEANUP * nconfl_to_reduce) nconfl_to_reduce += options.dbCleanInc.get();
        if (lpSolver) aux::timeCall<void>([&] { lpSolver->inProcess(); }, stats.LPTOTALTIME);
        return {SolveState::INPROCESSED, {}, lastSol};
      }
      Lit next = 0;
      if ((int)assumptions_lim.size() > decisionLevel() + 1) assumptions_lim.resize(decisionLevel() + 1);
      if (assumptions_lim.back() < (int)assumptions.size()) {
        for (int i = (decisionLevel() == 0 ? 0 : trail_lim.back()); i < (int)trail.size(); ++i) {
          Lit l = trail[i];
          if (assumptions.has(-l)) {  // found conflicting assumption
            if (isUnit(Level, l)) {   // negated assumption is unit
              backjumpTo(0);
              return {SolveState::INCONSISTENT, {cePools.take32()}, lastSol};
            } else {
              std::vector<CeSuper> result = aux::timeCall<std::vector<CeSuper>>(
                  [&] { return extractCore(ca[Reason[toVar(l)]].toExpanded(cePools), assumptions, -l); }, stats.CATIME);
              for ([[maybe_unused]] const CeSuper& core : result) {
                assert(core);
                assert(core->hasNegativeSlack(assumptions));
              }
              return {SolveState::INCONSISTENT, result, lastSol};
            }
          }
        }
      }
      while (assumptions_lim.back() < (int)assumptions.size()) {
        assert(decisionLevel() == (int)assumptions_lim.size() - 1);
        Lit l_assump = assumptions.getKeys()[assumptions_lim.back()];
        assert(!isFalse(Level, l_assump));  // otherwise above check should have caught this
        if (isTrue(Level, l_assump)) {      // assumption already propagated
          ++assumptions_lim.back();
        } else {  // unassigned assumption
          next = l_assump;
          assumptions_lim.push_back(assumptions_lim.back() + 1);
          break;
        }
      }
      if (next == 0) next = pickBranchLit(assumptions.isEmpty() && options.cgSolutionPhase);
      if (next == 0) {
        assert(order_heap.empty());
        assert((int)trail.size() == getNbVars());
        assert(checkSAT());
        lastSol.clear();
        // We want to keep track of the full solution for e.g. branching heuristics based on last found solution
        lastSol.resize(getNbVars() + 1);
        lastSol[0] = 0;
        for (Var v = 1; v <= getNbVars(); ++v) lastSol[v] = isTrue(Level, v) ? v : -v;
        backjumpTo(0);
        return {SolveState::SAT, {}, lastSol};
      }
      decide(next);
    }
  }
}

}  // namespace rs
