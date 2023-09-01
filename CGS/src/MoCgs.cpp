//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the MoCgs class.
////////////////////////////////////////////////////////////////////////////////

// INCLUDES
// #############################################################################
#include "MoCgs.hpp"

#include "TrainingPoint.hpp"
#include "utilities.hpp"

#include "Evaluators\include\Evaluator.hpp"

#include <algorithm>
#include <chrono>
#include <cassert>
#include <fstream>
#include <iostream>     // temp
#include <limits>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using bayesnet::TrainingPoint;

using std::chrono::system_clock;
using std::chrono::duration;
using std::mt19937;
using std::numeric_limits;
using std::string;
using std::unordered_set;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiate a Multi-Objective CGS solver object.
///
/// \param  evaluator       A multi-objective optimization problem evaluator.
/// \param  domains         The set of discrete variable domains of the
///                         design variables.
/// \param  rngSeed         A seed value for the random number generator.
/// \param  graphFilename   The filename of the initial classifier Bayesian
 ///                        network edge list.
////////////////////////////////////////////////////////////////////////////////
MoCgs::MoCgs(
    evaluators::Evaluator& evaluator,
    const vector<bayesnet::Domain>& domains,
    const mt19937::result_type& rngSeed,
    const string& graphFilename
    )
    :   Solver(evaluator, domains, rngSeed, graphFilename)
{
}   // MoCgs


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Main CGS algorithms multi-objective solver function
///
/// \param  nBest           Number of designs in the "best" or "good" class.
/// \param  nInitRandom     Number of random designs to generate at the
///                         beginning of a solution run (i.e., before solutions
///                         begin being generated with the guidance of the
///                         classifier distributions).
/// \param  batchSize       The number of solutions sampled and
///                         evaluated on each iteration between classifier
///                         updates.
/// \param  parentLimit     The upper limit on the number of
///                         parents a node is allowed to have in the 
///                         Bayesian network. If set to 0, no automatic
///                         network learning will occur.
////////////////////////////////////////////////////////////////////////////////
void
MoCgs::solve(int nBest, int nInitRandom, int batchSize, int parentLimit)
{
    // initializations
    unordered_set<MultiResult, KnownPoint::Hash> finalResults;

    long long evalCount = 0;
    long iterationCount = 0;
    system_clock::time_point startTime(system_clock::now());
    double elapsedTime(0.0);

    DesignSet newDesigns;


    // set classifier priors as a fixed uniform discrete distribution
    this->fixClassifierPriors();

    int iterationNo(0);

    // begin solution process
    while (!this->anyConvergersSatisfied(
            evalCount, iterationCount, elapsedTime
            )) {
        
        // clear the set of designs sampled for evaluation this iteration
        newDesigns.clear();

        // sample and screen a collection of new points to evaluate
        while (newDesigns.size() != batchSize) {

            // sample a new design
            vector<int> xNew(
                this->sampleNewX(nInitRandom, evalCount, newDesigns)
                );

            // accept or reject based on classifier posterior probability
            //if (true) {
            if (this->isAcceptedForEval(xNew)) {
                
                // add it to the set to be evaluated
                newDesigns.insert(xNew);
            }
        }

        // evaluate the new points
        DesignSet::iterator i;
        for (i = newDesigns.begin(); i != newDesigns.end(); ++i) {
            
            // evaluate the design
            vector<double> evalResult(this->evaluateX(*i));
            ++evalCount;
            //std::cout << evalCount << std::endl;
            int m = this->getEvaluator().getNObjectives();

            vector<double> fNew(evalResult.begin(), evalResult.begin() + m);
            vector<double> gNew(evalResult.begin() + m, evalResult.end());

            // add the new result to the set of best desings
            bestSet_.push_back(MultiResult(TrainingPoint(*i, 1), gNew, fNew));

            // update the violation stats
            if (!bestSet_.rbegin()->isFeasible()) {
                this->updateViolationStats(*(bestSet_.rbegin()));
            }
        }

        // update the violation metrics for all solutions in set of best
        this->updateBestSolutionViolations(bestSet_);

        // compute the domination rank of all solutions in current set of best
        this->updateDominationRank();

        // sort the set of best designs by domination rank (preferring feasible
        // to infeasible, less infeasible to more infeasible)
        std::sort(
            bestSet_.begin(), 
            bestSet_.end(), 
            MultiResult::DominationRankCompare()
            );
        
        // sort the "last front" via the crowding distance (i.e., the last front
        // whose members can't be fully accommodated in set of the top nBest 
        // solutions. only need to do this if the size of bestSet_ is greater 
        // than nBest and if the nBest'th element is feasible. if it's 
        // infeasible, the elements are already correctly sorted (i.e., sorted 
        // by amount of infeasibility)
        if (bestSet_.size() > nBest && bestSet_[nBest - 1].isFeasible()) {
            
            this->sortLastFrontByCrowdingDistance(nBest);
        }

        // update the classifier
        this->updateClassifier(bestSet_, nBest);

        // merge all non-dominated points in current best set with the master
        // set of all non-dominated designs
        this->updateFinalResultSet(finalResults);

        // trim any additional points off of set of best results
        if (bestSet_.size() > nBest) {

            bestSet_.erase(bestSet_.begin() + nBest, bestSet_.end());
        }

        // update tracked objectives
        // TODO: Implement this

        // if network learning is enabled, tell the classifier to learn the
        // network based on the latest training data and update the classifier
        if (parentLimit > 0) {
            
            this->learnBayesianNetwork(parentLimit);
        }

        // assert that the size of bestSet_ is equal to nBest
        #ifdef _DEBUG
            assert(bestSet_.size() == nBest);
        #endif

        // For Testing: save the best set objective values to a .csv file (each
        // iteration)
        iterationNo += 1;
        string filename("xParetoResults" + std::to_string(iterationNo) + ".csv");
        std::ofstream ofs(filename);
        unordered_set<MultiResult, KnownPoint::Hash>::const_iterator j;
        for (j = finalResults.begin(); j != finalResults.end(); ++j) {
            utilities::readwrite::writeVecToCsv(ofs, j->getObjectives());
        }

        // increment the number of iterations
        iterationCount += 1;

        // update elapsed time
        elapsedTime = duration<double>(system_clock::now() - startTime).count();
    }

    // For Testing: save the best set objective values to a .csv file (final
    // best)
    //string filename("ParetoResults.csv");
    //std::ofstream ofs(filename);
    //unordered_set<MultiResult, KnownPoint::Hash>::const_iterator i;
    //for (i = finalResults.begin(); i != finalResults.end(); ++i) {
    //    utilities::readwrite::writeVecToCsv(ofs, i->getObjectives());
    //}

}   // solve


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if any of the convergence criteria are satisfied. Returns
///         true when any one is satisfied.  I.e., the solver algorithm
///         terminates when the first converger is satisfied.
///
/// \param  evalCount       Current number of objective function evaluations
/// \param  iterationCount  Current number of CGS iterations.
/// \param  elapsedTime     Current solution time in seconds.
///
/// \return     Whether any stopping criterion is satisfied.
////////////////////////////////////////////////////////////////////////////////
bool
MoCgs::anyConvergersSatisfied(
    long long evalCount, long iterationCount, double elapsedTime
    )
{
    // MAXIMUM EVALUATIONS STOPPING CRITERION
    if (this->maxEvaluationsSatisfied(evalCount)) {

        return true;
    }

    // MAXIMUM ITERATIONS STOPPING CRITERION
    if (this->maxIterationsSatisfied(iterationCount)) {

        return true;
    }

    // MAXIMUM SOLUTION TIME STOPPING CRITERION
    if (this->maxTimeSatisfied(elapsedTime)) {

        return true;
    }

    // return false if none returned true
    return false;

}   // anyConvergersSatisfied


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the domination ranking of all feasible points in the current
///         set of best solutions. Uses the NSGA-II non-dominated sorting
///         approach.
///
/// \param  ignoreInfeasible    Whether to ignore infeasible designs
///                             when computing the domination rankings.
////////////////////////////////////////////////////////////////////////////////
void
MoCgs::updateDominationRank(bool ignoreInfeasible)
{   
    // in this function, the indices p and q correspond to the position of each
    // result in the current vector of best results
    typedef vector<MultiResult>::size_type Index;

    vector< vector<Index> > dominates(bestSet_.size());
    vector<long> nDominatedBy(bestSet_.size());
    vector< vector<Index> > fronts(1);

    // for each result in the current set of best results, figure out which
    // designs it dominates and how many designs it is dominated by.

    for (Index p = 0; p != bestSet_.size(); ++p) {
        
        if (ignoreInfeasible && !bestSet_[p].isFeasible()) {
            
            continue;
        }

        for (Index q = 0; q != bestSet_.size(); ++q) {
            
            if ((ignoreInfeasible && !bestSet_[q].isFeasible()) || p == q ) {

                continue;
            }

            if (isDominatedBy(bestSet_[q], bestSet_[p])) {
                
                // p dominates q
                dominates[p].push_back(q);

            } else if (isDominatedBy(bestSet_[p], bestSet_[q])) {
                
                // p is dominated by q
                nDominatedBy[p] += 1;
            }
        }

        // if result i is has domination count = 0, it is in the outer front
        if (nDominatedBy[p] == 0) {
            bestSet_[p].setDominationRank(0);
            fronts[0].push_back(p);
        }
    }


    // find the domination rank of each result.

    // initialize the front counter
    long i(0);

    while (fronts[i].size() != 0) {

        fronts.push_back(vector<Index>());

        for (Index p : fronts[i]) {
            
            if (ignoreInfeasible && !bestSet_[p].isFeasible()) {
                
                continue;
            }
            
            for (Index q : dominates[p]) {
                
                if ((ignoreInfeasible && !bestSet_[q].isFeasible()) || p == q) {

                    continue;
                }
                
                nDominatedBy[q] -= 1;
                
                if (nDominatedBy[q] == 0) {

                    bestSet_[q].setDominationRank(i + 1);
                    fronts[i + 1].push_back(q);
                }
            }
        }
        ++i;
    }

}   // updateDominationRank


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sorts the last front that can't be fully accommodated in the set
///         of the top nBest solutions. 
///
///         Sorting is done by crowding distance as defined by NSGA-II.
///         This function should ONLY be called when the following are
///         true of the current set of best designs:
///         1) it is sorted by domination rank and feasibility,
///         2) the nBest'th in it solution is feasible, and
///         3) there are more than nBest designs in it.
///
/// \param  nBest       The number of designs in the set of best solutions at
///                     (i.e., the number of designs that are assigned the
///                     class label of "good" for classifier training at each
///                     iteration).
////////////////////////////////////////////////////////////////////////////////
void
MoCgs::sortLastFrontByCrowdingDistance(int nBest)
{
    // get indices of the first and (one past) the last design that is a 
    // member of the last front
    int firstIndex = nBest - 1;
    int lastIndex = nBest - 1;
    int lastFrontRank = bestSet_[nBest - 1].getDominationRank();

    while (firstIndex != 0 && 
            bestSet_[firstIndex - 1].getDominationRank() == lastFrontRank) {

        firstIndex -= 1;
    }

    while (lastIndex != static_cast<int>(bestSet_.size()) &&
            bestSet_[lastIndex].getDominationRank() == lastFrontRank) {

        lastIndex += 1;
    }
    
    // if the last index did not advance, there are no designs in the last front
    // that will not fit in the set of nBest designs. if this is the case, there
    // is no reason to compute and sort by crowding distances
    if (lastIndex >= nBest) {
        
        // convert indices to iterators
        vector<MultiResult>::iterator first(bestSet_.begin() + firstIndex);
        vector<MultiResult>::iterator last(bestSet_.begin() + lastIndex);

        // compute the crowding distances for each member in range (first, last]
        this->computeCrowdingDistances(first, last);

        // sort the members by crowding distance
        sort(first, last, MultiResult::CrowdingCompare());
    }

}   // sortLastFrontByCrowdingDistance


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the crowding distances of the last front that can't be 
///         fully accommodated in the set of the top nBest solutions. 
///
///         Crowding distance as is computed as defined by NSGA-II.
///         Note that this function should ONLY be called when the following are
///         true of the current set of best designs:
///         1) it is sorted by domination rank and feasibility,
///         2) when the nBest'th solution is feasible, and
///         3) there are more that nBest designs in it.
///
/// \param  first   An iterator to the first design in bestSet_ that is a member
///                 of the domination rank set that can't be accommodated in 
///                 the top nBest designs.
/// \param  last    An iterator to (one past) the last design in bestSet_ that 
///                 is a member of the domination rank set that can't be 
///                 accommodated in the top nBest designs.
////////////////////////////////////////////////////////////////////////////////
void
MoCgs::computeCrowdingDistances(
    const vector<MultiResult>::iterator first, 
    const vector<MultiResult>::iterator last
    )
{   
    // function-scope typedef(s)
    typedef vector<MultiResult>::iterator Iterator;

    // initialize all crowding distances to 0.0
    for (Iterator i = first; i != last; ++i) {

        i->setCrowdingDistance(0.0);
    }

    // loop over each objective
    int nObj(this->getEvaluator().getNObjectives());

    for (int m = 0; m != nObj; ++m) {
        
        // sort by current objective value
        sort(first, last, MultiResult::ObjectiveCompare(m));

        // assign huge values to front and back of the sorted range (extrema of
        // the pareto front
        first->setCrowdingDistance(numeric_limits<double>::max());
        (last - 1)->setCrowdingDistance(numeric_limits<double>::max());

        // get the min and max of this objective
        double fmin(first->getObjectives()[m]);
        double fmax((last - 1)->getObjectives()[m]);

        for (Iterator i = first + 1; i != last - 1; ++i) {

            double fprev((i - 1)->getObjectives()[m]);
            double fnext((i + 1)->getObjectives()[m]);

            i->setCrowdingDistance(
                i->getCrowdingDistance() + (fnext - fprev) / (fmax - fmin)
                );
        }
    }

}   // computeCrowdingDistances


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the final result set based on the current set of best
///         designs. 
///         
///         The final result set will, at the end of each
///         iteration, contain all non-dominated designs found thus far. This 
///         quantity may be greater than that of nBest (the max allowed quantity
///         in the set of "best" designs). 
///         If no feasible designs have been identified, the final result set
///         will contain the set of non-dominated infeasible designs. It shall
///         never contain a mix of feasible and infesible designs.
///
/// \param  finalResults    A reference to the set of final results.
////////////////////////////////////////////////////////////////////////////////
void
MoCgs::updateFinalResultSet(
    unordered_set<MultiResult, KnownPoint::Hash>& finalResults
    )
{
    // vector of iterators to all new non-dominated desings
    vector<vector<MultiResult>::const_iterator> newNonDominatedDesigns;


    // if any feasible designs have been identified in bestSet_
    if (bestSet_.begin()->isFeasible()) {

        // if no feasible designs exist in the final result set, clear it out
        // because they'll all be replaced by the current set of non-dominated
        // feasible points
        if (finalResults.size() > 0 && !finalResults.begin()->isFeasible()) {

            finalResults.clear();
        }

    } else {
        
        // no feasible designs have been discovered. need to compute the
        // domination rank of the infeasible points
        this->updateDominationRank(false);
    }
    

    // process all non-dominated designs in the bestSet_
    vector<MultiResult>::const_iterator i = bestSet_.begin();

    while (i->getDominationRank() == 0 && i != bestSet_.end()) {
            
        // flag for whether to insert i after checking for domination by
        // all points currently in the final result set
        bool isNonDominated = true;

        // check if this design is not already in the set of final results
        if (finalResults.find(*i) == finalResults.end()) {
                
            unordered_set<MultiResult, KnownPoint::Hash>::const_iterator j;
            j = finalResults.begin();

            while (j != finalResults.end()) {
                    
                if (isDominatedBy(*j, *i)) {
                        
                    j = finalResults.erase(j);
                        
                } else if (isDominatedBy(*i, *j)) {
                        
                    isNonDominated = false;
                    ++j;

                } else {

                    ++j;
                }
            }
        }

        // if the design is non-dominated, add it to the set of new 
        // non-dominated designs
        if (isNonDominated) {
                
            newNonDominatedDesigns.push_back(i);
        }

        ++i;
    }


    // add all new non-dominated designs to the set of final results
    for (vector<MultiResult>::const_iterator i : newNonDominatedDesigns) {

        finalResults.insert(*i);
    }

}   // updateFinalResultSet


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs