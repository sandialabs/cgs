//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the SoCgs class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "SoCgs.hpp"

#include "Domain.hpp"
#include "Evaluators\include\Evaluator.hpp"
#include "TrainingPoint.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using bayesnet::Domain;
using bayesnet::TrainingPoint;

using std::chrono::system_clock;
using std::chrono::duration;
using std::cout;
using std::endl;
using std::max;
using std::numeric_limits;
using std::sort;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a Single Objective CGS (SoCgs) solver objecte.
///
/// \param  evaluator       Objective function evaluator object
/// \param  domains         Vector of variable domains
/// \param  rngSeed         Random number generator seed
/// \param  graphFilename   Filename which contains an initial graph edge list.
////////////////////////////////////////////////////////////////////////////////
SoCgs::SoCgs(
    evaluators::Evaluator& evaluator,
    const std::vector<Domain>& domains,
    const std::mt19937::result_type& rngSeed,
    const string& graphFilename
    )
    :   Solver(evaluator, domains, rngSeed, graphFilename),
        targetObjective_(-std::numeric_limits<double>::max())
        
{
}   // SoCgs


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Main CGS algorithm single-objective solver function
///
/// \param  nBest           Number of designs in the "best" or "good" class
/// \param  nInitRandom     Number of random designs to generate at the
///                         beginning of a solution run (i.e., before solutions
///                         begin being generated with the guidance of the
///                         classifier distributions)
/// \param  batchSize       The number of solutions sampled and evaluated on 
///                         each iteration between classifier updates.
/// \param  parentLimit     The upper limit on the number of
///                         parents a node is allowed to have in the 
///                         Bayesian network. If set to 0, no automatic
///                         network learning will occur.
////////////////////////////////////////////////////////////////////////////////
void
SoCgs::solve(int nBest, int nInitRandom, int batchSize, int parentLimit)
{
    // initializations
    SingleResult finalResult(
        TrainingPoint(), 
        vector<double>(1, -numeric_limits<double>::max()),
        numeric_limits<double>::max()
        );
                    
    long long evalCount(0);
    long iterationCount(0);
    system_clock::time_point startTime(system_clock::now());
    double elapsedTime(0.0);
    std::list<double> trackedObjectives;

    DesignSet newDesigns;
    

    // set classifier priors as a fixed uniform discrete distribution
    this->fixClassifierPriors();


    // begin solution process
    while (!this->anyConvergersSatisfied(
        evalCount, iterationCount, elapsedTime, finalResult, trackedObjectives
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

        // process the new points
        DesignSet::iterator i;
        for (i = newDesigns.begin(); i != newDesigns.end(); ++i) {
            
            // evaluate the design
            vector<double> evalResult(this->evaluateX(*i));
            ++evalCount;
            double fNew = evalResult[0];
            vector<double> gNew(evalResult.begin() + 1, evalResult.end());

            // construct a known point object and update violation stats
            SingleResult newResult(TrainingPoint(*i, 1), gNew, fNew);
            bestSet_.push_back(newResult);

            if (!newResult.isFeasible()) {
                this->updateViolationStats(newResult);
                this->updateBestSolutionViolations(bestSet_);
            }

            // sort the set of best results
            sort(
                bestSet_.begin(), 
                bestSet_.end(), 
                SingleResult::ObjectiveCompare()
                );

            // update classifier with the new evaluated point
            this->updateClassifier(bestSet_, nBest);

            // trim any additional points off of set of best results
            if (bestSet_.size() > nBest) {
                bestSet_.erase(bestSet_.begin() + nBest, bestSet_.end());
            }
            
            // update the current best known design
            finalResult = *(bestSet_.begin());

            // update tracked objectives
            if (this->getNTrackedEvaluations() != 0) {

                trackedObjectives.push_front(finalResult.getObjective());

                while (
                    trackedObjectives.size() > this->getNTrackedEvaluations()
                    ) {

                    trackedObjectives.pop_back();
                }
            }

            // print progress
            cout << "f: " << finalResult.getObjective() << ", g: ";
            vector<double>::size_type i = 0;
            for (i = 0; i != finalResult.getConstraints().size(); ++i) {
                cout << finalResult.getConstraints()[i] << ", ";
            }
            cout << "evaluations: " << evalCount << endl;
        }
        
        // if network learning is enabled, tell the classifier to learn the
        // network based on the latest training data and update the classifier
        if (parentLimit > 0) {
            
            this->learnBayesianNetwork(parentLimit);
        }

        // increment the iteration count
        iterationCount += 1;

        // update the elapsed time
        elapsedTime = duration<double>(system_clock::now() - startTime).count();
    }

    //classifier_.printBayesianNetwork();
    //getchar();

}   // solve


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets a target objective function value as a convergence criterion.
///
/// \param  targetObjective     The target objective function value. the main
///                             loop of the solve function terminates when a 
///                             design with this objective value is found.
////////////////////////////////////////////////////////////////////////////////
void
SoCgs::setTargetObjective(double targetObjective)
{
    targetObjective_ = targetObjective;

}   // setTargetObjective


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if any of the convergence criteria are satisfied. Returns
///         true when any one is satisfied.  I.e., the solver algorithm
///         terminates when the first converger is satisfied.
///
/// \param  evalCount           Current number of function evaluations.
/// \param  iterationCount      Current number of CGS iterations.
/// \param  elapsedTime         Current solution time in seconds.
/// \param  bestDesign          The current best known design.
/// \param  trackedObjectives   Objective function values tracked over a number
///                             of evaluations.
///
/// \return     Whether any stopping criterion is satisfied.
////////////////////////////////////////////////////////////////////////////////
bool
SoCgs::anyConvergersSatisfied(
    long long evalCount, 
    long iterationCount,
    double elapsedTime,
    const SingleResult& bestDesign,
    const std::list<double>& trackedObjectives
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

    // TARGET OBJECTIVE FUNCTION STOPPING CRITERION
    if (bestDesign.getObjective() <= targetObjective_ && 
        bestDesign.isFeasible()
        ) {
    
        return true;
    }


    // TRACKED PERCENT CHANGE STOPPING CRITERION

    // if number of tracked evalutions is zero, skip
    if (this->getNTrackedEvaluations() != 0) {

         if (trackedObjectives.size() < this->getNTrackedEvaluations()) {

            return false;

        } else {

            double newF = bestDesign.getObjective();
            double oldF = trackedObjectives.back();
            double percentChange;

            if (newF == oldF) {

                percentChange = 0.0;

            } else {

                percentChange = abs(newF - oldF) / max(abs(newF),abs(oldF));
            }

            if (percentChange <= this->getPercentChangeTarget()) {

                return true;
            }
        }
    }
   
    // return false if none returned true
    return false;
        
}   // anyConvergersSatisfied


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs