//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the Solver class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "Solver.hpp"

#include "Domain.hpp"
#include "KnownPoint.hpp"
#include "MultiResult.hpp"
#include "SingleResult.hpp"
#include "TrainingPoint.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using bayesnet::Classifier;
using bayesnet::Domain;
using bayesnet::TrainingPoint;

using std::max;
using std::mt19937;
using std::numeric_limits;
using std::set;
using std::string;
using std::uniform_int_distribution;
using std::unordered_set;
using std::vector;


// STATIC DATA MEMBER DEFINITIONS
// #############################################################################
template <class RT>
const int cgs::Solver<RT>::NCLASSES(2);

template <class RT>
const int cgs::Solver<RT>::WORST_CLASS(0);

template <class RT>
const int cgs::Solver<RT>::BEST_CLASS(1);


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a CGS Solver object with specified values.
///
/// \param  evaluator       Objective function evaluator object.
/// \param  domains         Vector of variable domains.
/// \param  rngSeed         Random number generator seed.
/// \param  graphFilename   Filename which contains an initial graph edge list.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
Solver<RT>::Solver(
    evaluators::Evaluator& evaluator, 
    const std::vector<Domain>& domains, 
    const std::mt19937::result_type& rngSeed,
    const string& graphFilename
    )
    :   classifier_(domains, NCLASSES, graphFilename, 1.0, true),
        evaluator_(evaluator),
        maxEvals_(std::numeric_limits<long long>::max()),
        maxIterations_(std::numeric_limits<long>::max()),
        maxTime_(std::numeric_limits<double>::max()),
        nTrackedEvaluations_(0),
        percentChangeTarget_(0.000),
        rng_(rngSeed),
        violationCounts_(evaluator.getNConstraints(), 0),
        violationSums_(evaluator.getNConstraints(), 0),
        domains_(domains)
{
};  // Solver


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the discrete variable domains.
///
/// \return The discrete variable domains.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
const vector<Domain>&
Solver<RT>::getDomains() const
{
    return domains_;

}   // getDomains


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets an upper limit on the number of function evaluations as a
///         termination criterion.
///
/// \param  maxEvals    Maximum number of function evaluations. The main 
///                     terminates when the limit is reached.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::setMaxEvals(long long maxEvals)
{
    maxEvals_ = maxEvals;

}   // setMaxEvals


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets an upper limit on the number of CGS iterations as a
///         termination criterion
///
/// \param  maxIterations       Maximum number of CGS iterations. The main loop 
///                             terminates when the limit is reached.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::setMaxIterations(long maxIterations)
{
    maxIterations_ = maxIterations;

}   // setMaxIterations


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets an upper limit on the solution wall-clock time as a termination
///         criterion.
///
/// \param  maxTime     The maximum time (in seconds) that the optimizer is
///                     allowed to run. The main loop terminates when the limit 
///                     is reached.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::setMaxTime(double maxTime)
{
    maxTime_ = maxTime;

}   // setMaxTime


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the percent change target and the number of tracked iterations
///         that are used to calculate the percent change convergence criterion.
///
///         The optimizer run function will terminate when the percent change is
///         less than or equal to the target.
///
/// \param  nTrackedEvaluations     The number of evaluations over which to 
///                                 track percent change. A values of zero
///                                 indicates that this converger is inactive.
/// \param  percentChangeTarget     The new target for percent change in 
///                                 objective value.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::setTrackedChange(
    int nTrackedEvaluations,
    double percentChangeTarget
    )
{
    nTrackedEvaluations_ = nTrackedEvaluations;
    
    percentChangeTarget_ = percentChangeTarget;

}   // setTrackedPercentChange


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the prior probability distributions in the classifier to fixed,
///         uniform, discrete distributions.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::fixClassifierPriors()
{
    classifier_.setFixedPriors(vector<double>(classifier_.getNClasses(), 1));

}   // fixClassifierPriors


// PROTECTED MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates a new candidate design by sampling the classifier
///         distributions. 
///
///         If the number of evaluations is less than nInitRandom, the new 
///         design is generated by sampling random uniform distributions. 
///         Otherwise, a design is generated by sampling the classifier's 
///         "good" class conditional distributions, i.e. P(x|'good').
///
/// \param  nInitRandom     Number of random designs to generate at the
///                         beginning of a solution run (i.e., before solutions
///                         begin being generated with the guidance of the
///                         classifier distributions).
/// \param  evalCount       The current count of function evaluations.
/// \param  newDesigns      The current set of new designs that have been
///                         sampled for evaluation on the current iteration.
///
/// \return     The new candidate design to propose for evalution.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
vector<int>
Solver<RT>::sampleNewX(
    int nInitRandom, size_t evalCount, const DesignSet& newDesigns
    )
{
    vector<int> x;
    bool isNew(false);

    while (!isNew) {

        if (evalCount < nInitRandom) {

            // generate a random design
            x = this->generateRandomDesign();

        } else {

            // sample a design from the classifier distributions
            x = classifier_.sampleXVal(BEST_CLASS, rng_);
        }
        
        // check if design has already been evaluated or if it has already been
        // sampled for evalution on this iteration
        if (!classifier_.containsPoint(TrainingPoint(x, WORST_CLASS)) &&
            !classifier_.containsPoint(TrainingPoint(x, BEST_CLASS)) &&
            newDesigns.find(x) == newDesigns.end()) {

                isNew = true;
        } 
    }

    return x;

}   // sampleNewX


////////////////////////////////////////////////////////////////////////////////
/// \brief  Determines whether to accept a proposed design for evaluation.
///
///         Determination is made by comparing the posterior probability of 
///         being "good", P('good'|x), to a random sample from a uniform 
///         distribution.
///
/// \param  xNew    The design being proposed for evalution.
///
/// \return     \b True if the design is accepted for evalution.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
bool
Solver<RT>::isAcceptedForEval(const vector<int>& xNew)
{
    // return value
    bool ret(true);

    // use classifier determine posterior probability of best class
    double pcx = classifier_.computePosteriors(xNew)[BEST_CLASS];

    // accept/reject by comparing to uniform sample
    if (pcx < uniformDistribution_(rng_)) {
        // reject
        ret = false;
    }

    return ret;

}   // isAcceptedForEval


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the objective and constraint values of a design.
///
/// \param  xi   The vector of indices that represents the candidate design.
///
/// \return     The objective and constraint values.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
vector<double>
Solver<RT>::evaluateX(const vector<int>& xi)
{
    // the solver (and the classifier) operate on indices, so we have to convert
    // x to the evaluator's domain values
    vector<double> x;
    for (size_t i = 0; i != xi.size(); ++i) {
        x.push_back(domains_[i].indexToXval(xi[i]));
    }

    return evaluator_.evaluateDesign(x);

}   // evaluateX


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the constraint violation counts and sums.
///
/// \param  newResult   The known (evaluated) point to update stats with.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::updateViolationStats(const KnownPoint& newResult)
{
    vector<double> gs = newResult.getConstraints();

    // assert that the size of violationCounts_  and violationSums_ is equal to 
    // the size of gs
#ifdef _DEBUG
    assert(gs.size() == violationCounts_.size());
    assert(gs.size() == violationSums_.size());
#endif

    for (int i = 0; i != gs.size(); i++) {

        if (gs[i] < 0) {

            violationCounts_[i] += 1;
            violationSums_[i] += (-1)*gs[i];
        }
    }

}   // updateViolationStats


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the overall violation of a KnownPoint object. 
///
///         If the object violates any constraints, the infeasibility is set to 
///         the sum of all "relative" violations (i.e. each violation is divided
///         by the average violation for that constraint of all points
///         evaluated thus far and multiplied by 100).  The result for all
///         constraints is then summed to get the total relative violation.
///
/// \param  result  The result type object whose violation is to be updated.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::updateOverallViolation(KnownPoint& result)
{
    if (!result.isFeasible()) {
        
        // design violates at least one constraint, so we need to compute a
        // non-zero violation score.
        double violation(0);
        vector<double>::const_iterator i = result.getConstraints().begin();
        vector<double> gs(result.getConstraints());
        int j = 0;
        
        for (; i != result.getConstraints().end(); ++i, ++j) {

            if (*i < 0) {

                double meanViolation = violationSums_[j] / violationCounts_[j];

                violation += (-1) * (*i) / meanViolation * 100.0;
            }
        }
        
        result.setOverallViolation(violation);

    } else {

        // design does not violate any constraints
        result.setOverallViolation(0.0);
    }

}   // udpateOverallViolation


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the overall violations of all of the infeasible points in 
///         the set of \a N best evaluated designs (if there are any).
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::updateBestSolutionViolations(vector<RT>& bestSet)
{
    // update the infeasibility metrics of all infeasible points
    vector<RT> infeasibleBests;

    for (RT& i : bestSet) {

        if (!i.isFeasible()) {
            
            this->updateOverallViolation(i);
        }
    }

}   // updateBestSolutionViolations


////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the classifier with the current set of best designs.
///
///         When this function is called, the current set of best designs may
///         be greater that the parameter \a N best. The top \a N best designs
///         are assigned a class label of 'good', and all others are assigned a
///         class label of bad. Design that already exist in the classifier
///         will either remain where they are or will be moved from the 'good'
///         side to the 'bad' side (if they are no longer among the top \a N 
///         best designs found thus far). This function assumes that the
///         incoming set of best designs is already sorted by some measure of 
///         goodness.
///
/// \param  bestSet     The sorted set of current designs to process.
/// \param  nBest       Number of points that may be classified as 'best'.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::updateClassifier(
    const vector<RT>& bestSet, const int& nBest
    )
{
    // if bestSet is less than full
    if (bestSet.size() < nBest) {
        
        // iterate over all points set of best results
        for (RT i : bestSet) {
            
            // if the point is already in the classifier, do nothing
            if (classifier_.containsPoint(i.getTrainingPoint())) {
                
                continue;

            } else {
                
                // otherwise, it is new so add it to the classifier
                classifier_.addPoint(i.getTrainingPoint());
            }
        }

    } else {
        
        // iterate over the first nBest points in set of best results
        typedef vector<RT>::const_iterator cit;

        for (cit i = bestSet.begin(); i != bestSet.begin() + nBest; ++i) {
            
            // if the point is already in the clasisifier, do nothing
            if (classifier_.containsPoint(i->getTrainingPoint())) {

                continue;

            } else {

                // otherwise, is is new so add it to the clasisifier
                classifier_.addPoint(i->getTrainingPoint());
            }
        }

        // iterate over the remaining point(s) in the set of best results
        for (cit i = bestSet.begin() + nBest; i != bestSet.end(); ++i) {

            // if the point is already in the classifier, move from good to bad
            if (classifier_.containsPoint(i->getTrainingPoint())) {

                classifier_.movePoint(
                    i->getTrainingPoint().getXValues(),
                    BEST_CLASS,
                    WORST_CLASS
                    );

            } else {
                
                // otherwise, it's a new bad point
                TrainingPoint newTp(i->getTrainingPoint());
                newTp.setClassLabel(WORST_CLASS);

                classifier_.addPoint(newTp);
            }
        }
    }
    
}   // updateClassifier


////////////////////////////////////////////////////////////////////////////////
/// \brief  Learns the Bayesian network structure from the current set of 
///         training points and updatest the classifier accordingly.
///
/// param   parentLimit     The upper bound on the number of parents a node can
///                         have in the Bayesian Network.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
void
Solver<RT>::learnBayesianNetwork(int parentLimit)
{
    vector<int> searchOrder(this->getDomains().size());
    std::iota(searchOrder.begin(), searchOrder.end(), 0);
    std::shuffle(
        searchOrder.begin(), searchOrder.end(), rng_
        );

    classifier_.learnNetworkAndUpdate(
        parentLimit, searchOrder, BEST_CLASS, BEST_CLASS
        );

}   // learnBayesianNetwork


////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether the number of function evaluations is greater than
///         or equal to the user-specified limit.
///
/// \param  evalCount   The current number of evalutions completed.
///
/// \return     \b True if the number of evalutions has reached the limit.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
bool
Solver<RT>::maxEvaluationsSatisfied(long long evalCount) const
{
    return evalCount >= maxEvals_;

}   // maxEvalutionsSatisfied


////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether the number of CGS iterations is greater than
///         or equal to the user-specified limit.
///
/// \param  iterationCount  The current number of iterations completed.
///
/// \return     \b True if the number of iterations has reached the max.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
bool
Solver<RT>::maxIterationsSatisfied(long iterationCount) const
{
    return iterationCount >= maxIterations_;

}   // maxEvalutionsSatisfied


////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether the solution time has met or exceeded the 
///         user-specified limit on max time (in seconds).
///
/// \param  elapsedTime     The current number of seconds elapsed.
///
/// \return     \b True if the elapsed time has reached the max.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
bool
Solver<RT>::maxTimeSatisfied(double elapsedTime) const
{
    return elapsedTime >= maxTime_;

}   // maxTimeSatisfied


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the number of tracked evalutions for the tracked percent change
///         converger.
///
/// \return     The number of tracked evaluations.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
int
Solver<RT>::getNTrackedEvaluations() const
{
    return nTrackedEvaluations_;

}   // getNTrackedEvaluations


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the target percent change convergence criterion.
///
/// \return     The target percent change.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
double
Solver<RT>::getPercentChangeTarget() const
{
    return percentChangeTarget_;

}   // getPercentChangeTarget


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the evaluator that lives in the base class
///
/// \return     The objective function evaluator.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
const evaluators::Evaluator&
Solver<RT>::getEvaluator() const
{
    return evaluator_;

}   // getEvaluator


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates a design by sampling from uniform a distributions.
///
/// \return     The new random design.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
vector<int>
Solver<RT>::generateRandomDesign()
{
    vector<int> ret;
    
    vector<Domain>::const_iterator i = domains_.begin();
    vector<Domain>::const_iterator end = domains_.end();

    for (i; i != end; ++i) {

        uniform_int_distribution<int> dist(0, i->getSize() - 1);

        ret.push_back(dist(rng_));
    }

    return ret;

}   // generateRandomDesign


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a hash value of a vector of integers. 
///
/// \param  x   The vector of ints to compute a hash value for.
///
/// \return     The hash value.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
size_t
Solver<RT>::XHash::operator() (const vector<int>& x) const
{
    return utilities::hash::hash_range(x.begin(), x.end());
    
}   //  XHash::operator()


// EXPLICIT TEMPLATE INSTANTIATIONS
// =============================================================================
template class Solver<SingleResult>;
template class Solver<MultiResult>;


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs