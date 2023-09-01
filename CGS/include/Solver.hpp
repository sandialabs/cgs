//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the definition of the Solver class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_SOLVER_HPP
#define CGS_INCLUDE_SOLVER_HPP


// INCLUDES
// #############################################################################
#include "Classifier.hpp"

#include "Evaluators\include\Evaluator.hpp"

#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class Domain;
}

namespace cgs {
    class KnownPoint;
}

// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A base class for the single- and multi-objective CGS solver classes.
///
///         This class is an abstract template class with explicit 
///         instantiations using the SingleResult and MultiResult result types 
///         (RTs). The derived solver classes (SoCgs and MoCgs) inherit
///         specific instantiations of the Solver template using the
///         appropriate RT: SingleResult for SoCgs and MultiResult for MoCgs.
////////////////////////////////////////////////////////////////////////////////
template <class RT>
class Solver
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    Solver(
        evaluators::Evaluator& evaluator, 
        const std::vector<bayesnet::Domain>& domains,
        const std::mt19937::result_type& rngSeed,
        const std::string& graphFilename
        );

    virtual
    ~Solver() {};

    
    // MEMBER FUNCTIONS
    virtual
    void
    solve(
        int nBest, 
        int nInitRandom,
        int batchSize,
        int parentLimit
        ) = 0;

    const std::vector<bayesnet::Domain>&
    getDomains() const;

    void
    setMaxEvals(long long maxEvals);

    void
    setMaxIterations(long maxIterations);

    void
    setMaxTime(double maxTime);

    void
    setTrackedChange(
        int nTrackedEvaluations,
        double percentChangeTarget
        );

    void
    fixClassifierPriors();


    // CONSTANTS

    /// \brief  The number of classes in the classifier.
    static const int NCLASSES;
    
    /// \brief  The label of the worst class in the classifier.
    static const int WORST_CLASS;

    /// \brief  The label of the best class in the classifier.
    static const int BEST_CLASS;


    // PROTECTED DECLARATIONS
    // =================================
protected:
    // NESTED STRUCTS

    /// \brief  A functor that computes the hash value of a feature vector.
    struct XHash {
        size_t operator() (const std::vector<int>& x) const;
    };


    // TYPEDEFS
    typedef std::unordered_set<std::vector<int>, XHash> DesignSet;


    // MEMBER FUNCTIONS
    std::vector<int>
    sampleNewX(int nInitRandom, size_t evalCount, const DesignSet& newDesigns);
    
    bool
    isAcceptedForEval(const std::vector<int>& xNew);

    std::vector<double>
    evaluateX(const std::vector<int>& xi);

    void
    updateViolationStats(const KnownPoint& newResult);

    void
    updateOverallViolation(KnownPoint& result);

    void
    updateBestSolutionViolations(std::vector<RT>& nBestResults);

    void
    updateClassifier(const std::vector<RT>& bestSet, const int& nBest);

    void
    learnBayesianNetwork(int parentLimit);

    bool
    maxEvaluationsSatisfied(long long evalCount) const;

    bool
    maxIterationsSatisfied(long iterationCount) const;

    bool
    maxTimeSatisfied(double elapsedTime) const;

    int
    getNTrackedEvaluations() const;

    double
    getPercentChangeTarget() const;

    const evaluators::Evaluator&
    getEvaluator() const;


    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    std::vector<int>
    generateRandomDesign();


    // DATA MEMBERS

    /// \brief  The Bayesian classifier for sampling and filtering new designs.
    bayesnet::Classifier classifier_;

    /// \brief  The domains of the discrete design variables.
    const std::vector<bayesnet::Domain>& domains_;

    /// \brief  The objective function evaluator.
    evaluators::Evaluator& evaluator_;

    /// \brief  The maximum function evalution termination criterion.
    long long maxEvals_;

    /// \brief  The maximum iterations termination criterion.
    long maxIterations_;

    /// \brief  The maximum wall-clock time termination criterion.
    double maxTime_;

    /// \brief  The number of tracked evaluations to use for the percent change
    ///         convergence criterion.
    int nTrackedEvaluations_;

    /// \brief  The percent change convergence criterion.
    double percentChangeTarget_;

    /// \brief  The random number engine used for any random number generation
    ///         needed throughout the optimization process.
    std::mt19937 rng_;

    /// \brief  The number of time each constraint has been violated over all
    ///         objective function evaluations.
    std::vector<long long> violationCounts_;

    /// \brief  The sum of all constraint violations tracked for each
    ///         constraint separately.
    std::vector<double> violationSums_;

    /// \brief  A random uniform distribution.
    std::uniform_real_distribution<double> uniformDistribution_;

    
};  // class Solver


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs


#endif  // CGS_INCLUDE_SOLVER_HPP