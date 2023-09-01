//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the definition of the MoCgs class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_MOCGS_HPP
#define CGS_INCLUDE_MOCGS_HPP


// INCLUDES
// #############################################################################
#include "Domain.hpp"
#include "MultiResult.hpp"
#include "Solver.hpp"

#include <random>
#include <string>
#include <unordered_set>
#include <vector>

// FORWARD DECLARATIONS
// #############################################################################
namespace evaluators {
    class Evaluator;
}

// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace cgs {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A multi-objective CGS solver class.
////////////////////////////////////////////////////////////////////////////////
class MoCgs : public Solver<MultiResult>
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS 
    MoCgs(
        evaluators::Evaluator& evaluator,
        const std::vector<bayesnet::Domain>& domains,
        const std::mt19937::result_type& rngSeed = std::mt19937::default_seed,
        const std::string& graphFilename = ""
        );

    // MEMBER FUNCTIONS
    virtual
    void
    solve(int nBest, int nInitRandom, int batchSize = 1, int parentLimit = 0);
    
    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    bool
    anyConvergersSatisfied(
        long long evalCount, long iterationCount, double elapsedTime
        );

    void
    updateDominationRank(bool ignoreFeasible = true);

    void
    sortLastFrontByCrowdingDistance(int nBest);

    void
    computeCrowdingDistances(
        const std::vector<MultiResult>::iterator first, 
        const std::vector<MultiResult>::iterator last
        );

    void
    updateFinalResultSet(
        std::unordered_set<MultiResult, KnownPoint::Hash>& finalResults
        );


    // DATA MEMBERS

    /// \brief  The set of \a N best solutions found thus far. The classifier
    ///         is trained with these solutions given the "best" class label.
    std::vector<MultiResult> bestSet_;

    
};  // class MoCgs


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs

#endif  // PROJECT_PATH_CLASSNAME_HPP