//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the definition of the SoCgs class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_SOCGS_HPP
#define CGS_INCLUDE_SOCGS_HPP


// INCLUDES
// #############################################################################
#include "Domain.hpp"
#include "SingleResult.hpp"
#include "Solver.hpp"

#include <list>
#include <random>
#include <string>
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
/// \brief  A single-objective CGS solver class.
////////////////////////////////////////////////////////////////////////////////
class SoCgs : public Solver<SingleResult>
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors 
    SoCgs(
        evaluators::Evaluator& evaluator,
        const std::vector<bayesnet::Domain>& domains,
        const std::mt19937::result_type& rngSeed = std::mt19937::default_seed,
        const std::string& graphFilename = ""
        );

    // Member Functions
    virtual
    void
    solve(
        int nBest,
        int nInitRandom,
        int batchSize = 1,
        int parentLimit = 0
        );
    
    void
    setTargetObjective(double targetObjective);

    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // Member Functions
    bool
    anyConvergersSatisfied(
        long long evalCount,
        long iterationCount,
        double elapsedSeconds,
        const SingleResult& bestDesign,
        const std::list<double>& trackedObjectives
        );

   
    // Data Members

    /// \brief  The set of \a N best solutions found thus far. The classifier
    ///         is trained with these solutions given the "best" class label.
    std::vector<SingleResult> bestSet_;

    /// \brief  The target objective function value. A converence criterion.
    double targetObjective_;
    
};  // class SoCgs



// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs

#endif  // CGS_INCLUDE_SOCGS_HPP