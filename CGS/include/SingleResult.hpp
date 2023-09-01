//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the definition of the SingleResult class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_SINGLERESULT_HPP
#define CGS_INCLUDE_SINGLERESULT_HPP


// INCLUDES
// #############################################################################
#include "KnownPoint.hpp"

#include <vector>


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class TrainingPoint;
}


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents a design that has been evaluated in a single-objective
///         optimization context.
///
///         A SingleResult contains the design and result of evaluating a 
///         single objective evaluator. This class is derived from KnownPoint 
///         and is intended for single objective optimization.
////////////////////////////////////////////////////////////////////////////////
class SingleResult : public KnownPoint
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS 
    SingleResult();

    SingleResult(
        const bayesnet::TrainingPoint& trainingPoint,
        const std::vector<double>& constraintValues,
        double objectiveValue
        );
    
    // MEMBER FUNCTIONS
    double
    getObjective() const;

    
    // FUNCTORS
    /// \brief  A functor for comparing SingleResult objects by their objective
    ///         function values.
    struct ObjectiveCompare
    {
        bool operator() (
            const SingleResult& lhs, const SingleResult& rhs
            ) const;
    };

    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // DATA MEMBERS

    /// \brief The objective function value of this design.
    double objectiveValue_;
    
    
};  // class SingleResult


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs

#endif  // CGS_INCLUDE_SINGLERESULT_HPP