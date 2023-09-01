//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_WELDBEAM_CPP
#define EVALUATORS_INCLUDE_WELDBEAM_CPP


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"

#include <map>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class Weldbeam : public Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors and Operators
    Weldbeam(bool trackingOn = false);

    
    // Member Functions
    virtual
    std::vector<double>
    evaluateDesign(const std::vector<double>& xIndices);

    virtual
    int
    getNObjectives() const;

    virtual
    int
    getNConstraints() const;


};  // class Weldbeam


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators


#endif  // EVALUATORS_INCLUDE_WELDBEAM_CPP