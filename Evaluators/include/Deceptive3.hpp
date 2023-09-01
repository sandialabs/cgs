//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_DECEPTIVE3_HPP
#define EVALUATORS_INCLUDE_DECEPTIVE3_HPP


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"

#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class Deceptive3 : public Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors 
    Deceptive3(bool trackingOn = false);

    // Member Functions
    virtual
    std::vector<double>
    evaluateDesign(const std::vector<double>& xValues);

    virtual
    int
    getNObjectives() const;

    virtual
    int
    getNConstraints() const;

    
};  // class Deceptive3


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators

#endif  // EVALUATORS_INCLUDE_DECEPTIVE3_HPP