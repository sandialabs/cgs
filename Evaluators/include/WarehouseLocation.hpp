//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_WAREHOUSELOCATION_HPP
#define EVALUATORS_INCLUDE_WAREHOUSELOCATION_HPP


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class WarehouseLocation : public Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors
    WarehouseLocation(bool trackingOn = false);
    

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

    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // Member Functions
    
    // Data Members
    std::vector<std::vector<int> > supplyCost_;

    std::vector<int> capacities_;

    int costToOpen_;
    
    
};  // class WarehouseLocation


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators

#endif  // EVALUATORS_INCLUDE_WAREHOUSELOCATION_HPP