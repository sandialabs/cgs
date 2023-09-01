//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_KNAPSACK_HPP
#define EVALUATORS_INCLUDE_KNAPSACK_HPP


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"

#include <string>
#include <vector>

// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class Knapsack : public Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors 
    Knapsack(const std::string& filename, bool trackingOn = false);
    
    
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
    void
    readDataFile(const std::string& filename);


    // Data Members
    size_t nItems_;
    
    long capacity_;

    std::vector<long> values_;

    std::vector<long> weights_;
    
    
};  // class Knapsack


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators

#endif  // EVALUATORS_INCLUDE_KNAPSACK_HPP