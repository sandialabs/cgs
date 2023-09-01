//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_EVALUATOR_HPP
#define EVALUATORS_INCLUDE_EVALUATOR_HPP


// INCLUDES
// #############################################################################
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors
    explicit
    Evaluator(bool trackingOn = false);

    virtual
    ~Evaluator() {};


    // Member Functions
    virtual
    std::vector<double>
    evaluateDesign(const std::vector<double>& xValues) = 0;

    virtual
    int
    getNObjectives() const = 0;

    virtual
    int
    getNConstraints() const = 0;

    const std::vector<double>&
    getRateOfConvergence() const;

    void
    clearRateOfConvergence();


    // PROTECTED DECLARATIONS
    // =================================
protected:
    // Member Functions
    void
    trackRateOfConvergence(const std::vector<double>& fout);
    

    // Data Members
    bool trackingOn_;

    std::vector<double> rateOfConv_;

        
};  // class Evaluator


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators


#endif  // EVALUATORS_INCLUDE_EVALUATOR_HPP