//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"

#include <algorithm>
#include <limits>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::any_of;
using std::numeric_limits;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace evaluators {


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       Evaluator

Function:   Instantiates an Evaluator object with a constructor arg.

Params:     trackingOn......default = false. boolean that sets whether rate
                            of convergence tracking is on for this object.
================================================================================
*/
Evaluator::Evaluator(bool trackingOn) : trackingOn_(trackingOn)
{ 
};  //  Evaluator


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       getRateOfConvergence

Function:   Gets the rateOfConv_ member variable.

Params:     (none)

Returns:    (vector<double>)    rateOfConv_
================================================================================
*/
const vector<double>&
Evaluator::getRateOfConvergence() const
{
    return rateOfConv_;

}   // getRateOfConvergence

/*
================================================================================
Name:       clearRateOfConvergence

Function:   Clears the vector that stores the rate of convergence data

Params:     (none)

Returns:    (void)
================================================================================
*/
void
Evaluator::clearRateOfConvergence()
{
    rateOfConv_.clear();

}   // clearRateOfConvergence


// PROTECTED MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       trackRateOfConvergence

Function:   Tracks the rate of convergence. If the evaluateDesign member 
            function is called repeatedly for a Evaluator object, the rate of 
            convergence towards optimal solutions can be tracked. If the new 
            objective function value is less than the the best found thus far, 
            the rate of convergence vector is appended with the new value. 
            Otherwise, the rate of conv. vector is appended with a copy of the 
            previously known best objective value.

Params:     fout....the output of evalX [f, g1, g2, ... , gn]

Returns:    (void)
================================================================================
*/
void
Evaluator::trackRateOfConvergence(const vector<double>& fout)
{
    // separate objective and constraints
    double f(fout.front());
    vector<double> g(fout.begin() + 1, fout.end());

    // assess feasibility
    bool isFeasible(!any_of(g.begin(), g.end(), [](double i){return i < 0;}));

    // if roc vector is empty...
    if (rateOfConv_.empty()) {

        // insert f as first element (if it is feasible)
        if (isFeasible) {

            rateOfConv_.push_back(f);

        } else {

            rateOfConv_.push_back(numeric_limits<double>::max());

        }

    } else {

        // if f is better than best currently found...
        if (f < rateOfConv_.back() && isFeasible) {

            // add it to roc vector
            rateOfConv_.push_back(f);

        } else {

            // append roc vector with best found so far
            rateOfConv_.push_back(rateOfConv_.back());

        }
    }

}   // trackRateOfConvergence


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators