//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "Deceptive3.hpp"


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       Deceptive3

Function:   Instantiates Deceptive3 evaluator class

Params:     trackingOn..default = false. sets whether to track rate of
                        convergence
================================================================================
*/
Deceptive3::Deceptive3(bool trackingOn)
    :   Evaluator(trackingOn)
{
}   // Deceptive3


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       evaluateDesign

Function:   Evaluates the objective and constraints of a design vector.

Params:     xValues.....the design to evaluate

Returns:    (vector<double>)  the objective and constraint values
================================================================================
*/
vector<double>
Deceptive3::evaluateDesign(const vector<double>& xValues)
{
    // convert the input design vector to an integral type
    vector<int> x(xValues.size());
    size_t itemNo(0);

    for (vector<double>::const_iterator i = xValues.begin(); 
            i != xValues.end(); 
            ++i) {

        x[itemNo++] = static_cast<int>(*i);
    }

    // compute the objective function value
    double f(0);

    size_t i = 0;
    while (i < x.size()) {
        
        // sum the next 3 values in the design
        int u(0);
        for (int j = 0; j != 3; ++j) {
            // prevent out of bounds errors
            if (i < x.size()) {
                u += x[i++];
            }
        }

        // add the subfunction value to the running total of the objective value
        if (u == 0) {
            f += 0.9;
        } else if (u == 1) {
            f += 0.8;
        } else if (u == 2) {
            f += 0.0;
        } else {
            f += 1.0;
        }
    }

    // invert f for a minimization solver
    f *= -1;

    // return vector of objective and constraint values [f, g0, g1, ...]
    vector<double> ret(1, f);

    // track rate of convergence
    if (trackingOn_) {
        
        this->trackRateOfConvergence(ret);
    }

    return ret;

}   // evaluateDesign

/*
================================================================================
Name:       getNObjectives

Function:   Returns the number of objectives of this optimization problem

Params:     (none)

Returns:    (int)   the number of objectives
================================================================================
*/
int
Deceptive3::getNObjectives() const
{
    return 1;

}   // getNObjectives

/*
================================================================================
Name:       getNConstraints

Function:   Returns the number of constraints of this optimization problem

Params:     (none)

Returns:    (int)   the number of constraints
================================================================================
*/
int
Deceptive3::getNConstraints() const
{
    return 0;

}   // getNConstraints


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators