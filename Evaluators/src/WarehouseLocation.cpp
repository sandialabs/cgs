//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "WarehouseLocation.hpp"

#include <algorithm>
#include <set>

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
Name:       WarehouseLocation

Function:   Instantiates a WarehouseLocation evaluator object.

Params:     trackingOn..default = false. sets whether to track rate of
                        convergence
================================================================================
*/
WarehouseLocation::WarehouseLocation(bool trackingOn)
    :   Evaluator(trackingOn)
{
    // supply cost matrix (costs for each warehouse to supply each store)
    int cost[10][5] = {    
                        {20, 24, 11, 25, 30},
                        {28, 27, 82, 83, 74},
                        {74, 97, 71, 96, 70},
                        {2, 55, 73, 69, 61},
                        {46, 96, 59, 83, 4},
                        {42, 22, 29, 67, 59},
                        {1, 5, 73, 59, 56},
                        {10, 83, 13, 43, 96},
                        {93, 35, 63,85, 46},
                        {47, 65, 55, 71, 95}
                    };

    for (size_t i = 0; i != 10; ++i) {
        supplyCost_.push_back(
            vector<int>(cost[i], cost[i] + sizeof(cost[i]) / sizeof(cost[i][0]))
            );
    }

    // capacity vector (number of stores each warehouse can serve)
    int cap[] = {1,4,2,1,3};
    capacities_ = vector<int>(cap, cap + sizeof(cap) / sizeof(cap[0]));

    // fixed cost to open each warehouse
    costToOpen_ = 30;

}   // WarehouseLocation


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
WarehouseLocation::evaluateDesign(const vector<double>& xValues)
{
    // CONVERT INPUT TO INTEGERS
    vector<int> x(xValues.size());
    size_t storeNo(0);

    for (vector<double>::const_iterator i = xValues.begin(); 
            i != xValues.end(); 
            ++i) {

        x[storeNo++] = static_cast<int>(*i);
    }


    // COMPUTE TOTAL COST OF THIS CONFIGURATION
    double cost(0);
    size_t nStores(supplyCost_.size());

    // cost to supply each store
    for (size_t i = 0; i != nStores; ++i) {
        cost += supplyCost_[i][x[i]];
    }

    // cost to open warehouses
    std::set<int> openWarehouses(x.begin(), x.end());
    cost += openWarehouses.size() * costToOpen_;


    // CHECK FOR CONSTRAINT SATISFACTION

    // a warehouse cannot serve more stores than its capacity
    size_t nWarehouses(supplyCost_[0].size());
    vector<double> g(nWarehouses);
    int nServed;

    for (size_t i = 0; i != nWarehouses; ++i) {

        nServed = static_cast<int>(std::count(x.begin(), x.end(), i));
        g[i] = capacities_[i] - nServed;
    }

    
    // return vector of objective and constraint values [f, g0, g1, ...]
    double f(cost);
    vector<double> ret;
    ret.push_back(f);
    ret.insert(ret.end(), g.begin(), g.end());


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
WarehouseLocation::getNObjectives() const
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
WarehouseLocation::getNConstraints() const
{
    return 5;

}   // getNConstraints


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators