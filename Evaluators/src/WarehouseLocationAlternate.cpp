//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "WarehouseLocationAlternate.hpp"

#include <algorithm>
#include <numeric>
#include <set>

// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::accumulate;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// STATIC DATA MEMBER DEFINITIONS
// #############################################################################


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       WarehouseLocationAlternate

Function:   Instantiates a WarehouseLocationAlternate evaluator object.

Params:     trackingOn..default = false. sets whether to track rate of
                        convergence
================================================================================
*/
WarehouseLocationAlternate::WarehouseLocationAlternate(bool trackingOn)
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
    costToOpen_ = 30.0;

}   // WarehouseLocationAlternate


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
WarehouseLocationAlternate::evaluateDesign(const vector<double>& xValues)
{
    // CONVERT INPUT TO 2-D DESIGN ARRAY INTEGERS
    // One row per store. One column per warehouse. I.e., a value of 1 in [i, j] 
    // indicates that store i is served by warehouse j.
    vector< vector<int> > X;
    size_t nStores(supplyCost_.size());
    size_t nWarehouses(capacities_.size());

    for (size_t i = 0; i != nStores; ++i) {

        X.push_back(vector<int>(nWarehouses));

        for (size_t j = 0; j != nWarehouses; ++j) {
            
            // map this i,j combo to the correct index of the 1-D input vector
            size_t k = i * nWarehouses + j;
            
            // cast the value to our new 2-D array
            X[i][j] = static_cast<int>(xValues[k]);
        }
    }


    // COMPUTE TOTAL COST OF THIS CONFIGURATION
    double cost(0);
    vector<double> costsToOpen(nWarehouses);

    // compute cost to supply each store
    for (size_t i = 0; i != nStores; ++i) {

        for (size_t j = 0; j != nWarehouses; ++j) {
            
            // add the cost only if the design array has a non-zero value
            // for this store/warehouse combo
            if (X[i][j] != 0) {

                cost += supplyCost_[i][j];

                costsToOpen[j] = costToOpen_;
            }
        }
    }
    
    // add the cost of opening the warehouses
    cost += accumulate(costsToOpen.begin(), costsToOpen.end(), 0.0);


    // CHECK FOR CONSTRAINT SATISFACTION

    // each store must be served by exactly one warehouse
    int nStoresNotServed(0);
    int nStoresOverServed(0);

    // a warehouse cannot serve more stores than its capacity
    vector<int> nServed(nWarehouses);

    for (size_t i = 0; i != nStores; ++i) {
        
        int servedBy(accumulate(X[i].begin(), X[i].end(), 0));

        if (servedBy == 0) {

            nStoresNotServed += 1;

        } else if (servedBy > 1) {

            nStoresOverServed += 1;
        }

        for (size_t j = 0; j != nWarehouses; ++j) {

            if (X[i][j] != 0) {

                nServed[j] += 1;
            }
        }
    }

    vector<double> capacityMinusServed(nWarehouses);

    for (size_t j = 0; j != nWarehouses; ++j) {

        capacityMinusServed[j] = capacities_[j] - nServed[j];
    }

    
    // RETURN VECTOR OF OBJECTIVE AND CONSTRAINT VALUES [f, g0, g1, ...]
    // f is the total cost of opening warehouses and serving stores
    double f(cost);

    // g0 - g4 are the warehouse capacity constriants (one per warehouse)
    // g5 and 6 are the store overserved and notserved, respectively
    vector<double> g(capacityMinusServed);
    g.push_back(-nStoresOverServed);
    g.push_back(-nStoresNotServed);

    vector<double> ret;
    ret.push_back(f);
    ret.insert(ret.end(), g.begin(), g.end());

    // TRACK RATE OF CONVERGENCE
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
WarehouseLocationAlternate::getNObjectives() const
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
WarehouseLocationAlternate::getNConstraints() const
{
    return 7;

}   // getNConstraints


// PRIVATE MEMBER FUNCTIONS
// #############################################################################


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators