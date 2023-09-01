//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "Knapsack.hpp"

#include <fstream>
#include <numeric>
#include <sstream>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::back_inserter;
using std::copy;
using std::ifstream;
using std::inner_product;
using std::istream_iterator;
using std::istringstream;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       Knapsack

Function:   Instantiates a Knapsack evaluator class

Params:     filename....a .txt file that contains the capacity, number of items,
                        and the values and weights of each item.
            trackingOn..default = false. sets whether to track rate of
                        convergence
================================================================================
*/
Knapsack::Knapsack(const string& filename, bool trackingOn)
    :   Evaluator(trackingOn)
{
    // read the data file and populate data members
    this->readDataFile(filename);

}   // Knapsack


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
Knapsack::evaluateDesign(const vector<double>& xValues)
{
    // convert the input design vector to an integral type
    vector<long> x(xValues.size());
    size_t itemNo(0);

    for (vector<double>::const_iterator i = xValues.begin(); 
            i != xValues.end(); 
            ++i) {

        x[itemNo++] = static_cast<long>(*i);
    }

    // compute the value and weight of this design
    long value = inner_product(x.begin(), x.end(), values_.begin(), 0);
    long weight = inner_product(x.begin(), x.end(), weights_.begin(), 0);

    // compute the objective and constraint values
    double f(value * (-1));

    // constraint 0: weight
    double g0 = static_cast<double>(capacity_ - weight);

    // return vector of objective and constraint values [f, g0, g1, ...]
    vector<double> ret;
    ret.push_back(f);
    ret.push_back(g0);

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
Knapsack::getNObjectives() const
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
Knapsack::getNConstraints() const
{
    return 1;

}   // getNConstraints


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       readDataFile

Function:   Reads a .txt file that contains the problem-specific parameters.
            The file must contain non-negative integers only and and have the 
            following format:

                [capacity][nItems]
                [value0][weight0]
                [value1][weight1]
                ...

Params:     filename....the name of the .txt file with the problem data

Returns:    (void)
================================================================================
*/
void
Knapsack::readDataFile(const string& filename)
{
    // read the data file
    ifstream ifs(filename);

    // read the capacity and the number of items
    ifs >> capacity_;
    ifs >> nItems_;
    ifs.get();

    // resize the value and weight vectors
    values_.resize(nItems_);
    weights_.resize(nItems_);

    // read the weights and the values
    string line;
    size_t itemNo(0);

    while (getline(ifs, line)) {

        istringstream iss(line);
        vector<long> temp;

        copy(
            istream_iterator<long>(iss),
            istream_iterator<long>(),
            back_inserter(temp)
            );

        values_[itemNo] = temp[0];
        weights_[itemNo] = temp[1];
          
        ++itemNo;
    }

}   // readFilename


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators