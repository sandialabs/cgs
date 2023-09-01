//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

// FILE CONTENTS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the definitions of general utility namespace functions. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "CGS\include\utilities.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::any_of;
using std::cout;
using std::endl;
using std::overflow_error;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace utilities {

// IN-NAMESPACE FORWARD DECLARATIONS
// #############################################################################
namespace {
    
    void 
    generatePermutations(
        const std::vector<int>& maxima,
        std::vector<int> minima,
        std::vector< std::vector<int> >& result,
        size_t index = 0
        );

}   // namespace (unnamed)


// FUNCTION DEFINITIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates all permutations (with repetition allowed) of a set of
///         integer variables with integer domains.
///
///         E.g.: Given the domains: x0 {0, 1}, x1 {0, 1}
///         This function produces the following: (0, 0) (0, 1) (1, 0) (1, 1)
///
/// \param  domainSizes     A vector of the sizes of each variable's domain.
///
/// \return     A vector of all permutations.
////////////////////////////////////////////////////////////////////////////////
vector< vector<int> >
algorithms::generatePermutationsWithRepetition(
    const std::vector<int>& domainSizes 
    )
{
    // warn the caller if input is invalid (zero or negative domain sizes)
    if (any_of(
        domainSizes.begin(), domainSizes.end(), [](int i){return i <= 0;})
        ) {
        cout << "Warning (utilities.cpp): zero or negative domain size given "
            "to permutation generator. Empty vector returned." << endl;

        return vector< vector<int> >();
    }
    
    // if no maxima are given, return an empty vector
    if (domainSizes.size() == 0) {

        return vector< vector<int> >();
    }
    
    // maxima, minima and result vector
    vector<int> maxima;
    vector<int>::const_iterator i;

    for (i = domainSizes.begin(); i != domainSizes.end(); ++i) {
        maxima.push_back(*i - 1);
    }

    vector<int> minima(domainSizes.size(), 0);
    vector< vector<int> > result;

    // call the recursive function
    generatePermutations(maxima, minima, result);
    

    return result;

}   // generatePermutationsWithRepetition


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the factorial of an integer
///
/// \param  n   The integer to compute a factorial of.
///
/// \return     The factorial.
////////////////////////////////////////////////////////////////////////////////
double
algorithms::factorial(size_t n)
{
    if (n == 0)
        return 1;

    double ret = 1;
    for (size_t i = 1; i <= n; ++i) {

        double temp = ret;
        ret *= i;

        if (ret < temp || ret == std::numeric_limits<double>::infinity()) {

            throw overflow_error("Error (utilities.cpp): factorial attempted "
                "to compute a floating point number that is larger than the "
                "maximum allowed.");
        }

    }

    return ret;

}   // factorial


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the quotient of two factorials of integers
///
/// \param  numArg      The numerator factorial argument.
/// \param  denomArg    The denominator factorial argument.
///
/// \return     The quotient of two factorials: (numArg)! / (denomArg)!
////////////////////////////////////////////////////////////////////////////////
double
algorithms::divideFactorials(size_t numArg, size_t denomArg) 
{
    if (numArg == 1 || numArg == 0) {

        // e.g. 1!/6!
        return 1 / factorial(denomArg);

    } else if (denomArg == 1 || denomArg == 0) {
        
        // e.g. 6!/1!
        return factorial(numArg);

    } else if (numArg == denomArg) {
        
        // e.g. 6!/6!
        return 1.0;

    } else if (numArg > denomArg) {

        // e.g. 6!/4! = 6*5 (the 4!'s cancel)
        size_t diff = numArg - denomArg;
        double ret(static_cast<double>(numArg));

        while (diff > 1) {
            
            ret *= ret - 1.0;
            --diff;
        }

        return ret;

    } else {
        
        // e.g. 4!/6! = 1/(6*5) (the 4!'s cancel)
        size_t diff = denomArg - numArg;
        double ret(static_cast<double>(denomArg));

        while (diff > 1) {
            
            ret *= ret - 1.0;
            --diff;
        }

        return 1 / ret;
    }

}   // divideFactorials


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the logarithm (base 2) of a factorial of an integer.
///
/// \param  n   The integer to compute a log factorial of.
///
/// \return     The log of the factorial
////////////////////////////////////////////////////////////////////////////////
double
algorithms::logFactorial(size_t n)
{
    // shortcut: if n = 0 or 1, return 0 (b/c log(1) = 0)
    if (n == 0 || n == 1) {
        
        return 0.0;
    }

    // log(n!) = log(1) + log(2) + ... + log(n)
    size_t i = 1;
    double ret(0);

    while (i <= n) {

        ret += std::log2(i);
        ++i;
    }

    return ret;

}   // logFactorial


////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates a random seed based on the current time.
///
/// \return     The random seed.
////////////////////////////////////////////////////////////////////////////////
unsigned int
algorithms::generateRandomSeed()
{
    typedef std::chrono::high_resolution_clock myclock;
    
    myclock::time_point now(myclock::now());
    
    return static_cast<unsigned>(now.time_since_epoch().count());
}


// FILE SCOPE FUNCTIONS
// #############################################################################
namespace {
////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates all permutations (with repetition allowed) of a set of
///         integer variables with integer domains.
///         This function is recursive.
///
/// \param  maxima          The maximal values of each variable.
/// \param  currPerm        the current permutation. original caller should
///                         supply the minimum values of the domain.
/// \param  result          ref to the container to which the result will
///                         be written ( vector<vector<int> >).
/// \param  index           The index of the current variable. 
///                         needed for the recursion to work properly, and
///                         should generally not be supplied by the user.
////////////////////////////////////////////////////////////////////////////////
void 
generatePermutations(
    const vector<int>& maxima, 
    vector<int> currPerm, 
    vector< vector<int> >& result, 
    size_t index
    )
{
    int max = maxima[index];
    size_t last = maxima.size() - 1;

    for (int i = currPerm[index]; i <= max; ++i) {

        currPerm[index] = i;

        if (index == last) {

            result.push_back(currPerm);

        } else {
            
            generatePermutations(
                maxima,  currPerm, result, index + 1
                );
        }
    }

}   // generatePermutations

}   // namespace (unnamed)

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace utilities