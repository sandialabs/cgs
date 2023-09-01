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
///  \brief Contains the implementation of the DiscreteDistribution class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "DiscreteDistribution.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::all_of;
using std::any_of;
using std::discrete_distribution;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Default constructor.
///
///         Instantiates a DiscreteDistribution that always produces a zero.
////////////////////////////////////////////////////////////////////////////////
DiscreteDistribution::DiscreteDistribution() 
{
}   // DiscreteDistribution


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a DiscreteDistribution object with a vector of weights.
///         
///         If a vector of zeros is given, a uniform distribution is created
///         with near-zero weights (smallest possible value for type double).
///         If an empty vector is given, a distribution with only one value (0)
///         and a weight of 1.0 is created.
///
/// \param  weights     Vector of initial weights that define the probabilities 
///                     of each value in the discrete distribution.
////////////////////////////////////////////////////////////////////////////////
DiscreteDistribution::DiscreteDistribution
    (const std::vector<double>& weights)
    :   weights_(weights)
{
    // check for empty vector
    if (weights_.size() == 0) {

        weights_.push_back(1.0);
    }

    // check the weights validity. if they aren't valid, clear them to near-zero
    // (but not zero) values
    if (!this->weightsAreValid()) {
        
        // if weights are all zero, clear them to near-zero. no need to warn
        // anyone about it. otherwise, clear and warn.
        if (this->weightsAreAllZero()) {

            this->clearWeights();

        } else {
            
            this->clearWeightsAndWarn();
        }
        
    }

    // update the distribution with the initial weights
    this->updateDistribution();

}   // DiscreteDistrubution


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a DiscreteDistribution object with a domain size and
///         pseudocount.
///
///         If a vector of zeros is given, a uniform distribution is created
///         with near-zero weights (smallest possible value for type double).
///         If an empty vector is given, a distribution with only one value (0)
///         and a weight of 1.0 is created.
///
/// \param  domainSize      The number of values in the discrete distribution.
/// \param  pseudocount     The initial weights that define the discrete
///                         distribution probabilities. All weights will be 
///                         initially set to the pseudocount when this 
///                         constructor is used (i.e., uniform distribution).
////////////////////////////////////////////////////////////////////////////////
DiscreteDistribution::DiscreteDistribution(int domainSize, double pseudocount)
    :   weights_(vector<double>(domainSize, pseudocount))
{
    // check for empty vector
    if (weights_.size() == 0) {

        weights_.push_back(1.0);
    }

    // check the weights validity. if they aren't valid, clear them to near-zero
    // (but not zero) values
    if (!this->weightsAreValid()) {
        
        // if weights are all zero, clear them to near-zero. no need to warn
        // anyone about it. otherwise, clear and warn.
        if (this->weightsAreAllZero()) {

            this->clearWeights();

        } else {
            
            this->clearWeightsAndWarn();
        }
        
    }

    // update the distribution with the initial weights
    this->updateDistribution();

}   // DiscreteDistrubution


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the probability of each integer in the distribution.
///
/// \return                     Vector of probabilities
////////////////////////////////////////////////////////////////////////////////
vector<double>
DiscreteDistribution::getProbabilities() const
{
    return distribution_.probabilities();

}   // getProbabilities


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets a specified probability from the distribution
///
/// \param  index   The index of the probability to return.
///
/// \return         The probability of the given integer.
////////////////////////////////////////////////////////////////////////////////
double
DiscreteDistribution::getProbability(int index) const
{
    return distribution_.probabilities()[index];

}   // getProbability


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adjusts the weight of one value in the distribution. 
///         
///         If a invalid weights are created, they are cleared.
///
/// \param  index   The index of the value whose weight is being adjusted.
/// \param  delta   The amount to adjust the weight.
////////////////////////////////////////////////////////////////////////////////
void
DiscreteDistribution::adjustWeight(int index, double delta)
{
    weights_[index] += delta;

    if (weights_[index] < 0) {
        
        // if the new weight is negative it is invalid. clear and warn
        this->clearWeightsAndWarn();

    } else if (weights_[index] == 0.0) {
        
        // if the new weight is zero, it's possible that all weights are zero,
        // which is invalid for std::discrete_distribution. clear to near-zero
        // but no need to warn
        if (this->weightsAreAllZero()) {

            this->clearWeights();
        }
    }

    this->updateDistribution();

}   // adjust weight


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets all weights to near-zero (but not zero) values. 
///
///         A weight set of all zeros is invalid.
////////////////////////////////////////////////////////////////////////////////
void
DiscreteDistribution::clearWeights()
{
    typedef vector<double>::iterator Iter;

    for (Iter i = weights_.begin(); i != weights_.end(); ++i) {

        *i = std::numeric_limits<double>::min();
    }

}   // clearWeights


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the distribution according to the current values in the
///         \link weights_ \endlink member variable.
////////////////////////////////////////////////////////////////////////////////
void
DiscreteDistribution::updateDistribution()
{
    // construct a param_type and then construct the distribution
    discrete_distribution<>::param_type params(
        weights_.size(), 0.0, static_cast<double>(weights_.size()),
        [&](const double& x) -> double {
            return weights_[static_cast<size_t>(x)];
        }
    );

    distribution_.param(params);

}   // updateDistribution


////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if the weights are valid. 
///
///         The rules for the determining the validity of a weight set are the 
///         same as those of the \c std::discrete_distribution: at least one 
///         weight must be positive and all weights be non non-negative.
///
/// \return         \b True if weights are valid.
////////////////////////////////////////////////////////////////////////////////
bool
DiscreteDistribution::weightsAreValid() const
{
    bool hasNoNegative(!any_of(
        weights_.begin(), weights_.end(), [](double i) {return i < 0.0;}
        ));

    bool hasPositive = any_of(
        weights_.begin(), weights_.end(), [](double i) {return i > 0.0;}
        );

    return hasNoNegative && hasPositive;

}   // wieghtsAreValid


////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if all weights are exactly equal to zero.
///
/// \return         \b True if all weights are zero.
////////////////////////////////////////////////////////////////////////////////
bool
DiscreteDistribution::weightsAreAllZero() const
{
    return (all_of(
        weights_.begin(), weights_.end(), [](double i) {return i == 0.0;}
        ));
}


////////////////////////////////////////////////////////////////////////////////
/// \brief  Clears the weights and displays a warning that the weights were
///         cleared (possibly due to being invalid).
////////////////////////////////////////////////////////////////////////////////
void
DiscreteDistribution::clearWeightsAndWarn()
{
    this->clearWeights();

    std::cout << "Warning (DiscreteDistribution.cpp): Weights were cleared "
        "automatically (possibly due to an invalid weight set)." << std::endl;
}

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet