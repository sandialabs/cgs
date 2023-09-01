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
///  \brief Contains the definition of the DiscreteDistribution class and the
///         implementation of DiscreteDistribution template function(s).
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_DISCRETEDISTRIBUTION_HPP
#define CGS_INCLUDE_DISCRETEDISTRIBUTION_HPP


// INCLUDES
// #############################################################################
#include <random>
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A probability distribution that generates random integers according
///         to a discrete distribution.
///
///         This class is a wrapper around the \c std::discrete_distribution. It
///         exists to enable users to easily manipulate the weights of the 
///         underlying distribution.
////////////////////////////////////////////////////////////////////////////////
class DiscreteDistribution
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    DiscreteDistribution();

    explicit
    DiscreteDistribution(const std::vector<double>& weights);
    
    explicit
    DiscreteDistribution(int domainSize, double pseudocount = 0.0);


    // OPERATORS
    template <class T>
    int
    operator() (T& generator);
    

    // MEMBER FUNCTIONS
    std::vector<double>
    getProbabilities() const;

    double
    getProbability(int index) const;

    void
    adjustWeight(int index, double delta);

    void
    clearWeights();

    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    void 
    updateDistribution();

    bool
    weightsAreValid() const;

    bool
    weightsAreAllZero() const;


    void
    clearWeightsAndWarn();


    // DATA MEMBERS

    /// \brief  The underlying discrete distribution that this wrapper class
    ///         manipulates and draws samples from.
    std::discrete_distribution<int> distribution_;

    /// \brief  The weights of the distribution.
    ///
    ///         The probability of each of the \a n integers represented by this
    ///         distribution is that integer's weight divided by the sum of all
    ///         of the weights.
    std::vector<double> weights_;


};  // class DiscreteDistribution


// TEMPLATE DEFINITIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Samples the distribution using a random engine from the standard
///         library.
///
/// \param  generator   Reference to one of the a standard random number 
///                     engines, e.g. mt19937, default_random_engine, etc. 
///
/// \return             The discrete value that has been sampled.
////////////////////////////////////////////////////////////////////////////////
template <class T>
int 
DiscreteDistribution::operator()(T& generator) 
{
    return distribution_(generator);

}   // operator()


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet


#endif  // CGS_INCLUDE_DISCRETEDISTRIBUTION_HPP