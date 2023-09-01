//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Contains the definition of the ProbabilityTable class and its
///         template function(s).
////////////////////////////////////////////////////////////////////////////////

// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_PROBABILITYTABLE_HPP
#define CGS_INCLUDE_PROBABILITYTABLE_HPP


// INCLUDES
// #############################################################################
#include "DiscreteDistribution.hpp"

#include <stdexcept>
#include <unordered_map>
#include <vector>


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class DirectedGraph;
    class Domain;
}


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {

    
// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents conditional probability table for a discrete variable.
///
///         A ProbabilityTable is a conditional probability table (CPT) of a
///         discrete random variable.  It stores the conditional probabilities
///         of a single variable with respect to the variables on which it is 
///         conditioned (i.e., its parent nodes in a Bayesian
///         network). The variable and its parent variables are identified by 
///         integer values in the range 0 to n-1, where n is the total number of
///         variables under consideration.
////////////////////////////////////////////////////////////////////////////////
class ProbabilityTable
{
	// PUBLIC DECLARATIONS
	// =================================
public:
	// STRUCTORS 
	ProbabilityTable(
        const std::vector<Domain>& domains,
        const DirectedGraph& graph,
        int xLabel,
        double pseudocount = 0.0
        );
        

	// MEMBER FUNCTIONS
    double
    getProbability(const std::vector<int>& x) const;

    void
    adjustWeight(const std::vector<int>& x, double delta);

    template <class T>
    int
    sampleDistribution(const std::vector<int>& x, T& generator);


	// PRIVATE DECLARATIONS
	// =================================
private:
    // MEMBER FUNCTIONS
    void
    makeKey(const std::vector<int>& x, std::vector<int>& key) const;
    

    // NESTED COMPOUNDS

    /// \brief  A functor that computes the hash value of a feature vector.
    struct XHash {
        size_t operator() (const std::vector<int>& x) const;
    };


	// DATA MEMBERS
    
    /// \brief  The parent variables of this table's primary variable.
    std::vector<int> parents_;
    
    /// \brief  The sizes of the domains (number of possible states) of the 
    ///         parent variables. 
    std::vector<int> parentDomainSizes_;
    
    /// \brief  The core structure for storing the conditional probabilities.
    ///
    ///         This member is an associative container of DiscreteDistributions
    ///         on the primary variable that are keyed by specific 
    ///         instantiations its parent variables. 
    std::unordered_map<std::vector<int>, DiscreteDistribution, XHash> prbTable_;

    /// \brief  The label of this table's primary variable.
    int xLabel_;


};	// class ProbabilityTable


// TEMPLATE DEFINITIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Samples a point from a distribution in the probability table, given
///         values of the parents. Assumes ancestral sampling (i.e. parent
///         variables have already been sampled prior to calling.)
///
/// \param  x           The current design vector in index space. contains
///                     the already-sampled values of the parent variables.
/// \param  generator   Reference to a standard random number engine, 
///                     e.g. mt19937, default_random_engine, etc.
///
/// \return     The discrete value that has been sampled.
////////////////////////////////////////////////////////////////////////////////
template <class T>
int
ProbabilityTable::sampleDistribution
    (const std::vector<int>& x, T& generator) 
{
    // make a key for this design vector
    std::vector<int> key;
    this->makeKey(x, key);

    try {
        if (prbTable_.find(key) == prbTable_.end()) {

            throw std::domain_error("Warning (ProbabilityTable.sample"
                "Distribution): Key value not found. Nearest in-range key value"
                " used.");
        }

    } catch (std::domain_error e) {
        std::cerr << e.what() << std::endl;

        // push out of range key values back in range
        for (size_t i = 0; i != key.size(); ++i) {

            size_t maxIndex = parentDomainSizes_[i] - 1;

            if (key[i] > maxIndex) {

                key[i] = static_cast<int>(maxIndex);

            } else if (key[i] < 0) {

                key[i] = 0;
            }
        }
    }

    // sample the distribution
    return prbTable_.at(key)(generator);

}   // sampleDistribution


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet


#endif	// CGS_INCLUDE_PROBABILITYTABLE_HPP