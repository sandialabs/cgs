//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Contains the implementation of the ProbabilityTable class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "ProbabilityTable.hpp"

#include "DirectedGraph.hpp"
#include "Domain.hpp"
#include "utilities.hpp"

#include <iostream>
#include <stdexcept>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::domain_error;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief Instantiates a ProbabilityTable with specified arguments.
///
/// \param  domains     The domains of all variables in the problem.
/// \param  graph       The directed graph representing node connections in
///                     the Bayesian network.
/// \param  xLabel      The label of the primary variable of this table.
/// \param  pseudocount An initial count that may be used to initialized the
///                     weights of all DiscreteDistributions in this probability
///                     table.
////////////////////////////////////////////////////////////////////////////////
ProbabilityTable::ProbabilityTable(
    const vector<Domain>& domains,
    const DirectedGraph& graph,
    int xLabel,
    double pseudocount
    )
    :   xLabel_(xLabel), parents_(graph.getParents(xLabel))
{
    // get the parent domain sizes
    for (size_t i = 0; i != parents_.size(); ++i) {

        parentDomainSizes_.push_back(domains[parents_[i]].getSize());
    }

    if (!hasParents(xLabel, graph)) {
        // this variable has no parents. there is only one entry in the
        // probability table map. its key is a vector with one zero in it
        vector<int> key(1, 0);

        prbTable_[key] = 
            DiscreteDistribution(
                domains[xLabel].getSize(), 
                pseudocount
                );

    } else {
        // this variable has parents. the probability table map keys will be all
        // permutations (with repetition) of the parent domain indices

        // generate permutations
        vector< vector<int> > keys =
            utilities::algorithms::generatePermutationsWithRepetition(
                parentDomainSizes_
                );

        // add all permutations as keys to the probability table container
        for (size_t i = 0; i != keys.size(); ++i) {

            prbTable_[keys[i]] = 
                DiscreteDistribution(
                    domains[xLabel].getSize(),
                    pseudocount
                    );
        }
    }

}   // ProbabilityTable


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets a conditional probability given a design vector. 
///
///         This function returns the probability of this table's variable, as 
///         identified by its index, given the states of the variables on which 
///         it is conditioned.
///
/// \param  x   The design vector in index space.
///
/// \return     The probability of xi given the values of its parents.
////////////////////////////////////////////////////////////////////////////////
double
ProbabilityTable::getProbability(const vector<int>& x) const
{
    // make the probability table key
    vector<int> key;
    this->makeKey(x, key);

    // get the probability from the table
    double p;

    // check if this key exists in the table
    if (prbTable_.find(key) != prbTable_.end()) {
        // key value was found. get probability...
        p = prbTable_.at(key).getProbability(x[xLabel_]);

    } else {
        // key value was not found. return probability of zero...
        p = 0.0;
    }
       
    return p;

}   // getProbability


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adjusts one of the distribution's weights in this ProbabilityTable.
///         Changes the appropriate distribution based on the value of this 
///         table's variable (xi) and the values of its parents.
///         Used when adding or removing a point from a bayesian network
///         classifier, for example.
///
/// \param  x       A design vector in index space.
/// \param  delta   The amount to adjust the distribution's weight.
////////////////////////////////////////////////////////////////////////////////
void
ProbabilityTable::adjustWeight(const vector<int>& x, double delta)
{
    // make the key for this design
    vector<int> key;
    this->makeKey(x, key);

    try {
        if (prbTable_.find(key) == prbTable_.end()) {

            throw domain_error("Warning (ProbabilityTable.adjustWeight): Key "
                "value not found. Weight not adjusted.");
        }

        // adjust the weight in the appropriate distribution
        prbTable_[key].adjustWeight(x[xLabel_], delta);

    } catch (domain_error e) {

        std::cout << e.what() << std::endl;
    }

}   // adjustWeight


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Makes a key for the \link prbTable_ \endlink data 
///         member based on the values of a design vector.
///
/// \param  x       The design vector in index space.
/// \param  key     The key that needs to be modified.
////////////////////////////////////////////////////////////////////////////////
void
ProbabilityTable::makeKey(const vector<int>& x, vector<int>& key) const
{
    long nParents(parents_.size());
    
    if (nParents == 0) {
        // xi has no parents
        key.push_back(0);

    } else {
        // xi has parents
        key.resize(nParents);

        for (long i = 0; i != nParents; ++i) {
            key[i] = x[parents_[i]];
        }
    }

}   // makeKey


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a hash value of a vector of integers.
///
/// \param  x   The vector of ints to compute a hash value for.
///
/// \return     The hash value.
////////////////////////////////////////////////////////////////////////////////
size_t
ProbabilityTable::XHash::operator() (const vector<int>& x) const
{
    return utilities::hash::hash_range(x.begin(), x.end());
    
}   //  XHash::operator()


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet