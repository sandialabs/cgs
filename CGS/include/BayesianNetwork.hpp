//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

// FILE CONTENTS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the definition of the BayesianNetwork class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_BAYESIANNETWORK_HPP
#define CGS_INCLUDE_BAYESIANNETWORK_HPP


// INCLUDES
// #############################################################################
#include "DirectedGraph.hpp"

#include <string>
#include <utility>
#include <vector>


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class Domain;
    class TrainingDataset;
}

// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents a Bayesian network of discrete random variables.
///
///         A BayesianNetwork is a DirectedGraph that represents a set of 
///         discrete random variables and their conditional dependencies. The 
///         structure of the network may be specified, or it may be estimated 
///         using a set of training data. BayesianNetwork is 
///         intended to be used in conjunction with the Classifier class in 
///         this namespace.
////////////////////////////////////////////////////////////////////////////////
class BayesianNetwork : public DirectedGraph
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS 
    BayesianNetwork(
        const std::vector<Domain>& domains, 
        const TrainingDataset& trainingPoints, 
        const std::string& edgeListFilename = ""
        );

    BayesianNetwork(
        const std::vector<Domain>& domains,
        const TrainingDataset& trainingPoints,
        const DirectedGraph& graph
        );
    
    
    // MEMBER FUNCTIONS
    bool
    learnNetworkFromData(
        int parentLimit, 
        const std::vector<int>& searchOrder,
        int learningClass = -1
        );
        

    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    double
    computeNetworkScore(const AdjacencyList& parentList) const;

    std::pair<int, double>
    findBestParentToAdd(
        int child, 
        const std::vector<int>& currentParents,
        const std::vector<int>& predecessors
        ) const;

    double
    computeNodeScore(
        int child, const std::vector<int>& parents
        ) const;

    
    // DATA MEMBERS

    /// \brief  The domains of all discrete variables that this BN represents.
    const std::vector<Domain>& domains_;

    /// \brief  The classifier training points from which to learn the BN 
    ///         structure.
    const TrainingDataset& trainingPoints_;

    /// \brief  The class whose training points shall be used to learn the BN 
    ///         structure. Learning the structure from only one class's data is 
    ///         optional. If this member is equal to -1, then training points 
    ///         from all classes are used.
    int learningClass_;
    
    
};  // class BayesianNetwork


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet

#endif  // CGS_INCLUDE_BAYESIANNETWORK_HPP