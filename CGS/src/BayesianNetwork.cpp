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
///  \brief Contains the implementation of the BayesianNetwork class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "BayesianNetwork.hpp"

#include "Domain.hpp"
#include "TrainingDataset.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_set>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::cout;
using std::endl;
using std::invalid_argument;
using std::is_permutation;
using std::mt19937;
using std::pair;
using std::string;
using std::unordered_set;
using std::vector;

using utilities::algorithms::logFactorial;
using utilities::algorithms::generatePermutationsWithRepetition;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a BayesianNetwork object. 
///
/// \param  domains             The domains of the variables in this network.
/// \param  trainingPoints      A dataset of training points that may be used to
///                             learn the structure of the network.
/// \param  edgeListFilename    A filename that contains the edges from which to
///                             build the network. Default behavior (no filename
///                             given) is to construct a "fully independent" 
///                             network (no edges between the nodes).
////////////////////////////////////////////////////////////////////////////////
BayesianNetwork::BayesianNetwork(
    const std::vector<Domain>& domains, 
    const TrainingDataset& trainingPoints, 
    const std::string& edgeListFilename
    ) 
    :   DirectedGraph(static_cast<int>(domains.size()), edgeListFilename), 
        domains_(domains), 
        trainingPoints_(trainingPoints),
        learningClass_(-1)
{
    // ensure that this graph is acyclic. if it isn't, delete all edges and warn
    //  the caller
    if (!this->isAcyclic()) {

        AdjacencyList emptyParentList;

        for (int i = 0; i != this->getNodeCount(); ++i) {

            emptyParentList[i] = vector<int>();
        }

        this->setParentList(emptyParentList);

        cout << "Warning (BayesianNetwork.cpp): The initial graph supplied to "
            "the constructor has one or more cycles. All edges were deleted." <<
            endl;
    }

}   // BayesianNetwork

////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a BayesianNetwork object using a DirectedGraph. 
///
/// \param  domains         The domains of the variables in this network.
/// \param  trainingPoints  A dataset of training points that may be used
///                         to learn the structure of the network.
/// \param  graph           A DirectedGraph object with the same number of nodes
///                         as there are variables in the domain.
////////////////////////////////////////////////////////////////////////////////
BayesianNetwork::BayesianNetwork(
    const std::vector<Domain>& domains,
    const TrainingDataset& trainingPoints,
    const DirectedGraph& graph
    ) 
    :   DirectedGraph(graph), 
        domains_(domains), 
        trainingPoints_(trainingPoints),
        learningClass_(-1)
{
    // ensure that this graph is acyclic. if it isn't, delete all edges and warn
    // the caller
    if (!this->isAcyclic()) {

        AdjacencyList emptyParentList;

        for (int i = 0; i != this->getNodeCount(); ++i) {

            emptyParentList[i] = vector<int>();
        }

        this->setParentList(emptyParentList);

        cout << "Warning (BayesianNetwork.cpp): The initial graph supplied to "
            "the constructor has one or more cycles. All edges were deleted." <<
            endl;
    }

}   // BayesianNetwork


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Attempts to learn the structure of the Bayesian Network based on the
///         data in the set of training points. 
///
///         Uses the "K2 Algorithm", as described by Cooper, G.F, Herskovits, 
///         E., 1992. "A Bayesian Method for the Induction of Probabilistic 
///         Networks from Data" Machine Learning, Vol. 9, pp. 309-347.
///
/// \param  parentLimit     A user-specified upper limit on the number of 
///                         parents any node may have.
/// \param  searchOrder     The order in which the K2 algorithm should traverse
///                         the nodes in its search for the best network.
///                         must be a permutation of integers 0 to nNodes - 1.
/// \param  learningClass   An optional parameter that enables the caller to 
///                         specify that the network be learned from data that 
///                         has the specified class label. Default behavior is 
///                         that data from all classes should be used.
///
/// \return     True if a network was identified that has a better score than 
///             the previous one.
////////////////////////////////////////////////////////////////////////////////
bool
BayesianNetwork::learnNetworkFromData(
    int parentLimit, const std::vector<int>& searchOrder, int learningClass
    ) 
{
    // check the validity of the provided search order
    vector<int> labels(this->getNodeCount());
    std::iota(labels.begin(), labels.end(), 0);

    if (!is_permutation(labels.begin(), labels.end(), searchOrder.begin())) {

        throw invalid_argument("searchOrder is not a permutation of 0 to n-1.");
    }
    

    // update the learning class member variable
    learningClass_ = learningClass;
    

    // compute the total network score of the current parent list and data set
    double oldNetworkScore(this->computeNetworkScore(this->getParentList()));

    // initialize the new adjacency list
    AdjacencyList newGraph;

    for (int i = 0; i != this->getNodeCount(); ++i) {

        newGraph[i] = vector<int>();
    }


    // process each node according to the specified order
    for (vector<int>::const_iterator i = searchOrder.begin(); 
            i != searchOrder.end(); 
            ++i) {
        
        vector<int> currentParents;
        bool okToProceed(true);
        vector<int> predecessors(searchOrder.cbegin(), i);

        double oldScore(computeNodeScore(*i, currentParents));

        while (okToProceed && currentParents.size() < parentLimit) {

            pair<int, double> bestToAdd = findBestParentToAdd(
                *i, currentParents, predecessors
                );

            vector<int> candidateParents(currentParents);
            candidateParents.push_back(bestToAdd.first);

            double newScore(bestToAdd.second);

            if (newScore > oldScore) {

                oldScore = newScore;
                currentParents = candidateParents;

            } else {

                okToProceed = false;
            }
        }

        // update the new graph to have the new parents
        for (vector<int>::const_iterator j = currentParents.begin();
            j != currentParents.end();
            ++j) {

            bayesnet::addEdge(*i, *j, newGraph);
        }
    }

    // check if the new network is better the the old one. if it is, keep it.
    // if not, keep the old one.
    double newNetworkScore(computeNetworkScore(newGraph));
    bool foundNewNetwork(false);

    if (newNetworkScore > oldNetworkScore) {
        
        this->setParentList(newGraph);
        foundNewNetwork = true;
    }

    // the new graph should never be cyclic by virtue of the algorithm
    #ifdef _DEBUG
        assert(this->isAcyclic());
    #endif
    
    return foundNewNetwork;

}   // learnNetworkFromData
	

// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a score that indicates how well the current network "fits" 
///         the current set of training points.
///
/// \return The score of the network.
////////////////////////////////////////////////////////////////////////////////
double
BayesianNetwork::computeNetworkScore(const AdjacencyList& parentList) const
{
    // initialize the score
    double score(0.0);

    // compute the total score as the product of individual node scores
    AdjacencyList::const_iterator i;

    for (i = parentList.begin(); i != parentList.end(); ++i) {

        score += computeNodeScore(i->first, i->second);
    }
    
    return score;

}   // computeNetworkScore

////////////////////////////////////////////////////////////////////////////////
/// \brief  Finds the best node to add to a child node's list of parents. 
///         
///         Candidates are chosen from the eligible parents, which are the set 
///         difference of the predecessors minus the node's current parents.
///
/// \param  child               Label of the node for which we are finding the 
///                             best parent to add.
/// \param  currentParents      The parents that have already been added to this
///                             child as part of the K2 algorithm.
/// \param  predecessors        The nodes that may preceed the child node in the
///                             network, as defined by the search order.
///
/// \return     The best parent to add to the child node and its resulting node 
///             score.
////////////////////////////////////////////////////////////////////////////////
pair<int, double>
BayesianNetwork::findBestParentToAdd(
    int child, 
    const std::vector<int>& currentParents, 
    const std::vector<int>& predecessors
    ) const
{
    // identify the predecessors that aren't already in the parent list
    unordered_set<int> eligiblePredecessors(
        predecessors.begin(), predecessors.end()
        );
        
    for (vector<int>::const_iterator i = currentParents.begin(); 
            i != currentParents.end(); 
            ++i) {

        eligiblePredecessors.erase(*i);
    }

    // try all eligible predecessors and pick the best one
    int bestToAdd = -1;
    double bestScore(-std::numeric_limits<double>::max());

    for (unordered_set<int>::const_iterator i = eligiblePredecessors.begin(); 
            i != eligiblePredecessors.end(); 
            ++i) {
        
        vector<int> candidateParents = currentParents;
        candidateParents.push_back(*i);

        double currentScore(
            computeNodeScore(child, candidateParents)
            );

        if (currentScore > bestScore) {

            bestScore = currentScore;
            bestToAdd = *i;
        }
    }

    return pair<int, double>(bestToAdd, bestScore);

}   // findBestParentToAdd

////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a score that indicates how well assigning a the given  
///         parents to this node enables it to represent a set of data. 
///
///         See Eq. 12 of [Cooper and Herskovits, 1992].
///
/// \param  child       Label of the node for which we are computing a score.
/// \param  parents     The candidate parents of the node.
///
/// \return     The score of the node if it were to have these parents.
////////////////////////////////////////////////////////////////////////////////
double
BayesianNetwork::computeNodeScore(
    int child, const std::vector<int>& parents
    ) const
{
    // generate vector of class labels that will be used for network learning
    vector<int> learningClasses;

    if (learningClass_ == -1) {
        
        // default. use all classes.
        for (int i = 0; i != trainingPoints_.getNClasses(); ++i) {
            learningClasses.push_back(i);
        }

    } else if (learningClass_ < 0 || 
                learningClass_ > trainingPoints_.getNClasses() - 1
                ) {

        // specified class is out of range. use all classes and print warning.
        for (int i = 0; i != trainingPoints_.getNClasses(); ++i) {
            learningClasses.push_back(i);
        }

        cout << "Warning (DirectedGraph.cpp): User-specified network training "
            "class is out of range.  All classes used." << endl;
        
    } else {

        // use only the specified class
        learningClasses.push_back(learningClass_);
    }


    // from the dataset, get counts of all possible child-node parent-set tuples
    // i.e. the number of times child node has value v and its parents are
    // instantiated as tuple w, for all possible tuples of v and w
    TrainingDataset::FeatureVectorCountMap instanceCounts;
    
    vector<int> childAndParents;
    childAndParents.push_back(child);
    childAndParents.insert(
        childAndParents.end(), parents.begin(), parents.end()
        );   
    
    vector<int> instanceKey(childAndParents.size());

    // for all classes of interest
    for (vector<int>::const_iterator i = learningClasses.begin(); 
            i != learningClasses.end(); 
            ++ i) {
        
        // process each training point with this class label
        typedef TrainingDataset::FeatureVectorCountMap::const_iterator cit;
        cit tpbegin(trainingPoints_.getClassTrainingPointCounts(*i).begin());
        cit tpend(trainingPoints_.getClassTrainingPointCounts(*i).end());

        for (cit j = tpbegin; j != tpend; ++j) {
        
            for (int k = 0; k != childAndParents.size(); ++k) {
                instanceKey[k] = j->first[childAndParents[k]];
            }

            size_t count = j->second;
        
            if (instanceCounts.find(instanceKey) != instanceCounts.end()) {
            
                // this key is already in the map
                instanceCounts[instanceKey] += count;

            } else {
            
                // this key is not yet in the map
                instanceCounts[instanceKey] = count;
            }
        }
    }

    // compute the score. the score is computed differently if this node has no
    // parents.
    double score = 0.0;
    
    if (parents.size() > 0) {

        // generate all of the unique instances of the parent tuples
        vector<int> parentDomainSizes;
        vector<int>::const_iterator i;

        for (i = parents.begin(); i != parents.end(); ++i) {

            parentDomainSizes.push_back(domains_[*i].getSize());
        }

        vector< vector<int> > allParentInstances =
            generatePermutationsWithRepetition(parentDomainSizes);


        // loop through all possible instantiations of the parents and update
        // the running score as a cumulative product
        for (vector< vector<int> >::const_iterator i = 
                allParentInstances.begin(); 
                i != allParentInstances.end(); 
                ++i) {
        
            // compute the subscore: product of factorials of the number of 
            // cases in the dataset in which the child node x has value vj and 
            // it parents are instantiated as the ith parent permutation

            // while we're at it, add the total occurrences in this dataset of 
            // this instantiation of the parents regardless of child node value
            double subScore = 0.0;
            size_t totalInstances = 0;
            int r = domains_[child].getSize();

            for (int j = 0; j != r; ++j) {
            
                vector<int> key;
                key.push_back(j);
                key.insert(key.end(), i->begin(), i->end());

                if (instanceCounts.find(key) != instanceCounts.end()) {

                    subScore += logFactorial(instanceCounts[key]);
                    totalInstances += instanceCounts[key];
                }
            }

            // update the running score
            score += logFactorial(r - 1) - 
                logFactorial(totalInstances + r - 1) + subScore;
        }

    } else {
        
        // this node has no parents. in this case, we don't care what the parent
        // instantiations are because there are no parents. we only care about
        // the number of cases in the dataset in which the "child" node x has
        // each of it's values vj. the score is not a cumulative product but a
        // single "parent agnostic" calculation.

        // compute the subscore: product of factorials of the number of 
        // cases in the dataset in which the child node x has value vj

        // while we're at it, add up all of the cases in the dataset
        double subScore = 0.0;
        int r = domains_[child].getSize();
        size_t totalInstances = 0;

        for (int j = 0; j != r; ++j) {
            
            vector<int> key(1, j);

            if (instanceCounts.find(key) != instanceCounts.end()) {

                subScore += logFactorial(instanceCounts[key]);
                totalInstances += instanceCounts[key];
            }
        }

        // compute the score
        score += logFactorial(r - 1) - 
            logFactorial(totalInstances + r - 1) + subScore;
    }

    return score;

}   // computeNodeScore

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet