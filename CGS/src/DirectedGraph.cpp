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
///  \brief Contains the implementation of the DirectedGraph class and
///         definitions of related non-member functions.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "DirectedGraph.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::string;
using std::unordered_set;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a DirectedGraph object.
///
/// \param  nNodes      The number of nodes in this graph.
/// \param  parentList  A valid parent-based adjacency list. Each node (numbered
///                     0, 1, ... , \a n - 1) must be represented.
////////////////////////////////////////////////////////////////////////////////
DirectedGraph::DirectedGraph(
    int nNodes, const DirectedGraph::AdjacencyList& parentList
    )
    :   topologicalOrder_(nNodes), parentList_(parentList), nNodes_(nNodes)
{
    // update the topological order
    this->updateTopologicalOrder();

}   // DirectedGraph


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a DirectedGraph object.
///
/// \param  nNodes              The number of nodes in this graph.
/// \param  edgeListFilename    The name of a file that contains 
///                             an edge list. Each row in the file represents 
///                             one edge, where the first integer is the parent 
///                             node and the second integer is the child node
///                             (numbered 0, 1, ... , \a n - 1).
////////////////////////////////////////////////////////////////////////////////
DirectedGraph::DirectedGraph(int nNodes, const std::string& edgeListFilename)
    :   topologicalOrder_(nNodes), nNodes_(nNodes)
{
    // add all nodes to the parent list
    for (int i = 0; i != nNodes; ++i) {

        parentList_[i] = vector<int>();
    }

    // read edges from the edgefile if one is given
    if (edgeListFilename != "") {

        ifstream ifs(edgeListFilename);
        string line;

        if (ifs.is_open()) {
            
            // warning flag
            bool hasNodesOutOfRange(false);

            // read each line of the edge list file
            while (std::getline(ifs, line)) {

                std::istringstream in(line);
                vector<int> temp;

                std::copy(
                    std::istream_iterator<int>(in),
                    std::istream_iterator<int>(),
                    std::back_inserter(temp)
                    );

                // check that the labels are in the range [0, n-1]
                if (temp[0] > nNodes - 1 || temp[1] > nNodes - 1 || 
                    temp[0] < 0 || temp[1] < 0
                    ) {

                    // don't add the edge and set warning flag to true
                    hasNodesOutOfRange = true;

                } else {

                    parentList_[temp[1]].push_back(temp[0]);
                }
            }

            // print warning if file contains out-of-range node labels
            if (hasNodesOutOfRange) {
                cout << "Warning (DirectedGraph.cpp): Edge list file "
                "contains labels outside the range [0, n-1]. Some "
                "edges not added to graph." << endl;
            }
            
        } else {
            
            cout << "Warning (DirectedGraph.cpp): Invalid filename. No edges "
                "added to parent list." << endl;
        }
    }

    // update the topological order
    this->updateTopologicalOrder();
    
}   // DirectedGraph


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the parents of a specified node.
///
/// \param  node    The node whose parents to return.
///
/// \return         The parents of the input node.
////////////////////////////////////////////////////////////////////////////////
const vector<int>&
DirectedGraph::getParents(int node) const
{
    return parentList_.at(node);

}   // getParents


////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the topological order of this graph if it is cyclic.
///         Otherwise, it returns an empty vector.
///
/// \return The topological order of the graph. Returns an empty set if the 
///         graph is cyclic.
////////////////////////////////////////////////////////////////////////////////
const std::vector<int>
DirectedGraph::getTopologicalOrder() const
{
    
    if (!isAcyclic_) {

        return vector<int>();
    }

    return topologicalOrder_;

}   // getTopologicalOrder


////////////////////////////////////////////////////////////////////////////////
/// \brief  Writes the parent list to the console window.
////////////////////////////////////////////////////////////////////////////////
void
DirectedGraph::printParentList() const
{
    typedef std::map< int, vector<int> > OrderedAdjacencyList;
    OrderedAdjacencyList orderedParentList(
        parentList_.begin(), parentList_.end()
        );

    OrderedAdjacencyList::const_iterator i;
    vector<int>::const_iterator j;
    
    for (i = orderedParentList.begin(); i != orderedParentList.end(); ++i) {
        
        cout << i->first << " <- ";

        for (j = i->second.begin(); j != i->second.end(); ++j) {
            cout << *j << " ";
        }

        cout << endl;
    }

}   // printParentList


////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether this directed graph has cycles. Should be checked
///         prior to getting the topological order since the topological order
///         has no meaning for cyclic directd graphs.
///
/// \return \b True if this graph has no directd cycles.
////////////////////////////////////////////////////////////////////////////////
bool
DirectedGraph::isAcyclic() const
{
    return isAcyclic_;

}   // isAcyclic


// PROTECTED MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the parent list of this graph.
///
/// \return The parent list of this graph.
////////////////////////////////////////////////////////////////////////////////
const DirectedGraph::AdjacencyList&
DirectedGraph::getParentList() const
{
    return parentList_;

}   // getParentList


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the parent list to a new AdjacencyList and updates the 
///         topological order.
///
/// \param  newParentList   The new parent list.
////////////////////////////////////////////////////////////////////////////////
void
DirectedGraph::setParentList(const AdjacencyList& newParentList)
{
    parentList_ = newParentList;

    this->updateTopologicalOrder();

}   // setParentList


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the number of nodes in this graph.
///
/// \return The number of nodes in this graph.
////////////////////////////////////////////////////////////////////////////////
int
DirectedGraph::getNodeCount() const
{
    return nNodes_;

}   //  nNodes


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the topological order of the graph and checks for cycles.
///         Updates the \link isAcyclic_ \endlink member variable accordingly.
////////////////////////////////////////////////////////////////////////////////
void
DirectedGraph::updateTopologicalOrder()
{
    // initializations
    unordered_set<int> exploredNodes;
    unordered_set<int> currentNodes;
    int currentLabel(0);
    bool foundCycle = false;

    // perform depth first search until all nodes have been explored
    int nNodes = static_cast<int>(parentList_.size());

    for (int startNode = 0; startNode != nNodes; ++startNode) {
        
        // clear the set of nodes visited on this recursion tree
        currentNodes.clear();

        // proceed if we haven't already explored this node
        if (exploredNodes.find(startNode) == exploredNodes.end()) {
            
            // call the recursive depth first search algorithm
            foundCycle = sortTopologicalDfs(
                startNode, currentLabel, exploredNodes, currentNodes
                );

            // break out of the loop if a cycle was found
            if (foundCycle) {
                break;
            }
        }
    }

    isAcyclic_ = !foundCycle;

    // check that updated ancestral order is permutation of all variable labels
#ifdef _DEBUG
    
    // topological ordering is only valid if the graph is acyclic. if it is,
    // assert that the ordering is a permutation of integers 0 to n-1.
    if (isAcyclic_) {

        vector<int> allLabels(this->getNodeCount());

        std::iota(allLabels.begin(), allLabels.end(), 0);

        assert(       
            std::is_permutation(
                topologicalOrder_.begin(), 
                topologicalOrder_.end(), 
                allLabels.begin()
            )
        );
    }
#endif

}   // updateTopologicalOrder


////////////////////////////////////////////////////////////////////////////////
/// \brief  Performs topological sort using depth first search (DFS). 
///
///         Since the structure of this graph is defined by a parent-based 
///         adjacency list, this function does the DFS backwards because the 
///         parents (incoming edges) are readily obtainable, while the outgoing 
///         edges are not. This function is recursive.
///
/// \param  startNode       The node that this recursion starts with.
/// \param  currentLabel    The current label for the topological order vector.
/// \param  exploredNodes   Set of all previously explored nodes.
/// \param  currentNodes    The set of nodes explored in the current branch of
///                         the recursion tree.
///
/// \return                 \b True if a cycle is detected.
////////////////////////////////////////////////////////////////////////////////
bool
DirectedGraph::sortTopologicalDfs(
    int startNode,
    int& currentLabel,
    std::unordered_set<int>& exploredNodes,
    std::unordered_set<int>& currentNodes
    )
{
    // mark this node as explored
    exploredNodes.insert(startNode);
    currentNodes.insert(startNode);

    // recursively call DFS along all of this node's edges
    for (vector<int>::iterator i = parentList_.at(startNode).begin(); 
            i != parentList_.at(startNode).end(); 
            ++i) {

        // check if this node has already been explored in this dive of DFS. 
        // if it has, this graph has a cycle
        if (currentNodes.find(*i) != currentNodes.end()) {

            return true;
        }

        // check if this node has been explored
        if (exploredNodes.find(*i) == exploredNodes.end()) {
            
            // perform dfs on this node. 
            // true is returned if a cycle was detected.
            if (sortTopologicalDfs(
                    *i, 
                    currentLabel, 
                    exploredNodes, 
                    currentNodes
                    )) {
                
                return true;
            }
        }
        // clear current nodes before exploring next branch up the parent list
        currentNodes.clear();
    }
    
    // update the topological order
    topologicalOrder_[currentLabel] = startNode;
    ++currentLabel;

    // if we got here, return false to indicate that no cycle was detected
    return false;

}   // sortTopologicalDfs


// NON-MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if a given node has parents in the given graph.
///
/// \param  node    The label of the node who is being checked for parents.
/// \param  graph   A DirectedGraph that this node is a member of.
///
/// \return         \b True if the node has parents in the graph.
////////////////////////////////////////////////////////////////////////////////
bool
hasParents(int node, const DirectedGraph& graph)
{
    return graph.getParents(node).size() != 0;

}   // hasParents


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds a directed edge to an AdjacencyList
///
/// \param  fromNode        The parent of the edge being added.
/// \param  toNode          The child of the edge being added.
/// \param  parentList      A reference to an adjacency list that the edge is
///                         being added to.
////////////////////////////////////////////////////////////////////////////////
void
addEdge(int toNode, int fromNode, DirectedGraph::AdjacencyList& parentList)
{
    parentList[toNode].push_back(fromNode);
   
}   // addEdge

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet