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
///  \brief Contains the definition of the DirectedGraph class and declarations 
///         of related non-member functions.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_DIRECTEDGRAPH_HPP
#define CGS_INCLUDE_DIRECTEDGRAPH_HPP


// INCLUDES
// #############################################################################
#include <string>
#include <unordered_map>
#include <unordered_set>
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
/// \brief  Represents a graph of nodes and directed edges.
////////////////////////////////////////////////////////////////////////////////
class DirectedGraph
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // TYPEDEFS

    /// \brief  An unordered adjacency list. There must be one unique key per
    ///         node in the graph (numbered 0, 1, ..., \a n - 1). The mapped
    ///         values are vectors of parent nodes to each key. A node key with
    ///         an empty vector mapped to it has no parents.
    typedef std::unordered_map< int, std::vector<int> > AdjacencyList;
    

    // STRUCTORS
    DirectedGraph(int nNodes, const AdjacencyList& parentList);

    explicit
    DirectedGraph(int nNodes, const std::string& edgeListFilename = "");

        
    // MEMBER FUNCTIONS
    const std::vector<int>&
    getParents(int node) const;

    const std::vector<int>
    getTopologicalOrder() const;

    void
    printParentList() const;

    bool
    isAcyclic() const;
    

    // PROTECTED DECLARATIONS
    // =================================
protected:
    // MEMBER FUNCTIONS
    const AdjacencyList&
    getParentList() const;

    void
    setParentList(const AdjacencyList& newParentList);

    int
    getNodeCount() const;


    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    void
    updateTopologicalOrder();
    
    bool
    sortTopologicalDfs(
        int startNode, 
        int& currentLabel, 
        std::unordered_set<int>& exploredNodes,
        std::unordered_set<int>& currentNodes
        );


    // DATA MEMBERS

    /// \brief  An unordered mapping of all nodes and their parents, i.e. the
    ///         nodes that have directed edges towards the keyed node.
    AdjacencyList parentList_;

    /// \brief  The topological order of the graph, i.e. a valid ordering for 
    ///         which for every directed edge \a uv from node \a u to node \a v,
    ///         \a u comes before \a v in the topological ordering. If \link
    ///         isAcyclic_ \endlink is \b false, the topological
    ///         order has no meaning.
    std::vector<int> topologicalOrder_;

    /// \brief  Whether this graph contains any cycles.
    bool isAcyclic_;

    /// \brief  The number of nodes (vertices) in this graph.
    const int nNodes_;
    
    
};  // class DirectedGraph


// NON-MEMBER FUNCTIONS
// #############################################################################
bool
hasParents(int node, const DirectedGraph& graph);

void
addEdge(int toNode, int fromNode, DirectedGraph::AdjacencyList& parentList);

    
// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet

#endif  // CGS_INCLUDE_DIRECTEDGRAPH_HPP