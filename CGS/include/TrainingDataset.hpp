//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the definition of the TrainingDataset class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_TRAININGDATASET_HPP
#define CGS_INCLUDE_TRAININGDATASET_HPP


// INCLUDES
// #############################################################################
#include "TrainingPoint.hpp"

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A training dataset that may be used to train a Classifier object.
////////////////////////////////////////////////////////////////////////////////
class TrainingDataset
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // NESTED COMPOUNDS

    /// \brief  A functor that computes the hash value of a feature vector.
    struct FeatureHash {
        size_t operator() (const std::vector<int>& x) const;
    };
    

    // TYPPEDEFS
    typedef std::unordered_map<std::vector<int>, size_t, FeatureHash> 
        FeatureVectorCountMap;
    

    // STRUCTORS
    explicit
    TrainingDataset(int nClasses, const std::string& filename = "");
    

    // MEMBER FUNCTIONS
    int
    addPoint(const TrainingPoint& tp);

    int
    removePoint(const TrainingPoint& tp);
    
    int
    movePoint(const std::vector<int>& xValues, int oldC, int newC);

    void
    addPointsFromFile(const std::string& filename);

    bool
    containsPoint(const TrainingPoint& tp) const;

    int
    getNClasses() const;

    const FeatureVectorCountMap&
    getClassTrainingPointCounts(int classLabel) const;

    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // DATA MEMBERS

    /// \brief The collection of classifier training data.
    ///
    ///         The training data is stored as a vector of associative
    ///         containers. The elements of the outer vector correspond to the
    ///         class labels that may exist in the training data. I.e., all 
    ///         training points of a particular class are placed in the same 
    ///         element of the outer vector. Within that vector, class-specific 
    ///         "count maps" contain the counts of each training point, where
    ///         the map key is a specific feature vector and the map value is
    ///         the number of times that specific feature vector with that class
    ///         label appears in the training data set.
    std::vector<FeatureVectorCountMap> trainingData_;
    
    
};  // class TrainingDataset


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet

#endif  // CGS_INCLUDE_TRAININGDATASET_HPP