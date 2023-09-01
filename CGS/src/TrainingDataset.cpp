//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the TrainingDataset class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "TrainingDataset.hpp"

#include "utilities.hpp"

#include <fstream>
#include <iostream>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::ifstream;
using std::map;
using std::string;
using std::unordered_map;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a TrainingDataset object with specified parameters.
///
/// \param  nClasses    The number of classes (categories) in the training data
/// \param  filename    A .csv filename from which to read in the
///                     initial set of training points. If no filename or an
///                     invalid filename is given, no data will be read.
////////////////////////////////////////////////////////////////////////////////
TrainingDataset::TrainingDataset(int nClasses, const string& filename) 
    :   trainingData_(nClasses)
{
    // add points from file to the training data set if one is given
    if (filename != "") {
        
        this->addPointsFromFile(filename);
    }

}   // TrainingDataSet


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds a new training point to the dataset
///
/// \param  tp      A new TrainingPoint object to add to the training set.
///
/// \return     The number of points added. Value of 1 indicates that the point 
///             was successfully added. Value of 0 indicates a failure.
////////////////////////////////////////////////////////////////////////////////
int
TrainingDataset::addPoint(const TrainingPoint& tp)
{
    // check to be sure the class label is in bounds
    if (tp.getClassLabel() >= 0 || tp.getClassLabel() < trainingData_.size()) {
        
        trainingData_[tp.getClassLabel()][tp.getXValues()] += 1;
        return 1;
    }

    return 0;

}   // addPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Removes a training point from the dataset
///
/// \param  tp      A TrainingPoint object to remove from the training set.
///
/// \return     The number of points removed. Value of 1 indicates that the 
///             point was successfully removed. Value of 0 indicates a failure.
////////////////////////////////////////////////////////////////////////////////
int
TrainingDataset::removePoint(const TrainingPoint& tp)
{   
    // get an iterator to the point in the dataset
    FeatureVectorCountMap::iterator i =
        trainingData_[tp.getClassLabel()].find(tp.getXValues());

    if (i == trainingData_[tp.getClassLabel()].end()) {
        
        // attribute vector does not exist in the map, do nothing and return 0
        return 0;

    } else if (i->second == 0) {
    
        // key exists but there are 0 instances, erase it from the map
        trainingData_[tp.getClassLabel()].erase(i);
        return 0;

    } else if (i->second == 1) {
        
        // there is exactly 1 instance, erase it from the map and return 1
        trainingData_[tp.getClassLabel()].erase(i);
        return 1;

    } else {
    
        // there is more than 1 instance, decrement it and return 1
        i->second -= 1;
        return 1;
    }

    return 0;
    
}   // removePoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Reassigns a class label by moving an instance of a specific 
///         attribute vector from one class bin to another. 
///
/// \param  xValues     An attribute vector of the point to reclassify.
/// \param  oldC        The old class of the point.
/// \param  newC        The new class of the point.
///
/// \return     The number of points moved. Value of 1 indicates that the 
///             point was successfully moved. Value of 0 indicates a failure.
////////////////////////////////////////////////////////////////////////////////
int
TrainingDataset::movePoint(const vector<int>& xValues, int oldC, int newC)
{
    // remove point from the old class, and only add it to the new class if the
    // remove was successful
    if (this->removePoint(TrainingPoint(xValues, oldC))) {

        return this->addPoint(TrainingPoint(xValues, newC));
    }

    return 0;

}   // movePoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds points to the data set from a .csv file
///
/// \param  filename    A .csv filename from which to read in the initial set of
///                     training points. If no filename or an invalid filename 
///                     is given, no data will be read.
////////////////////////////////////////////////////////////////////////////////
void
TrainingDataset::addPointsFromFile(const string& filename)
{
    ifstream ifs(filename);
    size_t nPointsRead(0);
    size_t nPointsAdded(0);

    // check to be sure the filename is valid
    if (ifs.is_open()) {
    
        TrainingPoint tp;


        while (tp.readCsv(ifs)) {

            nPointsRead += 1;
            nPointsAdded += this->addPoint(tp);
        }

    } else {

        std::cout << "Warning (TrainingDataset.cpp): Invalid filename. No "
            "points added to dataset." << std::endl;
    }

    // print warning if not all points in the datafile were added
    if (nPointsAdded < nPointsRead) {

        std::cout << "Warning (TrainingDataset.cpp): Not all points in the "
            "file were successully added." << std::endl;
    }

}   // addPointsFromFile


////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if a particular attribute vector / class label pair exists
///         in the training set.
///
/// \param  tp  The attribute vector / class label pair to look for.
///
/// \return     \b True if this specific attribute vector exists in the 
///             training set.
////////////////////////////////////////////////////////////////////////////////
bool
TrainingDataset::containsPoint(const TrainingPoint& tp) const
{
    
    return trainingData_[tp.getClassLabel()].find(tp.getXValues()) != 
            trainingData_[tp.getClassLabel()].end();

}   // containsPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the number of categories (classes) in the training dataset
///
/// \return     The number of classes in the training dataset.
////////////////////////////////////////////////////////////////////////////////
int
TrainingDataset::getNClasses() const
{
    return static_cast<int>(trainingData_.size());

}   // getNClasses


////////////////////////////////////////////////////////////////////////////////
/// \brief  Accesses all of the training points (feature-vector, instance-count
///         pairs) of a specified class.
///
/// \param  classLabel      The class label whose data to return.
///
/// \return     Map of the feature-vector, instance-count pairs.
////////////////////////////////////////////////////////////////////////////////
const TrainingDataset::FeatureVectorCountMap&
TrainingDataset::getClassTrainingPointCounts(int classLabel) const
{
    // check for validity class label argument
    if (classLabel < 0 || classLabel > this->getNClasses() - 1) {

        std::cout << "Warning (TrainingDataset.cpp): Class label out of range. "
            "Data accessed is from class 0." << std::endl;

        classLabel = 0;
    }

    return trainingData_[classLabel];

}   // getClassTrainingPointCounts


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a hash value of a vector of integers. 
///
/// \param  x   The vector of ints to compute a hash value for.
///
/// \return     The hash value.
////////////////////////////////////////////////////////////////////////////////
size_t
TrainingDataset::FeatureHash::operator() (const vector<int>& x) const
{
    return utilities::hash::hash_range(x.begin(), x.end());
    
}   //  XHash::operator()


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet