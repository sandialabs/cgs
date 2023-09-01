//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the definition of the TrainingPoint class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_TRAININGPOINT_HPP
#define CGS_INCLUDE_TRAININGPOINT_HPP


// INCLUDES
// #############################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A discrete variable attribute vector / class label pair. 
///
///         The attribute vector is expressed as a vector of integers that may 
///         map to other problem specific values. Class labels are also integers
///         and may map to other, more descriptive, value types.
////////////////////////////////////////////////////////////////////////////////

class TrainingPoint
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    TrainingPoint();

    TrainingPoint(const std::vector<int>& xValues, int classLabel);

    // OPERATORS
    bool operator== (const TrainingPoint& tp) const;
    bool operator!= (const TrainingPoint& tp) const {return !operator==(tp);};
    bool operator< (const TrainingPoint& tp) const;
    bool operator> (const TrainingPoint& tp) const {return tp < *this;};
    bool operator<= (const TrainingPoint& tp) const {return !operator> (tp);};
    bool operator>= (const TrainingPoint& tp) const {return !operator< (tp);};


    // MEMBER FUNCTIONS
    std::ifstream&
    readCsv(std::ifstream& ifs);
    
    int
    getClassLabel() const;

    void
    setClassLabel(int classLabel);

    const std::vector<int>&
    getXValues() const;


    // FUNCTORS

    /// \brief  A functor for computing a hash value of a TrainingPoint object.
    struct TrainingPointHash
    {
        size_t operator() (const TrainingPoint& tp) const;
    };


    // PRIVATE DECLARATIONS
    // =================================
private:
    // DATA MEMBERS

    /// \brief  The class label of this classifier training point.
    int classLabel_;

    /// \brief  The feature vector of this classifier training point.
    std::vector<int> xValues_;


};  // class TrainingPoint


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet

#endif  // CGS_INCLUDE_TRAININGPOINT_HPP