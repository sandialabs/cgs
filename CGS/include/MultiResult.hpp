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
///  \brief Contains the definition of the MultiResult class and declarations
///         of related non-member functions.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_MULTIRESULT_HPP
#define CGS_INCLUDE_MULTIRESULT_HPP


// INCLUDES
// #############################################################################
#include "KnownPoint.hpp"


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class TrainingPoint;
}


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents a design that has been evaluated in a multi-objective
///         optimization context.
///
///         A MultiResult object contains data that is associated with a
///         design that has been evaluated with a multi-objective optimization
///         evaluator.  In addition, it contains data members and comparison
///         operators for working with a collection of MultiResult objects
///         (such as a domination ranking comparitor).
////////////////////////////////////////////////////////////////////////////////
class MultiResult : public KnownPoint
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    MultiResult();

    MultiResult(
        const bayesnet::TrainingPoint& trainingPoint,
        const std::vector<double>& constraintValues,
        const std::vector<double>& objectiveValues
        );
    
    
    // MEMBER FUNCTIONS
    const std::vector<double>&
    getObjectives() const;

    long
    getDominationRank() const;

    double
    getCrowdingDistance() const;

    void
    setDominationRank(long rank);

    void
    setCrowdingDistance(double crowdingDist);


    // FUNCTORS

    /// \brief  A functor for comparing MultiResult objects by domination rank.
    struct DominationRankCompare
    {
        bool operator() (const MultiResult& lhs, const MultiResult& rhs) const;
    };

    /// \brief  A functor for comparing MultiResult objects by one of their
    ///         objective function values.
    struct ObjectiveCompare
    {
    public:
        
        ObjectiveCompare(int objToComp) {objToComp_ = objToComp;}

        bool operator() (const MultiResult& lhs, const MultiResult& rhs) const;

    private:
        int objToComp_;

    };

    /// \brief  A functor for comparing MultiResult objects by objective-space
    ///         crowding distance.
    struct CrowdingCompare
    {
        bool operator() (const MultiResult& lhs, const MultiResult& rhs) const;
    };
    

    // PRIVATE DECLARATIONS
    // =================================
private:
    // DATA MEMBERS

    /// \brief  The crowding distance of this design.
    ///
    ///         The crowding distance of a design is a scalar metric that
    ///         indicates how crowded it is by other designs of the same
    ///         domination ranking.
    double crowdingDistance_;

    /// \brief  The domination rank (i.e., Pareto "layer") of this design.
    ///
    ///         The domination rank of a design is an indicator of which Pareto
    ///         layer it lies in when compared to a set of designs in 
    ///         multi-objective space.  Given a collection of designs, all
    ///         non-dominated designs are given the best ranking.  If
    ///         all of those designs are removed, the remaining non-dominated
    ///         designs are then given the next-best ranking.  
    long dominationRank_;
    
    /// \brief  The objective function values of this design.
    std::vector<double> objectiveValues_;


};  // class MultiResult


// NON-MEMBER FUNCTIONS
// #############################################################################
bool
isDominatedBy(const MultiResult& mr1, const MultiResult& mr2);


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs

#endif  // CGS_INCLUDE_MULTIRESULT_HPP