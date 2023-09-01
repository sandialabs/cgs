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
///  \brief Contains the definition of the KnownPoint class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_KNOWNPOINT_HPP
#define CGS_INCLUDE_KNOWNPOINT_HPP


// INCLUDES
// #############################################################################
#include "TrainingPoint.hpp"

#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents a design that has been evaluated.
///
///         The KnownPoint class contains data and functional support for a
///         design that has been evaluated with an objective function and
///         constraint evaluator. It is intended to serve as a base class for
///         more specialized types of evaluated designs, such as those that are
///         developed specifically for single- and multi-objective optimization.
////////////////////////////////////////////////////////////////////////////////
class KnownPoint
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    KnownPoint();

    KnownPoint(
        const bayesnet::TrainingPoint& trainingPoint,
        const std::vector<double>& constraintValues
        );

    virtual
    ~KnownPoint();


    // OPERATORS
    bool operator== (const KnownPoint& kp) const;
    bool operator!= (const KnownPoint& kp) const {return !operator==(kp);};
    bool operator< (const KnownPoint& kp) const;
    bool operator> (const KnownPoint& kp) const {return kp < *this;};
    bool operator<= (const KnownPoint& kp) const {return !operator> (kp);};
    bool operator>= (const KnownPoint& kp) const {return !operator< (kp);};
    

    // MEMBER FUNCTIONS
    const bayesnet::TrainingPoint&
    getTrainingPoint() const;

    const std::vector<double>&
    getConstraints() const;

    double
    getOverallViolation() const;

    void
    setOverallViolation(double violation);
    
    bool
    isFeasible() const;


    // FUNCTORS

    /// \brief  A functor for computing a hash value of a KnownPoint object.
    struct Hash
    {
        size_t operator() (const KnownPoint& kp) const;
    };


    // PRIVATE DECLARATIONS
    // =================================
private:
    // DATA MEMBERS

    /// \brief  A classifier training point whose feature vector was evaluated
    ///         to create this KnownPoint.
    bayesnet::TrainingPoint trainingPoint_;
    
    /// \brief  The problem-specific constraint values that result from
    ///         evaluating this design (negative implies a constraint 
    ///         violation).
    std::vector<double> constraintValues_;

    /// \brief  The overal constraint violation of this design, which may be 
    ///         relative to those of a collection of other infeasible designs.
    double overallViolation_;

    /// \brief  Whether this KnownPoint violates any constraints.
    bool isFeasible_;


};  // class KnownPoint


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs


#endif  // CGS_INCLUDE_KNOWNPOINT_HPP