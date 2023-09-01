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
///  \brief Contains the implementation of the KnownPoint class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "KnownPoint.hpp"
#include "utilities.hpp"

#include <algorithm>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::any_of;
using std::vector;

using bayesnet::TrainingPoint;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Default constructor.
///
///         Instantiates a KnownPoint object with no design variable
///         information, no constraint value information, zero overall
///         violation, and a feasibility state of \b true.
////////////////////////////////////////////////////////////////////////////////
KnownPoint::KnownPoint()
    :   overallViolation_(0.0), isFeasible_(true)
{
}   // KnownPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a KnownPoint object with known data attributes.
///
///         The overall violation is set to zero regardless of constraint
///         values. To manipulate the overall violation of a KnownPoint, use the
///         \link setOverallViolation \endlink member function.
///
/// \param  trainingPoint       An x-vector / class label pair representing the
///                             evaluated design.
/// \param  constraintValues    The inequality constraint values. Negatives
///                             imply constraint violations.
////////////////////////////////////////////////////////////////////////////////
KnownPoint::KnownPoint(
    const TrainingPoint& trainingPoint,
    const std::vector<double>& constraintValues
    )
    :   trainingPoint_(trainingPoint), 
        constraintValues_(constraintValues),
        overallViolation_(0.0)
{
    // assign correct value to isFeasible base on constraintValues
    isFeasible_ = !any_of(
        constraintValues.begin(), 
        constraintValues.end(), 
        [](double i){return i < 0;}
        );

}   // KnownPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Default destructor.
////////////////////////////////////////////////////////////////////////////////
KnownPoint::~KnownPoint()
{
}   //  ~KnownPoint


// OPERATORS
// #############################################################################

////////////////////////////////////////////////////////////////////////////////
/// \brief  Determines if the xValues of the given KnownPoint are equal to those
///         of *this. 
///
/// \param  kp      A point pair to compare to *this.
///
/// \return \b True if the two KnownPoint's xValues are equal.
////////////////////////////////////////////////////////////////////////////////
bool
KnownPoint::operator== (const KnownPoint& kp) const
{
    return this->trainingPoint_.getXValues() == 
                kp.getTrainingPoint().getXValues();

}   // operator==


////////////////////////////////////////////////////////////////////////////////
/// \brief  Determines if a given TrainingPoint is less than *this. Returns 
///         true if the xValues_ data member of *this is less. 
///
/// \param  kp      A KnownPoint to compare to *this.
///
/// \return \b True of the given TrainingPoint is less than *this.
////////////////////////////////////////////////////////////////////////////////
bool
KnownPoint::operator< (const KnownPoint& kp) const
{
    return this->trainingPoint_.getXValues() < 
                kp.getTrainingPoint().getXValues();
    
}   // operator<


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the attribute vector / class label pair of this KnownPoint.
///
/// \return A TrainingPoint that is contained in this KnownPoint.
////////////////////////////////////////////////////////////////////////////////
const TrainingPoint&
KnownPoint::getTrainingPoint() const
{
    return trainingPoint_;

}   // getTrainingPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the constraint values of this KnownPoint.
///
/// \return A vector of raw constraint values contained in this KnownPoint.
////////////////////////////////////////////////////////////////////////////////
const std::vector<double>&
KnownPoint::getConstraints() const
{
    return constraintValues_;

}   // getConstraints


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the amount of the infeasibility for this design.  
///
///         Overall violation may be relative to collection of other infeasible 
///         designs and is therefore calculated outside of this class. 
///         It is intended to be a scalar value that represents the total extent
///         to which a design is infeasible.
///
/// \return The overal violation of this design. Returns a value of zero if the 
///         design is feasible.
////////////////////////////////////////////////////////////////////////////////
double
KnownPoint::getOverallViolation() const
{
    return overallViolation_;

}   // getOverallViolation


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the value of the overall violation of this KnownPoint.
///
/// \param  violation   The new overall constraint violation value.
////////////////////////////////////////////////////////////////////////////////
void
KnownPoint::setOverallViolation(double violation)
{
    overallViolation_ = violation;

}   // setOverallViolation


////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether this KnownPoint has any constraint violations.
///
///         Returns true only if ALL inequality constraints are non-negative.
///
/// \return \b True if the point is feasible.
////////////////////////////////////////////////////////////////////////////////
bool
KnownPoint::isFeasible() const
{
    return isFeasible_;

}   // isFeasible


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes a hash value for a KnownPoint object.
///
/// \param  kp  The known point for which to compute a hash value.
///
/// \return     The hash value of the given KnownPoint.
////////////////////////////////////////////////////////////////////////////////
size_t
KnownPoint::Hash::operator() (const KnownPoint& kp) const
{
    vector<int>::const_iterator first(
        kp.getTrainingPoint().getXValues().begin()
        );

    vector<int>::const_iterator last(
        kp.getTrainingPoint().getXValues().end()
        );

    return utilities::hash::hash_range(first, last);
    
}   //  KnownPoint::Hash::operator()


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs