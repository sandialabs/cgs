//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the implementation of the SingleResult class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "SingleResult.hpp"

#include "TrainingPoint.hpp"


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::vector;

using bayesnet::TrainingPoint;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Constructs a SingleResult object with default values.
////////////////////////////////////////////////////////////////////////////////

SingleResult::SingleResult()
{
}   // SingleResult


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a SingleResult object with specified parameters.
///
/// \param  trainingPoint       An x vector / class label pair representing the
///                             design.
/// \param  constraintValues    The inequality constraint values. Negatives
///                             imply constraint violations.
/// \param  objectiveValue      The objective function value.
////////////////////////////////////////////////////////////////////////////////
SingleResult::SingleResult(
    const bayesnet::TrainingPoint& trainingPoint,
    const std::vector<double>& constraintValues,
    double objectiveValue
    )
    :   KnownPoint(trainingPoint, constraintValues),
        objectiveValue_(objectiveValue)
{
}   // SingleResult


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the objective function value of this design.
///
/// \return The objective function value of this evaluated design.
////////////////////////////////////////////////////////////////////////////////
double
SingleResult::getObjective() const
{
    return objectiveValue_;

}   // getObjective


////////////////////////////////////////////////////////////////////////////////
/// \brief  Comparison operator for two SingleResult objects.
///
///         Feasible solutions are always preferred to infeasible. If both are
///         feasible, the one with lower objective value is preferred. If both
///         are infeasible, the one with the lesser overall constraint violation
///         is preferred. Ties are broken arbitrarily by doing a lexicographical 
///         compare on the xValues. 
///
/// \param  lhs     The first SingleResult to compare.
/// \param  rhs     The second SingleResult to compare.
///
/// \return \b True if the first argument is preferred to the second.
////////////////////////////////////////////////////////////////////////////////
bool
SingleResult::ObjectiveCompare::operator() 
    (const SingleResult& lhs, const SingleResult& rhs) const
{
    // initialize return value
    bool ret(false);

    if (lhs.isFeasible()) {

        if (rhs.isFeasible()) {

            // both are feasible
            if (lhs.getObjective() != rhs.getObjective()) {

                ret = lhs.getObjective() < rhs.getObjective();

            } else {

                ret = lhs.getTrainingPoint().getXValues() < 
                        rhs.getTrainingPoint().getXValues();
            }

        } else {

            // only lhs is feasible
            ret = true;

        }

    } else {

        if (rhs.isFeasible()) {

            // only rhs is feasible
            ret = false;

        } else {

            // both are infeasible
            if (lhs.getOverallViolation() != rhs.getOverallViolation()) {

                ret = lhs.getOverallViolation() < rhs.getOverallViolation();

            } else {

                ret = lhs.getTrainingPoint().getXValues() < 
                        rhs.getTrainingPoint().getXValues();
            }
        }
    }
    
    return ret;

}   // ObjectiveCompare::operator()


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs