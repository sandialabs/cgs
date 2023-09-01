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
///  \brief Contains the implementation of the MultiResult class and
///         definitions of related non-member functions.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "MultiResult.hpp"

#include "TrainingPoint.hpp"

#include <limits>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using bayesnet::TrainingPoint;

using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace cgs {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a MultiResult with default data values.
////////////////////////////////////////////////////////////////////////////////
MultiResult::MultiResult()
{
}   // MultiResult


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a MultiResult object with prescribed data values.
///
/// \param  trainingPoint       An x-vector / class-label pair representing the
///                             design.
/// \param  constraintValues    The inequality constraint values. Negatives
///                             imply constraint violations.
/// \param  objectiveValues     A vector of objective function values.
////////////////////////////////////////////////////////////////////////////////
MultiResult::MultiResult(
    const bayesnet::TrainingPoint& trainingPoint,
    const std::vector<double>& constraintValues,
    const std::vector<double>& objectiveValues
    )
    :   KnownPoint(trainingPoint, constraintValues),
        dominationRank_(std::numeric_limits<long>::max()),
        objectiveValues_(objectiveValues),
        crowdingDistance_(-std::numeric_limits<double>::max())
{
}   // MultiResult


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the objective function values of this design.
///
/// \return The evaluated objective function values of this design.
////////////////////////////////////////////////////////////////////////////////
const std::vector<double>&
MultiResult::getObjectives() const
{
    return objectiveValues_;

}   // getObjectives


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the domination ranking of this design.
///
/// \return The domination ranking of this design.
////////////////////////////////////////////////////////////////////////////////
long
MultiResult::getDominationRank() const
{
    return dominationRank_;

}   // getDominationRank


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the crowding distance of this design.
///
/// \return The crowding distance of this design.
////////////////////////////////////////////////////////////////////////////////
double
MultiResult::getCrowdingDistance() const
{
    return crowdingDistance_;

}   // getDominationRank


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the domination rank of this design.
///
/// \param  rank    The new domination rank of this design.
////////////////////////////////////////////////////////////////////////////////
void
MultiResult::setDominationRank(long rank)
{
    dominationRank_ = rank;

}   // setDominationRank


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the crowding distance of this design.
///
/// \param  crowdingDistance    The new crowding distance of this design.
////////////////////////////////////////////////////////////////////////////////
void
MultiResult::setCrowdingDistance(double crowdingDistance)
{
    crowdingDistance_ = crowdingDistance;

}   // setCrowdingDistance


////////////////////////////////////////////////////////////////////////////////
/// \brief  Operator for comparing MultiResult objects by domination rank and
///         feasibility.
///
///         Feasible solutions are always preferred to infeasible. If both are
///         feasible, the one with lower domination rank is preferred. If both
///         are infeasible, the one with the lesser overall constraint violation
///         is preferred. Ties are broken arbitrarily by doing a lexicographical 
///         compare on the x values.
///
/// \param  lhs The first MultiResult to compare.
/// \param  rhs The second MultiResult to compare.
///
/// \return \b True if the first argument is preferred to the second.
////////////////////////////////////////////////////////////////////////////////
bool
MultiResult::DominationRankCompare::operator() (
    const MultiResult& lhs, const MultiResult& rhs
    ) const
{
    // returns true if lhs is preferred to rhs
    bool ret(false);

    if (lhs.isFeasible()) {

        if (rhs.isFeasible()) {

            // both are feasible
            if (lhs.getDominationRank() != rhs.getDominationRank()) {

                ret = lhs.getDominationRank() < rhs.getDominationRank();

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

}   // DominationRankCompare::operator()


////////////////////////////////////////////////////////////////////////////////
/// \brief  Operator for comparing two MultiResult objects by a specified 
///         objective value
///         
///         The index of the objective function to compare on is specified 
///         during construction of this functor.  The MultiResult with lower 
///         objective function value is preferred.  Feasibility is ignored.
///
/// \param  lhs     The first MultiResult to compare.
/// \param  rhs     The second  MultiResult to compare.
///
/// \return \b True if the first argument is preferred to the second.
////////////////////////////////////////////////////////////////////////////////
bool
MultiResult::ObjectiveCompare::operator() (
    const MultiResult& lhs, const MultiResult& rhs
    ) const
{
    return lhs.getObjectives()[objToComp_] < rhs.getObjectives()[objToComp_];

}   // ObjectiveCompare::operator()


////////////////////////////////////////////////////////////////////////////////
/// \brief  Operator for comparing two MultiResult objects by crowding
///         distance.
///
///         The object with higher crowding distance is preferred. 
///         Feasibility is ignored.
///
/// \param  lhs     The first MultiResult to compare.
/// \param  rhs     The second  MultiResult to compare.
///
/// \return \b True if the first argument is preferred to the second.
////////////////////////////////////////////////////////////////////////////////
bool
MultiResult::CrowdingCompare::operator() (
    const MultiResult& lhs, const MultiResult& rhs
    ) const
{
    return lhs.getCrowdingDistance() > rhs.getCrowdingDistance();

}   // CrowdingCompare::operator()


// NON-MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Indicates whether one MultiResult object is dominated by another.
///
///         This function assumes that all objectives are sought to be
///         minimized.
///
/// \param  mr1     The first MultiResult object to compare.
/// \param  mr2     The second MultiResult object to compare.
///
/// \return \b True if the first argument is dominated by the second. \b False
///         otherwise.
////////////////////////////////////////////////////////////////////////////////
bool
isDominatedBy(const MultiResult& mr1, const MultiResult& mr2)
{
    size_t nBetterThan(0);
    size_t nBetterOrEqual(0);
    size_t nObj = mr1.getObjectives().size();

    for (size_t i = 0; i != nObj; ++i) {
        
        if (mr2.getObjectives()[i] < mr1.getObjectives()[i]) {
            nBetterThan += 1;
        }

        if (mr2.getObjectives()[i] <= mr1.getObjectives()[i]) {
            nBetterOrEqual += 1;
        }

    }

    if (nBetterOrEqual == nObj && nBetterThan > 0) {

        return true;
    }

    return false;
}

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace cgs