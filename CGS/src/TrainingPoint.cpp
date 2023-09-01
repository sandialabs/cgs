//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the TrainingPoint class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "TrainingPoint.hpp"

#include "utilities.hpp"


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::ifstream;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Default constructor. Instantiates a TrainingPoint object with an 
///         empty feature vector and class label of zero.
////////////////////////////////////////////////////////////////////////////////
TrainingPoint::TrainingPoint() : classLabel_(0) 
{
};  // TrainingPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a TrainingPoint object with specified values.
///
/// \param  xValues     The integers representing the feature vector states.
/// \param  classLabel  The categorical class label of this design.
////////////////////////////////////////////////////////////////////////////////
TrainingPoint::TrainingPoint(const std::vector<int>& xValues, int classLabel)
    :   xValues_(xValues), classLabel_(classLabel) 
{ 
};  //TrainingPoint


// OPERATORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Determines if a given TrainingPoint is equal to *this. Returns true 
///         if the feature vector and the class label are BOTH equal.
///
/// \param  tp      A TrainingPoint to compare to *this.
///
/// \return         \b True if the two TrainingPoints are equal.
////////////////////////////////////////////////////////////////////////////////
bool
TrainingPoint::operator== (const TrainingPoint& tp) const
{
    bool xIsEqual = (this->xValues_ == tp.xValues_);
    bool cIsEqual = (this->classLabel_ == tp.classLabel_);

    return (xIsEqual && cIsEqual);

}   // operator==


////////////////////////////////////////////////////////////////////////////////
/// \brief  Determines if a given TrainingPoint is less than *this. Returns 
///         true if the feature vector of *this is less. If the two 
///         feature vectors are equal, returns true if the class label
///         of *this is less.
///
/// \param  tp      A TrainingPoint to compare to *this.
///
/// \return         \b True of the given TrainingPoint is less than *this.
////////////////////////////////////////////////////////////////////////////////
bool
TrainingPoint::operator< (const TrainingPoint& tp) const
{
    bool xIsLess = (this->xValues_ < tp.xValues_);
    bool cIsLess = (this->classLabel_ < tp.classLabel_);
    bool xIsEqual = (this->xValues_ == tp.xValues_);

    bool ret = false;

    if (xIsLess)
        ret = true;
    else if (xIsEqual && cIsLess)
        ret = true;

    return ret;

}   // operator<

/*
================================================================================
Name:       operators (!=), (>), (<=), (>=)

Function:   NOTE: Defined in header file in terms of operator== and operator <.

Params:     ptp.....a point pair to compare to *this

Returns:    (bool)
================================================================================
*/

/*
================================================================================
Name:       TrainingPointHash::operator()

Function:   Computes a hash value of a TrainingPoint object. Used when storing a 
            set of TrainingPoints in an unordered_map<TrainingPoint, int>.

Params:     ptp.....a point pair to compute the hash value of.

Returns:    (size_t)    The hash value.
================================================================================
*/
////////////////////////////////////////////////////////////////////////////////
/// \brief
////////////////////////////////////////////////////////////////////////////////
size_t
TrainingPoint::TrainingPointHash::operator() (const TrainingPoint& ptp) const
{
    size_t ret = 
        utilities::hash::hash_range(ptp.xValues_.begin(), ptp.xValues_.end());

    return utilities::hash::hash_combine(ret, ptp.classLabel_);

}   // TrainingPointHash::operator()


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Reads a row of a .csv file that contains integers and updates the
///         feature vector and class label of this TrainingPoint. The last
///         integer in the row is assigned to the class label.
///
/// \param  ifs     The input file stream that is associated with the file to be
///                 read from.
///
/// \return         The ifstream that is being read from. It will be in an 
///                 error state if eof is reached or if something was read
///                 that could not be converted to an int.
////////////////////////////////////////////////////////////////////////////////
ifstream&
TrainingPoint::readCsv(ifstream& ifs)
{
    // check if ifs is in an error state
    if (ifs) {
        // clear current xValues_
        xValues_.clear();
        vector<int> v;

        // read ints from a .csv file
        if (utilities::readwrite::readCsvRow(ifs, v)) {

            vector<int>::const_iterator i;
            for (i = v.begin(); i != v.end(); ++i) {

                // if not at the class label, add int to x values. 
                if (i != v.end() - 1) {

                    xValues_.push_back(*i);

                } else {

                    classLabel_ = *i;
                }
            }
        }
    }
    // ifs will be in an error state if end of file reached or if something 
    // other than an integer was in the file

    return ifs;

}   // readCsv


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the class label of thie TrainingPoint.
///
/// \return     The class label.
////////////////////////////////////////////////////////////////////////////////
int
TrainingPoint::getClassLabel() const
{
    return classLabel_;

}   // getClassLabel


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the class label of this Training Point.
///
/// \param  classLabel  The new value of the class label.
////////////////////////////////////////////////////////////////////////////////
void
TrainingPoint::setClassLabel(int classLabel)
{
    classLabel_ = classLabel;

}   // setClassLabel


////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the feature vector of this TrainingPoint.
///
/// \return The feature vector.
////////////////////////////////////////////////////////////////////////////////
const vector<int>&
TrainingPoint::getXValues() const
{
    return xValues_;

}   // getXValues


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet