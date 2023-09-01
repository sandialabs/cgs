//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

// FILE CONTENTS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Contains the implementation of the Domain class.
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "Domain.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::invalid_argument;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a domain with integer values from 0 to (\b size - 1). 
///         
///         If \b size is < 2, an invalid_argument exception is thrown and 
///         handled by using the default constructor argument.
///
/// \param  size    The number of values in the domain.
////////////////////////////////////////////////////////////////////////////////
Domain::Domain(int size) : domainValues_(size)
{
    try {
        if (size < 2) {

            throw invalid_argument("Warning (Domain.cpp): Invalid constructor "
                "argument. Default constructor used.");
        }

        this->init(size);

    } catch (invalid_argument e) {
        std::cerr << e.what() << std::endl;
        domainValues_.resize(2);
        this->init(2);
    }
    
}   // Domain

////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a domain with \b size evenly spaced values beginning at 
///         min and ending at max. 
///
///         If \b size is < 2, an invalid_argument exception is thrown and 
///         handled by using the default constructor argument.
///
/// \param  size    The number of values in the domain.
/// \param  min     The minimum value in the domain.
/// \param  max     The maximum value in the domain.
////////////////////////////////////////////////////////////////////////////////
Domain::Domain(int size, double min, double max) : domainValues_(size)
{
    try {
        if (size < 2 || min >= max) {

            throw invalid_argument("Warning (Domain.cpp): Invalid constructor "
                "arguments. Default constructor used.");
        }

        size_ = size;
        min_ = min;
        max_ = max;

        double xVal = min;
        double delta = (max - min) / static_cast<double>(size - 1);

        for (size_t i = 0; i != size; ++i) {
            domainValues_[i] = xVal;
            xVal += delta;
        }

    } catch (invalid_argument e) {
        std::cerr << e.what() << std::endl;
        init(2);
    }
    
}   // Domain


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the number values in the domain.
///
/// \return The size of the domain.
////////////////////////////////////////////////////////////////////////////////
int
Domain::getSize() const
{
    return size_;

}   // getSize

////////////////////////////////////////////////////////////////////////////////
/// \brief  Converts a real value to an index representing the bin in which that 
///         value belongs. 
///
///         This function assumes domain values are evenly spaced. Bin size is 
///         equal to the difference between two adjecent real domain values. 
///         Domain values lie at the median value of each bin.
///
/// \param  x   The real value to convert to an index.
///
/// \return     The index of the bin in which the real value x belongs.
////////////////////////////////////////////////////////////////////////////////
int
Domain::xvalToIndex(double x) const
{
    double binSz = (max_ - min_) / static_cast<double>(size_ - 1);
    
    double i = ((x - min_ + binSz/2) * size_) / (max_ - min_ + binSz);

    if (i < 0)
        i = 0;

    if (i > size_ - 1)
        i = static_cast<double>(size_ - 1);

    return static_cast<int>(i);

}   // indexFromXval

////////////////////////////////////////////////////////////////////////////////
/// \brief      Converts an index i to the i'th domain's x value.
///
/// \param  i   The index to convert to an x value.
///
/// \return     The i'th domain's real value.
////////////////////////////////////////////////////////////////////////////////
double
Domain::indexToXval(int i) const
{
    return domainValues_[i];

}   // indexToXval


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Initializes a domain to to consist of \b size evenly spaced whole 
///         numbers beginning at 0 and ending at (\b size - 1).
///
/// \param  size    The size of the domain.
////////////////////////////////////////////////////////////////////////////////
void
Domain::init(int size)
{
    size_ = size;
    min_ = 0.0;
    max_ = static_cast<double>(size - 1);

    for (int i = 0; i != size; ++i) {
        domainValues_[i] = static_cast<double>(i);
    }

}   // init


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet