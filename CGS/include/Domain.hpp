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
///  \brief Contains the definition of the Domain class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_DOMAIN_HPP
#define CGS_INCLUDE_DOMAIN_HPP


// INCLUDES
// #############################################################################
#include <string>
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Represents the domain of a discrete variable.
///
///         A Domain represents the finite set of values that a discrete 
///         variable may take. The values must be evenly spaced. 
////////////////////////////////////////////////////////////////////////////////
class Domain
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    explicit
    Domain(int size = 2);

    Domain(int size, double min, double max);


    // MEMBER FUNCTIONS
    int
    getSize() const;
    
    int
    xvalToIndex(double x) const;

    double
    indexToXval(int i) const;
    
        
    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    void
    init(int size);


    // DATA MEMBERS
    
    /// \brief  A collection of all values in the domain.
    std::vector<double> domainValues_;

    /// \brief  The maximum value in the domain.
    double max_;

    /// \brief  The minimum value in the domain.
    double min_;

    /// \brief  The size of the domain.
    int size_;


};  // class Domain


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet


#endif  // CGS_INCLUDE_DOMAIN_HPP