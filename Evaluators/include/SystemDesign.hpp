//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDE GUARD
// #############################################################################
#ifndef EVALUATORS_INCLUDE_SYSTEMDESIGN_HPP
#define EVALUATORS_INCLUDE_SYSTEMDESIGN_HPP


// INCLUDES
// #############################################################################
#include "Evaluator.hpp"

#include <string>
#include <utility>
#include <vector>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################

namespace evaluators {


// CLASS DEFINITION
// #############################################################################
class SystemDesign : public Evaluator
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // Structors 
    SystemDesign(const std::string& filename, bool trackingOn = false);


    // Member Functions
    virtual
    std::vector<double>
    evaluateDesign(const std::vector<double>& xValues);

    virtual
    int
    getNObjectives() const;

    virtual
    int
    getNConstraints() const;
        
    
    // PRIVATE DECLARATIONS
    // =================================
private:
    // Typedefs
    typedef std::pair<size_t, size_t> TechOpt;

    typedef std::pair<TechOpt, TechOpt> Dependency;


    // Member Functions
    void
    readDataFile(const std::string& filename);


    // Data Members
    size_t nSubsystems_;

    double costLimit_;

    std::vector< std::vector<double> > costs_;

    std::vector< std::vector<double> > utilities_;

    std::vector<Dependency> necessitations_;

    std::vector<Dependency> obviations_;

    
};  // class SystemDesign


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluator

#endif  // EVALUATORS_INCLUDE_SYSTEMDESIGN_HPP