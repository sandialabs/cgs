//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "SystemDesign.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::getline;
using std::ifstream;
using std::istringstream;
using std::replace;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       SystemDesign

Function:   Instantiates a SystemDesign evaluator object.

Params:     filename....the name of a .csv file with the problem parameters.
            trackingOn..default = false. sets whether to track rate of
                        convergence
================================================================================
*/
SystemDesign::SystemDesign(const string& filename, bool trackingOn)
    :   Evaluator(trackingOn)
{
    // read the data file and populate data members
    this->readDataFile(filename);
        
}   // SystemDesign


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       evaluateDesign

Function:   Evaluates the objective and constraints of a design vector.

Params:     xValues.....the design to evaluate

Returns:    (vector<double>)  the objective and constraint values
================================================================================
*/
vector<double>
SystemDesign::evaluateDesign(const vector<double>& xValues)
{
    // convert the input design vector to an integral type
    vector<int> x(xValues.size());
    size_t subSysNo(0);

    for (vector<double>::const_iterator i = xValues.begin(); 
            i != xValues.end(); 
            ++i) {

        x[subSysNo++] = static_cast<int>(*i);
    }

    // compute the total cost and utility of this design
    double totalCost(0.0);
    double totalUtility(0.0);

    for (size_t i = 0; i != nSubsystems_; ++i) {
        
        totalCost += costs_[i][x[i]];
        totalUtility += utilities_[i][x[i]];
    }

    // compute the number of necessitation violations
    double nNecessViolations(0.0);

    for (vector<Dependency>::const_iterator i = necessitations_.begin(); 
            i != necessitations_.end(); 
            ++i) {

        TechOpt driver(i->first);
        TechOpt target(i->second);

        if (x[driver.first] == driver.second && 
            x[target.first] != target.second) {

            nNecessViolations += 1.0;
        }
    }

    // compute the number of obviation violations
    double nObviationViolations(0.0);

    for (vector<Dependency>::const_iterator i = obviations_.begin();
            i != obviations_.end();
            ++i) {

        TechOpt driver(i->first);
        TechOpt target(i->second);

        if (x[driver.first] == driver.second &&
            x[target.first] == target.second) {

            nObviationViolations += 1.0;
        }
    }

    // compute the objective value
    double f(totalUtility * (-1));

    // constraint 0: cost
    double g0 = costLimit_ - totalCost;

    // constraints 1 and 2: necessitation and obvation
    double g1 = 0.0 - nNecessViolations;
    double g2 = 0.0 - nObviationViolations;

    // return vector of objective and constraint values [f, g0, g1, ...]
    vector<double> ret;
    ret.push_back(f);
    ret.push_back(g0);
    ret.push_back(g1);
    ret.push_back(g2);

    // track rate of convergence
    if (trackingOn_) {

        this->trackRateOfConvergence(ret);
    }

    return ret;

}   // evaluateDesign

/*
================================================================================
Name:       getNObjectives

Function:   Returns the number of objectives of this optimization problem

Params:     (none)

Returns:    (int)   the number of objectives
================================================================================
*/
int
SystemDesign::getNObjectives() const
{
    return 1;

}   // getNObjectives

/*
================================================================================
Name:       getNConstraints

Function:   Returns the number of constraints of this optimization problem

Params:     (none)

Returns:    (int)   the number of constraints
================================================================================
*/
int
SystemDesign::getNConstraints() const
{
    return 3;

}   // getNConstraints
	

// PRIVATE MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       readDataFile

Function:   Reads a .csv file that contains the problem-specific parameters.
            See documentation for the correct file format

Params:     filename....the name of the .csv file with the problem data

Returns:    (void)
================================================================================
*/
void
SystemDesign::readDataFile(const string& filename)
{
    ifstream ifs(filename);
    string line;
    istringstream iss;
    
    // read the number of subsystems and the cost limit
    getline(ifs, line);
    replace(line.begin(), line.end(), ',', ' ');
    iss.str(line);

    iss >> nSubsystems_;
    iss >> costLimit_;

    // read the tech opt costs of each subsystem
    costs_.resize(nSubsystems_);

    for (size_t i = 0; i != nSubsystems_; ++i) {

        getline(ifs, line);
        replace(line.begin(), line.end(), ',', ' ');
        iss.str(line);
        double cost;

        while (iss >> cost) {
            costs_[i].push_back(cost);
        }

        iss.clear();
    }

    // read the tech opt utilities of each subsystem
    utilities_.resize(nSubsystems_);

    for (size_t i = 0; i != nSubsystems_; ++i) {

        getline(ifs, line);
        replace(line.begin(), line.end(), ',', ' ');
        iss.str(line);
        double utility;

        while (iss >> utility) {
            utilities_[i].push_back(utility);
        }

        iss.clear();
    }

    // read the necessitations
    getline(ifs, line);
    replace(line.begin(), line.end(), ',', ' ');
    iss.str(line);

    size_t nNecessitations;
    iss >> nNecessitations;

    for (size_t i = 0; i != nNecessitations; ++i) {

        getline(ifs, line);
        replace(line.begin(), line.end(), ',', ' ');

        iss.str(line);

        size_t driverSubSysId, driverTecOptId, targetSubSysId, targetTecOptId;

        iss >> driverSubSysId;
        iss >> driverTecOptId;
        iss >> targetSubSysId;
        iss >> targetTecOptId;

        necessitations_.push_back(
            Dependency(
                TechOpt(driverSubSysId, driverTecOptId), 
                TechOpt(targetSubSysId, targetTecOptId)
                )
            );
    }

    // read the obviations
    getline(ifs, line);
        replace(line.begin(), line.end(), ',', ' ');

    //iss.clear();
    iss.str(line);
    size_t nObviations;
    iss >> nObviations;

    for (size_t i = 0; i != nObviations; ++i) {

        getline(ifs, line);
        replace(line.begin(), line.end(), ',', ' ');

        iss.str(line);

        size_t driverSubSysId, driverTecOptId, targetSubSysId, targetTecOptId;

        iss >> driverSubSysId;
        iss >> driverTecOptId;
        iss >> targetSubSysId;
        iss >> targetTecOptId;

        obviations_.push_back(
            Dependency(
                TechOpt(driverSubSysId, driverTecOptId), 
                TechOpt(targetSubSysId, targetTecOptId)
                )
            );
    }

}   // readDataFile


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators