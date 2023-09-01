//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// INCLUDES
// #############################################################################
#include "Weldbeam.hpp"

#include <cmath>


// USING DECLARATIONS
// #############################################################################
using std::map;
using std::pow;
using std::sqrt;
using std::vector;


// #############################################################################
//  BEGIN NAMESPACE
// #############################################################################
namespace evaluators {


// STRUCTORS
// #############################################################################
/*
================================================================================
Name:       Weldbeam

Function:   Instantiates a Weldbeam evaluator object.

Params:     trackingOn......default = false. sets whether to track rate of 
                            convergence.
================================================================================
*/
Weldbeam::Weldbeam(bool trackingOn)
    :   Evaluator(trackingOn)
{
}   // Weldbeam


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
/*
================================================================================
Name:       evaluateDesign

Function:   Evaluates the objective function and constraints for a given design

Params:     x....the design to evaluate.

Returns:    (vector<double>)    the objective and constraint values
================================================================================
*/
vector<double>
Weldbeam::evaluateDesign(const std::vector<double>& x)
{
    // Welded beam problem from adapted from Deb and Goyal, 1997. 
    // "Optimization Engineering Using a Combined Genetic Search".
    // Stress, buckling, and deflection equations from "Engineering Optimization
    // Methods and Applications" Reklaitis 1983. (Example 1.2)


    // CONSTANT PARAMETERS
    // ===========================================

    // material properties and costs - listed in the following order: 
    //  0-steel, 1-cast iron, 2-brass, 3-aluminmum
    // material strengths (psi)
    static const double strength[] = {30e3, 8e3, 5e3, 8e3};
    // Young's modulus (psi)
    static const double yMod[] = {30e6, 14e6, 10e6, 16e6};
    // shear modulus(psi)
    static const double sMod[] = {12e6, 6e6, 4e6, 6e6};
    // weld material cost ($/in^3)
    static const double wCost[] = {0.1047, 0.0489, 0.5235,  0.5584};
    // bar material cost ($/in^3)
    static const double bCost[] = {0.0481, 0.0224, 0.2405, 0.2566};

    // other constant parameters
    static const double F = 6000;       //force on end of beam (lb.)
    static const double L = 14;         //extended length of beam (in.)
    static const double dmax = 0.25;    //maximum allowed beam deflection (in.)


    // DESIGN VARIABLE-DEPENDENT PARAMETERS
    // ===========================================

    // material properties and costs
    int material(static_cast<int>(x[1]));
    double S = strength[material];      
    double E = yMod[material];      
    double G = sMod[material];      
    double c1 = wCost[material];        
    double c2 = bCost[material];

    // geometric parameters
    double h = x[2];    // weld bead thickness (in.)
    double t = x[3];    // bar thickness (in.)
    double b = x[4];    // bar width (in.)
    double l = x[5];    // welded length (in.)


    // OBJECTIVE AND CONSTRAINT CALCULATIONS
    // ===========================================

    // objective function: cost ($)
    double f = (1 + c1) * pow(h, 2.0) * (l + x[0] * t) + c2 * t * b * (L + l);

    // constraint 0: bending stress in beam
    double sig = 6.0 * F * L / (b * pow(t, 2.0));
    double g0 = S - sig;

    // constraint 1: beam buckling load
    double tmp = 1.0 / 3.0;
    double a = 1.0 / 3.0 * G * t * pow(b, 3.0);
    double I = 1.0 / 12.0 * t * pow(b, 3.0);
    double Pc = 4.013 * sqrt(E * I * a) / pow(L, 2.0) * (1.0 - t / (2.0 * L) * 
        sqrt(E * I / a));
    double g1 = Pc - F;

    // constraint 2: beam deflection
    double delta = 4.0 * F * pow(L, 3.0) / (E * pow(t, 3.0) * b);
    double g2 = dmax - delta;

    // constraint 3: weld stress
    double A = (1.0 - x[0]) * 1.414 * h * l + x[0] * 1.414 * h * (t + l);
    double M = F * (L + l / 2.0);
    double R = sqrt(pow(l, 2.0) / 4.0 + pow( ( (t + h) / 2.0), 2.0));
    double J = (1.0 - x[0]) * 1.414 * h * l * (pow( (h + t), 2.0) / 4.0 + 
        pow(l, 2.0) / 12.0) + x[0] * 1.414 * h * pow( (h + t + l), 3.0) / 12.0;

    double t1 = F / A;  double t2 = M * R / J;
    double tau = sqrt(pow(t1, 2.0) + 2.0 * t1 * t2 * l / (2.0 * R) + 
        pow(t2, 2.0));

    double g3 = 0.577 * S - tau;

    
    // return vector of objective and constraint values [f, g0, g1, g2, ..., gx]
    vector<double> ret;
    ret.push_back(f);
    ret.push_back(g0);
    ret.push_back(g1);
    ret.push_back(g2);
    ret.push_back(g3);

    // track rate of convergence
    if (trackingOn_)
        trackRateOfConvergence(ret);

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
Weldbeam::getNObjectives() const
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
Weldbeam::getNConstraints() const
{
    return 4;

}   // getNConstraints


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace evaluators