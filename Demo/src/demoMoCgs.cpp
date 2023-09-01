//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


#include "demo.hpp"

#include "CGS\include\Domain.hpp"
#include "CGS\include\MoCgs.hpp"
#include "CGS\include\MultiResult.hpp"
#include "CGS\include\TrainingPoint.hpp"

#include "Evaluators\include\Knapsack2D.hpp"

#include <iostream>
#include <string>
#include <vector>

using bayesnet::Domain;
using bayesnet::TrainingPoint;
using cgs::MoCgs;
using cgs::MultiResult;

using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace demo {

void
demoMoCgs()
{
    // Begin demo 
    cout << "Begin MoCgs demo:" << endl;
    cout << "Pareto results will be saved to working directory..." << endl;

    // declare random seed
    unsigned int seed;

    // 2D KNAPSACK PROBLEM
    // =================================
    seed = 1;
    cout << "Seed: " << seed << endl;

    // construct vector of domains
    vector<Domain> ks2d_20_doms(20, Domain(2));

    // construct the evaluator object
    string filename20 = "..\\..\\Evaluators\\data\\knapsack20.txt";
    evaluators::Knapsack2D ks2d_20_eval(filename20);

    // construct the solver object
    MoCgs ks2d_20_solver(ks2d_20_eval, ks2d_20_doms, seed);

    // set solver parameters
    int ks2d_20_maxEvals(2000);
    int ks2d_20_nBest = 20;
    int ks2d_20_nInitRandom = 20;
    int ks2d_20_batchSize = 20;
    int ks2d_20_parentLimit = 2;

    // set converger(s)
    ks2d_20_solver.setMaxEvals(ks2d_20_maxEvals);

    // run the solver
    ks2d_20_solver.solve(
        ks2d_20_nBest,
        ks2d_20_nInitRandom,
        ks2d_20_batchSize,
        ks2d_20_parentLimit
        );

}   // demoMoCgs

}   // namespace demo