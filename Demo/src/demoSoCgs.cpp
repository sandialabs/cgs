//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


#include "demo.hpp"

#include "CGS\include\Domain.hpp"
#include "CGS\include\SoCgs.hpp"
#include "CGS\include\utilities.hpp"

#include "Evaluators\include\Deceptive3.hpp"
#include "Evaluators\include\Knapsack.hpp"
#include "Evaluators\include\SystemDesign.hpp"
#include "Evaluators\include\WarehouseLocation.hpp"
#include "Evaluators\include\Weldbeam.hpp"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>


using std::cout;
using std::endl;
using std::mt19937;
using std::string;
using std::vector;

using bayesnet::Domain;
using cgs::SoCgs;
using utilities::algorithms::generateRandomSeed;


namespace demo {

void
demoSoCgs()
{
    // begin demo
    cout << "Begin SoCgs demo:" << endl;

    // declare random seed
    unsigned int seed;
    
    // 20 ITEM KNAPSACK PROBLEM
    // =================================
    
    // initialize random seed
    seed = 1;
    std::cout << "Seed: " << seed << std::endl;

    // construct vector of domains
    vector<Domain> ks20_doms(20, Domain(2));

    // construct the solver object
    string filename20("..\\..\\Evaluators\\data\\knapsack20.txt");
    evaluators::Knapsack ks20_eval(filename20);
    SoCgs ks20_solver(ks20_eval, ks20_doms, seed);

    // set solver parameters
    double ks20_targetF = -795;
    int ks20_maxEvals = 5000;
    int ks20_nBest = 20;
    int ks20_nInitRandom = 20;
    int ks20_batchSize = 20;
    int ks20_parentLimit = 0;

    // set convergers
    ks20_solver.setTargetObjective(ks20_targetF);
    ks20_solver.setMaxEvals(ks20_maxEvals);
    
    // run the solver
    ks20_solver.solve(
        ks20_nBest, 
        ks20_nInitRandom, 
        ks20_batchSize,
        ks20_parentLimit
        );
    

    // 50 ITEM KNAPSACK PROBLEM
    // =================================
    /*
    // initialize random seed
    //seed = time(NULL);
    seed = 1;
    std::cout << "Seed: " << seed << std::endl;

    // construct vector of domains
    vector<Domain> ks50_doms(50, Domain(2));

    // construct the solver object
    string filename50 = "..\\..\\Evaluators\\data\\knapsack50.txt";
    evaluators::Knapsack ks50_eval(filename50);
    SoCgs ks50_solver(ks50_eval, ks50_doms, seed);

    // run the solver
    double ks50_targetF = -1900;
    int ks50_maxEvals = 10000;
    int ks50_nBest = 50;
    int ks50_nInitRandom = 50;
    int ks50_batchSize = 50;
    int ks50_parentLimit = 0;

    // set the convergers
    ks50_solver.setTargetObjective(ks50_targetF);
    ks50_solver.setMaxEvals(ks50_maxEvals);
    
    // run the solver
    ks50_solver.solve(
        ks50_nBest, 
        ks50_nInitRandom,
        ks50_batchSize,
        ks50_parentLimit
        );
    */

    // SYSTEM DESIGN PROBLEM
    // =================================
    /*
    // initialize random seed
    seed = time(NULL);
    //seed = 1;
    std::cout << "Seed: " << seed << std::endl;

    // construct vector of domains
    vector<Domain> sd_doms;
    vector<int> sd_domainSizes({4,6,5,5,2,2,8,6,5,7,5,5,2,4,2});

    for (size_t i = 0; i != sd_domainSizes.size(); ++i) {

    sd_doms.push_back(Domain(sd_domainSizes[i]));
    }


    // construct the solver object
    string sdFilename("..\\..\\Evaluators\\data\\systemDesignInputFile.csv");
    evaluators::SystemDesign sd_eval(sdFilename);
    SoCgs sd_solver(sd_eval, sd_doms, seed);

    // run the solver
    double sd_targetF = -1022;
    int sd_maxEvals = 5000;
    int sd_nBest = 20;
    int sd_nInitRandom = 20;
    int sd_batchSize = 20;
    int sd_parentLimit = 2;

    // set convergers
    sd_solver.setTargetObjective(sd_targetF);
    sd_solver.setMaxEvals(sd_maxEvals);

    // run the solver
    sd_solver.solve(
    sd_nBest,
    sd_nInitRandom,
    sd_batchSize,
    sd_parentLimit
    );
    */

    
    // WAREHOUSE LOCATION PROBLEM
    // =================================
    /*
    // initialize random seed
    seed = time(NULL);
    //seed = 1;
    std::cout << "Seed: " << seed << std::endl;

    // construct vector of domains
    vector<Domain> wh_doms(10, Domain(5));

    // construct the solver object
    evaluators::WarehouseLocation wh_eval;
    SoCgs wh_solver(wh_eval, wh_doms, seed);

    // run the solver
    double wh_targetF = 383;
    int wh_maxEvals = 5000;
    int wh_nBest = 20;
    int wh_nInitRandom = 20;
    int wh_batchSize = 20;
    int wh_parentLimit = 0;

    // set convergers
    wh_solver.setTargetObjective(wh_targetF);
    wh_solver.setMaxEvals(wh_maxEvals);

    // run the solver
    wh_solver.solve(
        wh_nBest,
        wh_nInitRandom,
        wh_batchSize,
        wh_parentLimit
        );
    */

    // WELDED BEAM PROBLEM
    // =================================
    /*
    // initialize random seed
    seed = time(NULL);
    //seed = 1;
    std::cout << "Seed: " << seed << std::endl;

    // construct vector of domains
    vector<Domain> wb_doms;
    wb_doms.push_back(Domain(2));
    wb_doms.push_back(Domain(4));
    wb_doms.push_back(Domain(8, 0.0625, 0.5));
    wb_doms.push_back(Domain(21, 7.5, 10.0));
    wb_doms.push_back(Domain(16, 0.0625, 1));
    wb_doms.push_back(Domain(24, 0.125, 3.0));

    // construct the solver object
    evaluators::Weldbeam wb_eval;
    SoCgs wb_solver(wb_eval, wb_doms, seed);

    // solver parameters
    double wb_targetF = 1.9509;
    //double wb_targetF = 1.98;
    int wb_maxEvals = 5000;
    int wb_nBest = 20;
    int wb_nInitRandom = 20;
    int wb_batchSize = 20;
    int wb_parentLimit = 0;

    // set the convergers
    wb_solver.setTargetObjective(wb_targetF);
    wb_solver.setMaxEvals(wb_maxEvals);

    // run the solver
    wb_solver.solve(
        wb_nBest, 
        wb_nInitRandom,
        wb_batchSize,
        wb_parentLimit
        ); */ 
    
}   // demoSoCgs

}   // namespace demo