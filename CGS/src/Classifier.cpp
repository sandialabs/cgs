//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//


// FILE CONTENTS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the implementation of the Classifier class. 
////////////////////////////////////////////////////////////////////////////////


// INCLUDES
// #############################################################################
#include "Classifier.hpp"

#include "Domain.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iterator>
#include <numeric>


// USING DECLARATIONS, DIRECTIVES, AND NAMESPACE ALIASES
// #############################################################################
using std::accumulate;
using std::cout;
using std::distance;
using std::endl;
using std::ifstream;
using std::inner_product;
using std::string;
using std::vector;


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// STRUCTORS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a Classifier object.
///
/// \param  domains             The attribute variable domains.
/// \param  nClasses            The number of categories in this classifier.
/// \param  graphFilename       A file with the initial graph edge list.
/// \param  pseudocount         The initial DiscreteDistribution weights.
/// \param  storeTrainingPoints Tells the classifier whether to store all of 
///                             the training points. Required when the 
///                             containsPoint or learnNetwork functions are to 
///                             be used.
////////////////////////////////////////////////////////////////////////////////
Classifier::Classifier(
    const std::vector<Domain>& domains, 
    int nClasses,
    const std::string& graphFilename,
    double pseudocount, 
    bool storeTrainingPoints
    )
    :   domains_(domains),
        hasFixedPriors_(false),
        nClasses_(nClasses),
        nVariables_(static_cast<int>(domains.size())),
        pseudocount_(pseudocount),
        priorProbabilities_(DiscreteDistribution(nClasses, pseudocount)),
        storeTrainingPoints_(storeTrainingPoints),
        trainingPts_(nClasses),
        bayesianNetwork_(domains, trainingPts_, graphFilename)
{
    // prior probabilities P(c) initialized to uniform distribution with 
    // user-defined pseudocount in the initialization list
    
    // initialize class cond. probabilities P(x|c) to uniform distributions with
    // user-defined pseudocount
    this->initializeConditionals();

}   // Classifier


////////////////////////////////////////////////////////////////////////////////
/// \brief  Instantiates a Classifier object.
///
/// \param  domains             The attribute variable domains.
/// \param  nClasses            The number of categories in this classifier.
/// \param  graph               A DirectedGraph object from whicht to construct
///                             the initial BayesianNetwork graph.
/// \param  pseudocount         The initial DiscreteDistribution weights.
/// \param  storeTrainingPoints Tells the classifier whether to store all of 
///                             the training points. Required when the 
///                             containsPoint learnNetwork functions are to be 
///                             used.
////////////////////////////////////////////////////////////////////////////////
Classifier::Classifier(
    const std::vector<Domain>& domains, 
    int nClasses,
    const DirectedGraph& graph,
    double pseudocount, 
    bool storeTrainingPoints
    )
    :   domains_(domains),
        hasFixedPriors_(false),
        nClasses_(nClasses),
        nVariables_(static_cast<int>(domains.size())),
        pseudocount_(pseudocount),
        priorProbabilities_(DiscreteDistribution(nClasses, pseudocount)),
        storeTrainingPoints_(storeTrainingPoints),
        trainingPts_(nClasses),
        bayesianNetwork_(domains, trainingPts_, graph)
{
    // prior probabilities P(c) initialized to uniform distribution with 
    // user-defined pseudocount in the initialization list
    
    // initialize class cond. probabilities P(x|c) to uniform distributions with
    // user-defined pseudocount
    this->initializeConditionals();

}   // Classifier


// PUBLIC MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Gets the number of classes in the classifier object.
///
/// \return The number of classes (categories) modeled by this classifier.
////////////////////////////////////////////////////////////////////////////////
int
Classifier::getNClasses() const
{
    return nClasses_;

}   // getNClasses


////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the prior probabilities to fixed values. Once this function is
///         called, the priors remain fixed according to the input parameter
///         and are not updated when new training points are added.
///
/// \param  weights The weights of the prior probability discrete distribution.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::setFixedPriors(const std::vector<double>& weights)
{
    hasFixedPriors_ = true;

    priorProbabilities_ = DiscreteDistribution(weights);

}   // setFixedPriors


////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the posterior probabilities, P(c|x), of a test point.
///
/// \param  testPt  An unobserved / unevaluated attribute vector.
///
/// \return The posterior probabilities in P(c|x).
////////////////////////////////////////////////////////////////////////////////
vector<double>
Classifier::computePosteriors(const std::vector<int>& testPt) const
{
    // compute posteriors using to Bayes formula: P(c|x) = P(c) * P(x|c) / P(x)
    // also can be stated as Posterior = Prior * Likelihood / Evidence

    // compute likelihoods
    vector<double> likelihoods;
    
    for (int i = 0; i != nClasses_; ++i) {
        
        double likelihood = 1.0;

        for (int j = 0; j != nVariables_; ++j) {
            likelihood *= 
                conditionalProbabilities_.at(j).at(i).getProbability(testPt);
        }

        likelihoods.push_back(likelihood);
    }

    // compute evidence
    vector<double> priorPs = priorProbabilities_.getProbabilities();

    double evidence = inner_product(
        likelihoods.begin(), likelihoods.end(),
        priorPs.begin(), 0.0
        );

    // compute posteriors
    vector<double> Pcx;

    for (int i = 0; i != nClasses_; ++i) {
        double prior = priorProbabilities_.getProbabilities()[i];
        Pcx.push_back(likelihoods[i] * prior / evidence);
    }
            
    return Pcx;

}   // computePosteriors


////////////////////////////////////////////////////////////////////////////////
/// \brief  Predicts the class label of a test point.
///
/// \param  testPt  An unobserved / unevaluated attribute vector.
///
/// \return The class label of the given test point.
////////////////////////////////////////////////////////////////////////////////
int
Classifier::classifyPoint (const std::vector<int>& testPt) const
{
    // compute posterior probabilities
    vector<double> Pxc = this->computePosteriors(testPt);

    // classLabel is the index of the max element in the posterior
    
    return static_cast<int>(
        distance(Pxc.begin(), max_element(Pxc.begin(), Pxc.end()))
        );

}   // classifyPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds a training point to the classifier.
///
/// \param  tp    A TrainingPoint object (attribute vector / class label pair).
////////////////////////////////////////////////////////////////////////////////
void
Classifier::addPoint(const TrainingPoint& tp)
{
    // if storeTrainingPoints is true, need to add the point to the training set
    if (storeTrainingPoints_) {

        trainingPts_.addPoint(tp);
    }

    // adjust the conditionals for this point
    double delta(1.0);
    this->adjustConditionals(tp, delta);

    // if priors are not fixed, adjust them 
    if (!hasFixedPriors_) {

        this->adjustPriors(tp.getClassLabel(), delta);
    }
    
}   // addPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Removes a training point from the classifier.
///
/// \param  tp  A TrainingPoint object (attribute vector / class label pair).
////////////////////////////////////////////////////////////////////////////////
void
Classifier::removePoint(const TrainingPoint& tp)
{
    // if storeTrainingPoints is true, need to remove the point from the 
    // training set
    if (storeTrainingPoints_) {

        trainingPts_.removePoint(tp);
    }

    // adjust the conditional probability distributions
    double delta(-1.0);
    this->adjustConditionals(tp, delta);

    // if priors are not fixed, adjust them
    if (!hasFixedPriors_) {

        this->adjustPriors(tp.getClassLabel(), delta);
    }

}   // removePoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Reclassifies a point. Effectively, this function removes a training
///         point from the classifier, assigns it a new class label, and then 
///         re-adds it to the classifier. 
///
/// \param  x       The attribute vector of the point to reclassify.
/// \param  oldC    The old class of the point.
/// \param  newC    The new class of the point.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::movePoint(const std::vector<int>& x, int oldC, int newC)
{

    this->removePoint(TrainingPoint(x, oldC));
    this->addPoint(TrainingPoint(x, newC));

}   // movePoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Reads a set of training points from a .csv file.
///
/// \param  filename    The name of a .csv file from which to read the points.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::addPointsFromFile(const std::string& filename)
{
    TrainingPoint tp;
    ifstream ifs(filename);

    while (tp.readCsv(ifs)) {

        this->addPoint(tp);
    }

}   // addPointsFromFile


////////////////////////////////////////////////////////////////////////////////
/// \brief  Checks if a training point exists in the training dataset.
///
/// \param  tp      The TrainingPoint to test for existence.
///
/// \return         Whether the point is in the training set. Note that this 
///                 function always returs false if the storeTrainingPoints
///                 flag is set to false.
////////////////////////////////////////////////////////////////////////////////
bool
Classifier::containsPoint(const TrainingPoint& tp) const
{
    return trainingPts_.containsPoint(tp);

}   // containsPoint


////////////////////////////////////////////////////////////////////////////////
/// \brief  Attempts to learn the bayesian network structure from the training
///         data and updates all distributions accordingly.
///
/// \param  parentLimit     An upper limit on the number of parents any node 
///                         may have.
/// \param  searchOrder     The order in which the K2 algorithm should traverse
///                         the nodes in its search for the best network.
///                         must be a permutation of integers 0 to nNodes - 1.
/// \param  learningClass   An optional parameter that enables the caller to 
///                         specify that the network be learned only from data 
///                         that has the specified class label. Default behavior
///                         is that data from all classes should be used.
/// \param  updatingClass   An optional parameter that enables the
///                         caller to specify that only the specified class's 
///                         class-conditional distributions be retrained when
///                         a new Bayesian network structure is learned. default
///                         behavior that all classes be retrained according to
///                         the latest Bayesian network.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::learnNetworkAndUpdate(
        int parentLimit, 
        const std::vector<int>& searchOrder, 
        int learningClass,
        int updatingClass
    )
{
    // learn the network from the current training dataset
    bool foundNewNetwork = bayesianNetwork_.learnNetworkFromData(
        parentLimit, searchOrder, learningClass
        );

    // retrain the distributions if a new better network was found
    if (foundNewNetwork) {

        // create a vector of the classes whose distributions should be updated
        vector<int> updatingClasses;

        if (updatingClass == -1) {
            
            // default. use all classes.
            for (int i = 0; i != nClasses_; ++i) {
                updatingClasses.push_back(i);
            }

        } else if (updatingClass < 0 || updatingClass > nClasses_ - 1) {

            // specified class is out of range. update all classes and warn
            for (int i = 0; i != nClasses_; ++i) {
                updatingClasses.push_back(i);
            }

            cout << "Warning (Classifier.cpp): User-specified post-network "
                "learning class is out of range. All classes updated." << endl;

        } else {

            // use only the specified class
            updatingClasses.push_back(updatingClass);
        }


        // retrain the specified distributions using the training data
        for (vector<int>::const_iterator j = updatingClasses.begin(); 
                j != updatingClasses.end(); 
                ++j) {
            
            // reset the conditional distributions for this class
            this->resetClassDistributions(*j);

            // add each training point with this class label
            TrainingDataset::FeatureVectorCountMap cTrainingPts(
                trainingPts_.getClassTrainingPointCounts(*j)
                );
                    
            TrainingDataset::FeatureVectorCountMap::const_iterator i;
            for (i = cTrainingPts.begin(); i != cTrainingPts.end(); ++i) {
            
                double delta = static_cast<double>(i->second);    
                TrainingPoint tp(i->first, *j);

                // adjust the conditionals for this point
                this->adjustConditionals(tp, delta);

                // if priors are not fixed, adjust them 
                if (!hasFixedPriors_) {

                    this->adjustPriors(tp.getClassLabel(), delta);
                }
            }
        }
    }

}   // learnNetworkAndUpdate


////////////////////////////////////////////////////////////////////////////////
/// \brief  Prints a representation of the Bayesian network.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::printBayesianNetwork()
{
    bayesianNetwork_.printParentList();

}   // printBayesianNetwork


// PRIVATE MEMBER FUNCTIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Initializes the conditionalProbabilities_ member variable to 
///         uniform distributions with a user-defined pseudocount.
///         The \link conditionalProbabilities_ \endlink member variable being 
///         initialized is a vector of vectors of ProbabilityTables. The number 
///         of elements is equal to the number of attribute variables.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::initializeConditionals()
{
    // clear the current contents of the conditional distributions
    conditionalProbabilities_.clear();

    // initialize conditional probability tables with pseudocount
    for (int i = 0; i != nVariables_; ++i) {

        ProbabilityTable pt(domains_, bayesianNetwork_, i, pseudocount_);

        vector<ProbabilityTable> Px(nClasses_, pt);

        conditionalProbabilities_.push_back(Px);
    }

}   // initializeConditionals

////////////////////////////////////////////////////////////////////////////////
/// \brief  Adjusts the weight of one class variable in the prior probability
///         distribution P(c). Called when a new training point is added to or
///         removed from the classifier, unless \link hasFixedPriors_ \endlink 
///         is set to true.
///
/// \param  cLabel      The class label whose weight is being changed.
/// \param  delta       The amount to adjust the weight.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::adjustPriors(int cLabel, double delta)
{
    priorProbabilities_.adjustWeight(cLabel, delta);

}   // adjustPriors


////////////////////////////////////////////////////////////////////////////////
/// \brief  Adjusts the weights of the attribute variables in the conditional
///         probability distributions P(x|c). Called when a new point is added
///         to the classifier.
///
/// \param  tp      The attribute vector / class label point pair.
/// \param  delta   The amount to adjust the weight.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::adjustConditionals(const TrainingPoint& tp, double delta)
{
    vector<int> x = tp.getXValues();
    int cl = tp.getClassLabel();
    
    for (int i = 0; i != nVariables_; ++i) {

        conditionalProbabilities_[i][cl].adjustWeight(x, delta);
    }

}   // adjustConditionals


////////////////////////////////////////////////////////////////////////////////
/// \brief  Resets all distributions to their initial states (i.e., all
///         distribution weights set to the pseudocount).
////////////////////////////////////////////////////////////////////////////////
void
Classifier::resetAllDistributions()
{
    // reset priors if they have not been declared fixed
    if (!hasFixedPriors_) {

        priorProbabilities_ = DiscreteDistribution(nClasses_, pseudocount_);
    }

    // reset conditional probability distributions
    this->initializeConditionals();


}   // resetAllDistributions


////////////////////////////////////////////////////////////////////////////////
/// \brief  Resets the distributions that are conditioned on the specified class
///         to their initial states (i.e., all distribution weights set to the 
///         pseudocount).
///
/// \param  cLabel      The label of the class whose distributions are to be
///                     reset.
////////////////////////////////////////////////////////////////////////////////
void
Classifier::resetClassDistributions(int cLabel) 
{
    // for each design variable
    for (int i = 0; i != nVariables_; ++i) {
        
        // construct a fresh probability table
        ProbabilityTable pt(domains_, bayesianNetwork_, i, pseudocount_);

        // assign its values to probability table that we wish to reset
        conditionalProbabilities_[i][cLabel] = pt;
    }
}

// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet