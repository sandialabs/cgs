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
///  \brief Contains the definition of class Classifier.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_CLASSIFIER_HPP
#define CGS_INCLUDE_CLASSIFIER_HPP


// INCLUDES
// #############################################################################
#include "BayesianNetwork.hpp"
#include "DiscreteDistribution.hpp"
#include "ProbabilityTable.hpp"
#include "TrainingPoint.hpp"
#include "TrainingDataset.hpp"

#include <random>
#include <string>
#include <unordered_map>
#include <vector>


// FORWARD DECLARATIONS
// #############################################################################
namespace bayesnet {
    class Domain;
}


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace bayesnet {


// CLASS DEFINITION
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  A Bayesian network classifier.
///
///         A Classifier, in the context of the \link bayesnet \endlink 
///         namespace, is a probabilistic classification model that is based on 
///         Bayes theorem. 
///         The Classifier is trained with a set of attribute vector / class
///         label pairs (TrainingDataset). Once trained, it can provide an 
///         estimate of the probability that an attribute vector belongs
///         to a certain class. The class which has the highest posterior
///         probability is the one in which that attribute vector is most likely
///         to be a member of.
///
///         The structure of the Bayesian network may be learned using the
///         appropriate member function, specified upon construction, or assumed
///         to be naive Bayes (default).
////////////////////////////////////////////////////////////////////////////////
class Classifier
{
    // PUBLIC DECLARATIONS
    // =================================
public:
    // STRUCTORS
    Classifier(
        const std::vector<Domain>& domains, 
        int nClasses,
        const std::string& graphFilename = "",
        double pseudocount = 1.0, 
        bool storeTrainingPoints = false
        );


    Classifier(
        const std::vector<Domain>& domains, 
        int nClasses,
        const DirectedGraph& graph,
        double pseudocount = 1.0, 
        bool storeTrainingPoints = false
        );


    // MEMBER FUNCTIONS
    int
    getNClasses() const;

    void
    setFixedPriors(const std::vector<double>& weights);

    std::vector<double>
    computePosteriors(const std::vector<int>& testPoint) const;
    
    int
    classifyPoint(const std::vector<int>& testPoint) const;

    void
    addPoint(const TrainingPoint& pointPair);

    void
    removePoint(const TrainingPoint& pointPair);

    void
    movePoint(const std::vector<int>& xValues, int oldClass, int newClass);

    void
    addPointsFromFile(const std::string& filename);

    bool
    containsPoint(const TrainingPoint& ptp) const;

    void
    learnNetworkAndUpdate(
        int parentLimit, 
        const std::vector<int>& searchOrder, 
        int trainingClass = -1,
        int updatingClass = -1
        );

    void
    printBayesianNetwork();

    template <class T>
    std::vector<int>    
    sampleXVal(int classLabel, T& generator);


    // PRIVATE DECLARATIONS
    // =================================
private:
    // MEMBER FUNCTIONS
    void
    initializeConditionals();

    void
    adjustPriors(int classLabel, double delta);

    void
    adjustConditionals(const TrainingPoint& ptp, double delta);

    void
    resetAllDistributions();

    void
    resetClassDistributions(int cLabel);

 
    // DATA MEMBERS

    /// \brief  The Bayesian network, which models variable dependencies.
    BayesianNetwork bayesianNetwork_;

    /// \brief  The conditional probabilities of the classifier. 
    ///
    ///         The outer vector has one element per class. The inner vectors
    ///         have one ProbabilityTable per attribute variable. Each table
    ///         holds the likelihoods of each variable given the class and its
    ///         parents in the Bayesian network P(x|c, pa).
    std::vector< std::vector<ProbabilityTable> > conditionalProbabilities_;

    /// \brief  Whether this Classifier has fixed prior probabilities P(c).
    ///         
    ///         If set to true, the prior probabilities will not be altered
    ///         when training points are added, moved, or deleted from the
    ///         TrainingDataset.
    bool hasFixedPriors_;

    /// \brief  The prior probabilitieis of the Classifier, P(c).
    DiscreteDistribution priorProbabilities_;

    /// \brief  Whether all training points are stored in memory. 
    ///
    ///         Points must be stored for checking if the classifier contains a 
    ///         particular point or the BN is to be learned from the dataset.
    bool storeTrainingPoints_;

    /// \brief  The set of attribute vector / class label pairs.
    TrainingDataset trainingPts_;

    /// \brief  The domains of the attribute variables.
    const std::vector<Domain>& domains_;
    
    /// \brief  The number classes modeled.
    const int nClasses_;

    /// \brief  The number of attribute variables modeled.
    const int nVariables_;

    /// \brief  An initial count for initializing/resetting all distributions.
    const double pseudocount_;


};  // class Classifier


// TEMPLATE DEFINITIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Generates a new point by sampling the conditional probability
///         distributions P(x|c).
///
/// \param  classLabel      The class from which to sample a new point.
/// \param  generator       A pseudo-random number engine from the c++ <random> 
///                         header (e.g. default_random_engine, mt19937, etc.).
///
/// \return                 The newly generated attribute vector.
////////////////////////////////////////////////////////////////////////////////
template <class T>
std::vector<int>
Classifier::sampleXVal(int classLabel, T& generator)
{
    // in general, samples must happen in ancestral order
    std::vector<int> ret(nVariables_);
    std::vector<int> sampleOrder(bayesianNetwork_.getTopologicalOrder());

    typedef std::vector<int>::const_iterator iter;

    for (iter i = sampleOrder.begin(); i != sampleOrder.end(); ++i) {

        ret[*i] = conditionalProbabilities_.
            at(*i).at(classLabel).sampleDistribution(ret, generator);
    }

    return ret;

}   // sampleXVal


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace bayesnet


#endif  // CGS_INCLUDE_CLASSIFIER_HPP