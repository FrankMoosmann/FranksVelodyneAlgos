#ifndef SVMCLASSIFIER_H
#define SVMCLASSIFIER_H

/*!
    \file   SvmClassifier.h
    \brief  Provides a support vector machine classifier
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    14.01.2006
*/

#include "Classifier.hpp"
#include "libsvm/svm.h"
#include <list>
#include <map>
#include <cmath>
#include <boost/program_options/options_description.hpp>

/*!
 * \class SvmClassifier
 * \brief Implements a SupportVectorMachine
 *
 * note: especially for the linear kernel, all the features should be
 * in the same range (usually -1 to 1). if not, training may never
 * terminate (you'll see a lot of '.' on the screen)
 */
class SvmClassifier : public virtual Classifier {

public:

  struct SvmParams : public svm_parameter { //!< svm_parameter supplemented by boost::program_options parsing
    //class WeightMap : public std::map<ClassLabel, float> {}; // needed for program_options to work
    //WeightMap classWeights; //!< specify if classes should get different weights (default 1)
    std::vector<float> classWeights; //!< specify if classes should get different weights (default 1)
    bool weightBalanceTrainingSet; //!< if true determines weights above automatically from class distribution
    boost::program_options::options_description parameterDescription; //!< object allowing to retrieve all parameters via boost program options
    SvmParams(); //!< initializes default values and creates parameterDescription object
  };

  /*!
    \brief Constructor to learn the SVM from the examples
    \param FeatureVectorBegin iterator to the first training-FeatureVector
    \param FeatureVectorEnd iterator past the last training-FeatureVector
    \param params structure containing the parameters
  */
  template<class InputIterator>
  SvmClassifier(InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const SvmParams &params);
  /*!
    \brief Method for creating a SVM from the given SVM-Problem
    This method creates a SVM model based on the given problem data-structure
    \param problem the given training files as "SVM-problem" (see libsvm/svm.h)
    \param classHistogram a histogram of the classLabels
    \param params structure containing the parameters
  */
  SvmClassifier(struct svm_problem *problem, unsigned int featureDimensionality, ClassLabelHistogram *classHistogram, svm_parameter *params);
  /*!
    \brief Method to read SvmClassifier from a file
    \param filename file to create classifier from
  */
  SvmClassifier(std::string filename);
  /*!
    \brief Destructor
  */
  virtual ~SvmClassifier();
  /*!
    \brief Method to save SvmClassifier to file
    \param filename path/file to save to
  */
  virtual void writeToFile(std::string filename) const;

  struct hyperplane { // hyperplane for deciding about a class {pos,neg} for any point x by y(x) = w^T * x + b
    std::vector<double> w;
    double b;
    ClassLabel pos;
    ClassLabel neg;
  };

  //! \brief Method for getting the hyperplane parameters of linear SVM
  std::vector<hyperplane> getHyperplanes() const;

protected:

  virtual void classify (FeatureVector &featureVec, ClassLabelHistogram &hist);

private:

  struct svm_model *svm;     // pointer to support vector machine
  struct svm_problem prob;   // contains pointer to training data etc. Not possible to free memory if svm_model produced by svm_train() still in use !!!
  unsigned int dimensionality;
  struct svm_node *testNode; // used for classification
  double *decValues; // used for classification

  void convertFeatureVectors (std::list<FeatureVector*>* featureVec);
  void convertFeatureVector (FeatureVector* fVec, struct svm_node *x);

};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////  template implementation  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template<class InputIterator>
SvmClassifier::SvmClassifier(InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const SvmParams &params_)
{
  SvmParams params(params_); // create a copy to get a non-const pointer
  std::list<FeatureVector*> featVecList; // save pointers to list for "convert" algorithm
  ClassLabelHistogram totalHistogram;
  unsigned int count = 0;
  dimensionality = 0;
  // setup initial statistics
  while (FeatureVectorBegin != FeatureVectorEnd) {
    FeatureVector &fv = *FeatureVectorBegin;
    dimensionality = std::max(dimensionality,fv.getSize());
    featVecList.push_back(&fv);
    ClassLabel cl = fv.getTag<TagTrueClassLabel>().label;
    totalHistogram.increaseCountOfLabel(cl);
    ++FeatureVectorBegin;
    ++count;
  }
  totalHistogram.getListOfClassLabels(inserter(classLabelList,classLabelList.begin())); // set up classLabelList
  unsigned int classCount = classLabelList.size();
  // set weighting in order to balance training set (if activated)
  if ((params.weightBalanceTrainingSet) || (params.classWeights.size() != 0)) {
    params.nr_weight = classCount;
    params.weight_label = (int*)malloc(params.nr_weight * sizeof(int));
    params.weight = (double*)malloc(params.nr_weight*sizeof(double));
    ClassLabel maxLabel = totalHistogram.getRealLabelWithMaxOccurance();
    float maxVal = totalHistogram.getCountOfLabel(maxLabel);
    std::list<ClassLabel>::iterator it = classLabelList.begin();
    for (int i=0; i<params.nr_weight; ++i, ++it) {
      params.weight_label[i] = *it;
      if (params.weightBalanceTrainingSet) {
        params.weight[i] = maxVal/(float)(totalHistogram.getCountOfLabel(*it));
      } else {
        params.weight[i] = (i >= (int)params.classWeights.size()) ? 1.0 : params.classWeights[i];
      }
      std::cout << std::endl << "SVM-weight[" << params.weight_label[i] << "] = " << params.weight[i] << std::flush;
    }
  }
  convertFeatureVectors(&featVecList);
  svm = svm_train(&prob, &params);
  testNode = new svm_node[dimensionality+1]; // allocate space for classification. +1 because sparse list ends with additional entry "-1"
  decValues = new double[classCount*(classCount-1)/2]; // allocate space for classification.
}


#endif //SVMCLASSIFIER_H
