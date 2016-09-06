#ifndef CLASSIFIER_H
#define CLASSIFIER_H

/*!
    \file   Classifier.h
    \brief  Provides a base class for classification
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    15.01.2006
*/

#include <FeatureVector.hpp>
#include <ClassLabel.hpp>
#include <ClassLabelHistogram.hpp>
#include <TagClassificationResult.hpp>
#include <TagTrueClassLabel.hpp>
#include <list>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/function.hpp>

/*
  The following code is a tiny example on how Classifiers can be constructed from an array of already existing data:
  
  float data[] = {0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.9,1.0, 0.4,0.23,0.45,0.56,0.98, 0.01,0.002,0.513,0.34,0.123};
  int32_t label[] = {1,1,2,2};
  int data_dims[] = {5, 4};
  int label_dims[] = {1, 4};

  int data_dim = data_dims[0];
  int data_num = data_dims[1];
  int label_dim = label_dims[0];
  int label_num = label_dims[1];

  assert(label_dim==1);
  assert(data_num==label_num);

  FeatureSet features;
  for (int i=0; i<data_num; i++) {
    FeatureVector feature(data+i*data_dim,data_dim);
    ClassLabel c = label[i];
    feature.setTag(TagTrueClassLabel(c));
    //cout << feature.getTag<TagTrueClassLabel>().label << endl;
    features.push_back(feature);
  }

  ExtraTreeClassifier::ExtraTreeParams extraTreeParams;
  Classifier *myclassifier = NULL;
  cout << "training..." << endl;
  myclassifier = new ExtraTreeClassifier(features.begin(),features.end(),extraTreeParams);
  cout << "training finished" << endl;
  delete myclassifier;

  return 0;
*/

/*!
 * \class Classifier
 * \brief This is an abstract base class for implementing classifiers
 */
class Classifier {

public:

  virtual ~Classifier() {};
  virtual void writeToFile(std::string filename) const = 0;

  /*!
    \brief Method for classifying new data
    \param featureVec the FeatureVector to classify
    \return the estimated class. More details will be stored within the featureVec as TagClassificationResult
  */
  ClassLabel classify (FeatureVector &featureVec);

  /*!
    \brief Method for classifying new data into only two classes with a variable threshold
    \param featureVec the FeatureVector to classify
    \param threshold if confidence value of normal classification is greater that this threshold, estimated class will be positive
    \param negativeClass class to be treated as "positive". Influences return value
    \param negativeClass class to be treated as "negative". Influences return value
    \return true if the estimated class is the positive class, false otherwise. More details will be stored within the featureVec as TagClassificationResult
    \throws runtime_error if the classifier was learned for more than two classes
  */
  bool binaryClassify (FeatureVector &featureVec, float threshold, ClassLabel positiveClass, ClassLabel negativeClass);
  
  typedef boost::function<void (void)> Callable;
  /*!
    \brief Method for classifying several features into only one class estimate. The calculated error will be available if the Features are tagged with "TrueClassLabel"
    \param first Iterator to first FeatureVector to be classified
    \param last Iterator beyond last FeatureVector to be classified
    \param error  optional pointer to a variable where classification error will be stored
    \param totalIntegratedHist  optional pointer to the histogram where all votes (to retrieve confidence values) will be stored
    \param perFeatIntegratedHist optional pointer to the histogram where FeatureVector based classification result votes (to retrieve confidence values) will be stored
    \param status a callable (stuct/class having operator() defined) that is called after each classification. Can be used to update a progress bar.
  */
  template<class InputIterator>
  void classify (InputIterator first, InputIterator last, double *error = NULL, ClassLabelHistogram *totalIntegratedHist = NULL, ClassLabelHistogram *perFeatIntegratedHist = NULL, Callable call = Callable()) {
    double err = 0.0;
    double cnt = 0.0;
    for (InputIterator k=first; k != last; k++) {
      ClassLabel featureClassLabel;
      FeatureVector &fVec = *k;
      try {
        featureClassLabel = fVec.getTag<TagTrueClassLabel>().label;
      } catch (std::exception &e) {
        featureClassLabel = NULL_CLASS_LABEL;
      }
      ClassLabel classifiedLabel = classify(fVec);
      assert(classifiedLabel == fVec.getTag<TagClassificationResult>().label);
      if (totalIntegratedHist != 0) {
        ClassLabelHistogram h = fVec.getTag<TagClassificationResult>().histogram;
        totalIntegratedHist->integrate(h);
      }
      if (perFeatIntegratedHist != 0)
        perFeatIntegratedHist->increaseCountOfLabel(classifiedLabel);
      if (classifiedLabel != featureClassLabel)
        err += 1.0;
      cnt += 1.0;
//      std::cout << "." << std::flush;
      if (!call.empty())
        call();
    }
    if (error) *error = err/cnt;
  };

  template <class OutputIterator>
  void getTrainedClassLabels(OutputIterator i) const { //!< get all class labels that occurred during training
    BOOST_FOREACH(ClassLabel c, classLabelList) {
      ++i = c;
    }
  }

protected:
  /*!
    \brief Method for classifying new data
    \param featureVec the FeatureVector to classify
    \param hist histogram to store the "votings". this will be used to assign a class label
  */
  virtual void classify (FeatureVector &featureVec, ClassLabelHistogram &hist) = 0;

  //! list of all classLabels that occurred in the training data (should be set up in the constructor)
  std::list<ClassLabel> classLabelList;

};
#endif //CLASSIFIER_H
