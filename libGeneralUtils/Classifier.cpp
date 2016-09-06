#include "Classifier.hpp"

ClassLabel Classifier::classify (FeatureVector &featureVec) {
  ClassLabelHistogram hist;
  classify(featureVec, hist);
  ClassLabel label = hist.getRealLabelWithMaxWeight();
  featureVec.setTag(TagClassificationResult(label, hist));
  return label;
};

bool Classifier::binaryClassify (FeatureVector &featureVec, float threshold, ClassLabel positiveClass, ClassLabel negativeClass) {
  if (classLabelList.size() != 2) throw std::runtime_error("Classifier::binaryClassify: class count is not 2");
  ClassLabelHistogram hist;
  classify(featureVec, hist);
  ClassLabel classifiedLabel; float confidence;
  hist.classify(classifiedLabel, confidence); // returns an always-positive confidence value
  if (classifiedLabel == negativeClass) confidence = -confidence;
  classifiedLabel = (confidence <= threshold) ? negativeClass : positiveClass;
  featureVec.setTag(TagClassificationResult(classifiedLabel, hist));
  return (classifiedLabel == positiveClass);
};
