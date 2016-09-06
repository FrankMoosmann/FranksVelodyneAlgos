#ifndef CLASSLABELHISTOGRAM_H
#define CLASSLABELHISTOGRAM_H
/*!
    \file   ClassLabelHistogram.h
    \brief  Provides a class for voting in classification problems
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    16.11.2005

    Copyright: Frank Moosmann
*/

#include "ClassLabel.hpp"
#include "FeatureVector.hpp"
#include <deprecated.hpp>
#include <map>
#include <list>

/*!
 * \class ClassLabelHistogram
 * \brief This class can be used in classification problems to do voting
 *
 * It is used by classifiers to return more a detailed classification outcome
 * as solely the class label. Further, several ClassLabelHistograms can be
 * combined.
 */
class ClassLabelHistogram {

public:
  
  ClassLabelHistogram(); //!< constructs an empty histogram
  ~ClassLabelHistogram();
  
  void increaseCountOfLabel (ClassLabel label, float weight = 1.0); //!< the central function: increases votes for a specific label
  void integrate(ClassLabelHistogram &hist); //!< integrate another histograms votes into this one
  void reset(); //!< reset all to zero

  unsigned int getTotalCount(); //!< returns the number of votes performed in total (i.e. sum of voted for all classes)
template <class OutputIterator>
  void getListOfClassLabels (OutputIterator inserter); //!< returns all ClassLabels that were already voted for
  DEPRECATED(void getListOfClassLabels (std::list<ClassLabel> *classLabelList)); //!< \deprecated returns all ClassLabels that were already voted for

  int getCountOfLabel (ClassLabel label); //!< returns the number of counts for a specific class
  float getWeightOfLabel (ClassLabel label); //!< returns the accumulated weights for a specific class
  float getAvgWeightOfLabel (ClassLabel label); //!< returns weight/count for a specific class
  void classify(ClassLabel &label, float &confidence); //!< returns best class label according to weights and positive confidence value i.e. getWeightDiffOfTopTwoLabels()
  float getWeightDiffOfTopTwoLabels(); //!< always returns a positive weight. this can be interpreted as confidence value for class "getRealLabelWithMaxWeight()"
  DEPRECATED(ClassLabel getLabelWithMaxOccurance ()); //!< should not be used!!! may return NULL_CLASS_LABEL
  ClassLabel getRealLabelWithMaxOccurance (); //!< get second best label, if best class was NULL_CLASS_LABEL
  ClassLabel getRealLabelWithMaxWeight (); //!< get second best label, if best class was NULL_CLASS_LABEL

  friend std::ostream &operator<< (std::ostream &ostr, const ClassLabelHistogram &clh) {
    for (std::map<ClassLabel, CountWeightPair>::const_iterator i = clh.histogram.begin(); i != clh.histogram.end(); ++i) {
      ostr << i->first << ":" << i->second.first << ";" << i->second.second << "  ";
    }
    return ostr;
  };

private:
  unsigned int dataCount;
  typedef std::pair<int,float> CountWeightPair;
  typedef std::pair<ClassLabel, CountWeightPair> LabelCountWeightPair;
  std::map<ClassLabel, CountWeightPair> histogram; //!< stores for each class its count and accumulated weight
};


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////  template implementation  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template <class OutputIterator>
void ClassLabelHistogram::getListOfClassLabels(OutputIterator inserter)
{
  for (std::map<ClassLabel, CountWeightPair>::iterator p=histogram.begin(); (p!=histogram.end()); ++p)  {
    ++inserter = p->first;
  }
}


#endif //CLASSLABELHISTOGRAM_H
