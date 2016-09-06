/*!
    \file   ClassLabelHistogram.cpp
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    16.11.2005
*/

#include <iostream>
#include <cfloat>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include "ClassLabelHistogram.hpp"

using namespace std;

 
//**************************************************************************
//                         Constructors / Destructors
//**************************************************************************
ClassLabelHistogram::ClassLabelHistogram()
{
  dataCount = 0;
}

ClassLabelHistogram::~ClassLabelHistogram()
{
}


//**************************************************************************
//                         public Methods
//**************************************************************************
void ClassLabelHistogram::increaseCountOfLabel(ClassLabel label, float weight)
{
  dataCount++;
  BOOST_AUTO(it,histogram.find(label));
  if (it != histogram.end()) {
    it->second.first += 1;
    it->second.second += weight;
  } else {
    histogram.insert(LabelCountWeightPair(label,CountWeightPair(1, weight)));
  }
}

void ClassLabelHistogram::integrate(ClassLabelHistogram &other)
{
  dataCount += other.dataCount;
  for (std::map<ClassLabel, CountWeightPair>::iterator i=other.histogram.begin(); i!=other.histogram.end(); ++i) {
    ClassLabel label   = i->first;
    unsigned int count = i->second.first;
    float weight       = i->second.second;
    BOOST_AUTO(it,histogram.find(label));
    if (it != histogram.end()) {
      it->second.first += count;
      it->second.first += weight;
    } else {
      histogram.insert(LabelCountWeightPair(label,CountWeightPair(count, weight)));
    }
  }
}
  
void ClassLabelHistogram::reset ()
{
  histogram.clear();
  dataCount = 0;
}


unsigned int ClassLabelHistogram::getTotalCount()
{
  return dataCount;
}

void ClassLabelHistogram::getListOfClassLabels(list<ClassLabel>* classLabelList)
{
  getListOfClassLabels(inserter(*classLabelList, classLabelList->end()));
}

int ClassLabelHistogram::getCountOfLabel (ClassLabel label)
{
  if (histogram.find(label) == histogram.end())
    return 0;
  return histogram[label].first;
}

float ClassLabelHistogram::getWeightOfLabel (ClassLabel label)
{
  if (histogram.find(label) == histogram.end())
    return 0.0;
  return histogram[label].second;
}

float ClassLabelHistogram::getAvgWeightOfLabel (ClassLabel label)
{
  BOOST_AUTO(it,histogram.find(label));
  if (it == histogram.end())
    return 0.0;
  //cout << endl << "weightsum " << it->second.second << " div count " << it->second.first;
  float result = it->second.second/((float)it->second.first);
  //cout << " return " << result;
  return result;
}

void ClassLabelHistogram::classify(ClassLabel &label, float &confidence)
{ //1:2;1.31016  2:2;-1.69018  3:2;0.38002
  float maxWeight1 = -FLT_MAX;
  float maxWeight2 = -FLT_MAX;
  label = NULL_CLASS_LABEL;
  for(map<ClassLabel, CountWeightPair>::iterator p=histogram.begin(); (p!=histogram.end()); ++p) {
    // iterators allow access to:
    // p->first = label
    // p->second.first = accumulated votes
    // p->second.second = accumulated confidences
    if (p->first != NULL_CLASS_LABEL) { // only process the classes different from NULL
      if (p->second.second > maxWeight1) {
        label = p->first;
        maxWeight2 = maxWeight1;
        maxWeight1 = p->second.second;
      } else if (p->second.second > maxWeight2) {
        maxWeight2 = p->second.second;
      }
    }
  }
  confidence = maxWeight1 - maxWeight2;
  if (maxWeight2 == -FLT_MAX) // no vote for any other class
    confidence = fabs(maxWeight1);
}

float ClassLabelHistogram::getWeightDiffOfTopTwoLabels()
{
  ClassLabel label; float confidence;
  classify(label, confidence);
  return confidence;
} 


ClassLabel ClassLabelHistogram::getLabelWithMaxOccurance ()
{
  int max = 0;
  ClassLabel maxLabel = NULL_CLASS_LABEL;
  for(map<ClassLabel, CountWeightPair>::iterator p=histogram.begin(); (p!=histogram.end()); ++p) {
    if ((*p).second.first > max) {
      max = (*p).second.first;
      maxLabel = (*p).first;
    }
  }
  return maxLabel;
}


ClassLabel ClassLabelHistogram::getRealLabelWithMaxOccurance ()
{
  int max = 0;
  ClassLabel maxLabel = NULL_CLASS_LABEL;
  for(map<ClassLabel, CountWeightPair>::iterator p=histogram.begin(); (p!=histogram.end()); ++p) {
    if (p->first != NULL_CLASS_LABEL) { // only process the classes different from NULL
      if (p->second.first > max) {
        max = (*p).second.first;
        maxLabel = (*p).first;
      }
    }
  }
  return maxLabel;
}


ClassLabel ClassLabelHistogram::getRealLabelWithMaxWeight ()
{
  ClassLabel label; float confidence;
  classify(label, confidence);
  return label;
}


//**************************************************************************
//                         private Methods
//**************************************************************************

// --- none ---
