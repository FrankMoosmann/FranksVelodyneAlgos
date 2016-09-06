/*!
    \file   ClassLabelHistogram.cpp
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    07.11.2005
*/

#include "FeatureVector.hpp"


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>

using namespace std;


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////     Class  FeatureVectorDataBase      /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

FeatureVectorDataBase::FeatureVectorDataBase()
{
  refCount = 0;
}

FeatureVectorDataBase::~FeatureVectorDataBase()
{
  BOOST_FOREACH(TagMap::value_type t, tags) {
    delete t.second;
  }
}

FeatureVectorDataBase* FeatureVectorDataBase::clone() const
{
  FeatureVectorDataBase* target = cloneData();
  cloneTags(target);
  return target;
}

void FeatureVectorDataBase::cloneTags(FeatureVectorDataBase* target) const
{
  BOOST_FOREACH(const TagMap::value_type &t, tags) {
    TagMap::value_type tc(t);
    tc.second = t.second->clone(); // make a copy of the tag
    target->tags.insert(tc);
  }
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////        Class  FeatureVector           /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

/********************************
 *         Constructors
 *******************************/

FeatureVector::FeatureVector() {
  fvdata = new FeatureVectorDenseSTLData();
  fvdata->refCount = 1;
}

FeatureVector::FeatureVector(const FeatureVector &other) {
  fvdata = other.fvdata;
  fvdata->refCount++;
}

FeatureVector::FeatureVector(FeatureVectorDataBase *data) {
  fvdata = data;
  fvdata->refCount = 1;
}

FeatureVector::~FeatureVector() {
  deinit();
};

FeatureVector& FeatureVector::operator= (const FeatureVector &other) {
  // check for self-assignment  
  if (this == &other)  
    return *this;
  deinit();
  fvdata = other.fvdata;
  fvdata->refCount++;
  return *this;
}
FeatureVector& FeatureVector::operator+= (const FeatureVector &other) {
  unsigned int size = fvdata->size();
  if (size != other.fvdata->size())
    throw length_error("FeatureVector::operator+= Both vectors must be of the same length");
  for (unsigned int i=0; i<size; ++i){
    fvdata->setElement(i, fvdata->getElement(i) + other.fvdata->getElement(i));
  }
  return *this;
}
FeatureVector& FeatureVector::operator/ (const float value) {
  unsigned int size = fvdata->size();
  for (unsigned int i=0; i<size; ++i){
    fvdata->setElement(i, fvdata->getElement(i) / value);
  }
  return *this;
}
bool FeatureVector::operator== (const FeatureVector &other) const {
  if (fvdata == other.fvdata) { // both point to the same data object
    return true;
  } else { // check every data element
    unsigned int size = fvdata->size();
    if (size != other.fvdata->size())
      return false;
    for (unsigned int i=0; i<size; ++i){
      if (fvdata->getElement(i) != other.fvdata->getElement(i))
        return false;
    }
  }
  return true;
}

//std::ostream &operator<<(std::ostream &ostr, const FeatureVector &f) {
//  BOOST_FOREACH(FeatureData d, *(f.data))
//    ostr << d << " ";
//  return ostr;
//}

void FeatureVector::deinit() {
   fvdata->refCount--;
   if (fvdata->refCount == 0)
    delete fvdata;
}

/********************************
 *       public Members
 *******************************/

//FeatureVector::iterator FeatureVector::begin() {
//  return data->begin();
//}
//
//FeatureVector::iterator FeatureVector::end() {
//  return data->end();
//}

bool FeatureVector::hasTag(const TagType type) const {
  BOOST_AUTO(i,fvdata->tags.find(type));
  return (i!=fvdata->tags.end());
}

Tag* FeatureVector::getTag(const TagType type) const {
  BOOST_AUTO(i,fvdata->tags.find(type));
  if (i != fvdata->tags.end()) {
    return i->second;
  } else {
    return NULL;
  }
}

/********************************
 *       static Members
 *******************************/

float FeatureVector::cityblockDistance(const FeatureVector &f1, const FeatureVector &f2)
{
  if (f1.fvdata == f2.fvdata) return 0.0; // both point to the same data in memory
  if (f1.size() != f2.size()) throw out_of_range("FeatureVector::cityblockDistance: features have different size");
  float sum = 0.0;
  unsigned int size = f1.size();
  for (unsigned int i=0; i<size; ++i) {
    sum += fabs(f1(i) - f2(i));
  }
//  the following can be activated again, when iterators are available (but rewrite, so that it is faster for sparse data):
//  vector<FeatureData>::const_iterator i1=f1.data->begin();
//  vector<FeatureData>::const_iterator i2=f2.data->begin();
//  for (; i1 != f1.data->end(); ++i1,++i2) {
//    sum += fabs(*i1 - *i1);
//  }
  return sum;
}

float FeatureVector::euclideanDistance(const FeatureVector &f1, const FeatureVector &f2)
{
  if (f1.fvdata == f2.fvdata) return 0.0; // both point to the same data in memory
  if (f1.size() != f2.size()) throw out_of_range("FeatureVector::euclideanDistance: features have different size");
  float sum = 0.0;
  unsigned int size = f1.size();
  for (unsigned int i=0; i<size; ++i) {
    float diff = (f1(i) - f2(i));
    sum += diff*diff;
  }
  //  the following can be activated again, when iterators are available (but rewrite, so that it is faster for sparse data):
//  vector<FeatureData>::const_iterator i1=f1.data->begin();
//  vector<FeatureData>::const_iterator i2=f2.data->begin();
//  for (; i1 != f1.data->end(); ++i1,++i2) {
//    float diff = (*i1 - *i1);
//    sum += diff*diff;
//  }
  return sqrt(sum);
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////   Class  FeatureVectorDenseSTLData    /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

FeatureVectorDenseSTLData::FeatureVectorDenseSTLData()
  : FeatureVectorDataBase()
{
}

FeatureVectorDenseSTLData::FeatureVectorDenseSTLData(const FeatureVector &fv)
  : FeatureVectorDataBase()
{
  unsigned int s = fv.size();
  data.reserve(s);
  for (unsigned int idx=0; idx<s;++idx)
    data.push_back(fv(idx));
}

FeatureVectorDenseSTLData::FeatureVectorDenseSTLData(unsigned int size, FeatureData initValue)
  : FeatureVectorDataBase()
   ,data(size,initValue)
{
}

FeatureVectorDenseSTLData::~FeatureVectorDenseSTLData()
{
}

FeatureVectorDataBase* FeatureVectorDenseSTLData::cloneData() const
{
  FeatureVectorDenseSTLData *newd = new FeatureVectorDenseSTLData();
  newd->data = vector<FeatureData>(this->data); // deep-copy data
  return newd;
}

