/*
 * FeatureMaxRangeNormalizer.cpp
 *
 *  Created on: Aug 9, 2010
 *      Author: sauerland
 */

#include "FeatureMaxRangeNormalizer.hpp"

using namespace std;


FeatureMaxRangeNormalizer::FeatureMaxRangeNormalizer(const FeatureMaxRangeNormalizer& other)
  : rangeMin(other.rangeMin), rangeMax(other.rangeMax), cropAtLimits(other.cropAtLimits)
{
}

FeatureMaxRangeNormalizer::FeatureMaxRangeNormalizer(std::vector<FeatureData> rangeMin_, std::vector<FeatureData> rangeMax_, bool cropAtLimits_)
  : rangeMin(rangeMin_), rangeMax(rangeMax_), cropAtLimits(cropAtLimits_)
{
}

FeatureMaxRangeNormalizer& FeatureMaxRangeNormalizer::operator=(FeatureMaxRangeNormalizer const& other)
{
  rangeMin = other.rangeMin;
  rangeMax = other.rangeMax;
  cropAtLimits = other.cropAtLimits;
  return *this;
}

FeatureMaxRangeNormalizer::~FeatureMaxRangeNormalizer() {
}

FeatureData FeatureMaxRangeNormalizer::denormalize(FeatureData data, size_t dim) const
{
  if (rangeMin.size() == 0) return data;
  // cannot really undo if cropping is activated!
  return data * (rangeMax[dim] - rangeMin[dim]) + rangeMin[dim];
}

FeatureData FeatureMaxRangeNormalizer::normalize(FeatureData data, size_t dim) const
{
  if (rangeMin.size() == 0) return data;
  if (cropAtLimits) data = min(max(data,rangeMin[dim]),rangeMax[dim]);
  return (data - rangeMin[dim]) / (rangeMax[dim] - rangeMin[dim]);
}

void FeatureMaxRangeNormalizer::normalize(FeatureSet::iterator begin, FeatureSet::iterator end, vector<float> gainList) const {
  if (rangeMin.size() == 0) return;
  if (gainList.size() < begin->size())
    throw out_of_range(" FeatureMaxRangeNormalizer::normalize; gainList is smaller than feature vector");
  unsigned int nDim = begin->size();
  FeatureSet::iterator currentFvec = begin;
  while (currentFvec != end) {
    //scale feature
    currentFvec->setElement(0, currentFvec->getElement(0) * gainList[0]);
    //normalized features
    for (unsigned int i = 1; i < nDim; i++) {
      float m = (1.0 * gainList[i]) / (rangeMax[i] - rangeMin[i]);
      float c = -m * rangeMin[i];
      currentFvec->setElement(i, m * currentFvec->getElement(i) + c);
      if (cropAtLimits) currentFvec->setElement(i, min(max(currentFvec->getElement(i),0.0f),1.0f));
    }
    currentFvec++;
  }
}

std::ostream &operator<<(std::ostream &ostr, const FeatureMaxRangeNormalizer &f)
{
  ostr << f.rangeMin.size() << "\t";
  for (unsigned int i = 0; i < f.rangeMax.size(); i++) {
    ostr << f.rangeMin[i] << " " << f.rangeMax[i] << "\t";
  }
  ostr << (f.cropAtLimits ? "1" : "0") << "\t";
  return ostr;
}

std::istream &operator>>(std::istream &istr, FeatureMaxRangeNormalizer &f)
{
  size_t s;
  istr >> s; // size
  f.rangeMin.resize(s);
  f.rangeMax.resize(s);
  for (size_t i=0; i<s; ++i) {
    istr >> f.rangeMin[i];
    istr >> f.rangeMax[i];
  }
  istr >> f.cropAtLimits;
  return istr;
}
