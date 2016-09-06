#ifndef FEATUREMAXRANGENORMALIZER_H_
#define FEATUREMAXRANGENORMALIZER_H_
/*
 * FeatureMaxRangeNormalizer.h
 *
 *  Created on: Aug 9, 2010
 *      Author: sauerland
 */

#include <climits>
#include <iostream>

#include "FeatureNormalizer.hpp"
#include <FeatureSet.hpp>
#include <FeatureVector.hpp>
#include <TagInvalidEntries.hpp>


/*!
 * \class FeatureMaxRangeNormalizer
 * \brief This class can be used to linearly normalize Features
 *
 * This class provides the functionality to normalize each feature dimension independently by: f' = f * m + c
 * Especially when learned from training data, the ranges [min, max] are mapped to [0, 1]
 *
 * Outlier handling can be activated during training
 */
class FeatureMaxRangeNormalizer : public FeatureBijectiveNormalizer
{
private:
  std::vector<FeatureData> rangeMin;
  std::vector<FeatureData> rangeMax;
  bool cropAtLimits; // if true data beyond rangeMin/Max is simply cropped to Min/Max

public:
  struct OutlierMethod {
    FeatureData lowQuantile; //!< ratio (percentage*100) of data that is ignored at the lower bounds
    FeatureData highQuantile; //!< ratio (percentage*100) of data that is ignored at the upper bounds
    size_t minNbFVals; //!< ignore quantiles only to the extent that this number of different values remain
    FeatureData maxFVal; //!< training values higher than this are ignored
    OutlierMethod() //!< default constructor, de-activates outlier handling
    : lowQuantile(0.0)
    , highQuantile(0.0)
    , minNbFVals(0)
    , maxFVal(std::numeric_limits<FeatureData>::max())
    {}
    OutlierMethod(FeatureData lowQuantile_, FeatureData highQuantile_, size_t minNbFVals_, FeatureData maxFVal_)
    : lowQuantile(lowQuantile_)
    , highQuantile(highQuantile_)
    , minNbFVals(minNbFVals_)
    , maxFVal(maxFVal_)
    {}
  };
  FeatureMaxRangeNormalizer() : cropAtLimits(false) {};
  template <class InputIterator>
  FeatureMaxRangeNormalizer(InputIterator featureVecBegin, InputIterator featureVecEnd, OutlierMethod method = OutlierMethod(), bool cropAtLimits = false); //!< determines the normalization ranges from sample data
  FeatureMaxRangeNormalizer(const FeatureMaxRangeNormalizer &other); //!< copy constructor
  FeatureMaxRangeNormalizer(std::vector<FeatureData> rangeMin, std::vector<FeatureData> rangeMax, bool cropAtLimits = false); //!< normalization ranges can be specified by hand
  FeatureMaxRangeNormalizer &operator=(FeatureMaxRangeNormalizer const& source); //!< assignment operator

  virtual ~FeatureMaxRangeNormalizer();

  virtual FeatureData denormalize(FeatureData data, size_t dim) const;
  virtual FeatureData normalize(FeatureData data, size_t dim) const;

  void normalize(FeatureSet::iterator begin, FeatureSet::iterator end, std::vector<float> gainList) const; //normalize data [0..1*gain], but first entry is normalized: new = old*gain

  friend std::ostream &operator<<(std::ostream &ostr, const FeatureMaxRangeNormalizer &f);
  friend std::istream &operator>>(std::istream &istr, FeatureMaxRangeNormalizer &f);
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////  template implementation  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
template <class InputIterator>
FeatureMaxRangeNormalizer::FeatureMaxRangeNormalizer(InputIterator featureVecBegin, InputIterator featureVecEnd, OutlierMethod omethod, bool cropAtLimits_)
  : cropAtLimits(cropAtLimits_)
{
  if (featureVecBegin == featureVecEnd)
    throw std::range_error("FeatureMaxRangeNormalizer: no input data given");

  unsigned int nDim = featureVecBegin->size(); // assume that all features have same size
  rangeMax.resize(nDim, -std::numeric_limits<FeatureData>::max());
  rangeMin.resize(nDim, std::numeric_limits<FeatureData>::max());
  std::vector<FeatureData> values; // list of values for a specific dimension
  std::vector<size_t> valCount; // used to calculate how many different values exist

  for (unsigned int d = 0; d < nDim; d++) {
    // setup list of values for this dimension
    values.clear();
    for (InputIterator featureVecIt=featureVecBegin; featureVecIt != featureVecEnd; ++featureVecIt) {
      FeatureVector &fv = *featureVecIt;
      if ( (fv.hasTag<TagInvalidEntries>())
           && (!(fv.getTag<TagInvalidEntries>().isValid(d))) )
        continue;
      if (fv[d] > omethod.maxFVal)
        continue;
      values.push_back(fv[d]);
    }
    // sort list so that values[0] holds smallest value and values[size-1] the highest
    std::sort(values.begin(), values.end());
    if ((values.size() == 0) || (values[0] == values[values.size()-1])) {
      rangeMin[d] = 0;
      rangeMax[d] = 1;
    } else {
      // create statistic how many values are used
      // since values are sorted, a sequential pass is sufficient and finally
      // valCount[i] holds the number of different values up to position i
      valCount.resize(values.size());
      if (omethod.minNbFVals > 0) {
        valCount[0] = 1;
        for (unsigned int f = 1; f < values.size(); f++) {
          valCount[f] = valCount[f-1] + (values[f]==values[f-1] ? 0 : 1);
        }
      }
      // set range based on specified quantiles
      size_t pLow = 0 + (size_t)(omethod.lowQuantile*(float)values.size());
      size_t pHig = values.size()-1 - (size_t)(omethod.highQuantile*(float)values.size());
      //std::cout << "0 " << pLow << " " << pHig << " " << values.size()-1;
      // check if range violates "minNbFVals" parameter and adjust if necessary
      if (omethod.minNbFVals > 0) {
        while ((valCount[pHig]-valCount[pLow]+1 < omethod.minNbFVals) && (pLow > 0) && (pHig < values.size()-1)) {
          // move pLow/pHig towards end
          if (pLow > 0) --pLow;
          if (pHig < values.size()-1) ++pHig;
        }
      }
      //std::cout << " --> " << pLow << " " << pHig << "   " << valCount[pLow] << " " << valCount[pHig] << " " << valCount[values.size()-1] << std::endl;
      rangeMin[d] = values[pLow];
      rangeMax[d] = values[pHig];
    }
  } // for each dimension
}

#endif /* FEATUREMAXRANGENORMALIZER_H_ */
