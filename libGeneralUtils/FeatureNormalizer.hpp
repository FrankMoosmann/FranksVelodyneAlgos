#ifndef FEATURENORMALIZER_H_
#define FEATURENORMALIZER_H_

/*!
    \file   FeatureNormalizer.h
    \brief  Provides a base class for normalizing features
    \author Frank Moosmann (<frank.moosmann@kit.edu>)
    \date   29.11.2011

    Copyright: Frank Moosmann
*/

#include "FeatureVector.hpp"
#include <vector>
#include <list>

/*!
 * \class FeatureNormalizer
 * \brief Base class for normalizing features
 *
 * This class provides base functionality to normalize Features. Normalization is in principle independent for each dimension.
 * Derive from this base class if normalization is not bijective (cannot be undone), otherwise derive from FeatureBijectiveNormalizer (below)
 */
class FeatureNormalizer
{
public:
  virtual ~FeatureNormalizer() {};

  //! normalize data at given dimension
  virtual FeatureData normalize(FeatureData data, size_t dim) const = 0;

  //! normalize one feature vector. might throw depending on implemented normalize method
  void normalize(FeatureVector &fVector) const {
    for (size_t i = 0; i < fVector.size(); i++) {
      fVector[i] = normalize(fVector[i], i);
    }
  }

  //! normalize all feature vectors. might throw depending on implemented normalize method
  template <class InputIterator>
  void normalize(InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd) const {
    while (FeatureVectorBegin != FeatureVectorEnd) {
      normalize(*FeatureVectorBegin);
      ++FeatureVectorBegin;
    }
  }
};

/*!
 * \class FeatureBijectiveNormalizer
 * \brief Base class for normalizing features
 *
 * Derive from this class of normalization is bijective, i.e. if it can be undone.
 */
class FeatureBijectiveNormalizer : public FeatureNormalizer
{
public:
  virtual ~FeatureBijectiveNormalizer() {};

  //! undo normalization of (normalized) data at given dimension
  virtual FeatureData denormalize(FeatureData data, size_t dim) const = 0;
};



#endif /* FEATURENORMALIZER_H_ */
