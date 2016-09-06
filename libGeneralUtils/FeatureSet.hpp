#ifndef FEATURESET_H_
#define FEATURESET_H_
/*!
    \file   FeatureSet.h
    \brief  Provides a class for storing (classification) features
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    26.1.2010

    Copyright: Frank Moosmann
*/

#include <list>
#include <set>
#include <map>
#include <boost/foreach.hpp>
#include <boost/function.hpp>

#include "FeatureVector.hpp"
#include "Tag.hpp"

/*!
 * \class FeatureSet
 *
 * \brief Container class for several features
 *
 * Derive from this class to implement database-specific loading/access.
 * It can be used to serve as input to the Classifiers library.
 */
class FeatureSet {
public:
  FeatureSet();
  virtual ~FeatureSet();

  //////////////////////////////////////////////
  ////         access functions             ////
  //////////////////////////////////////////////
  typedef std::list<FeatureVector>::iterator iterator;
  typedef std::list<FeatureVector>::const_iterator const_iterator;
  typedef std::insert_iterator< std::list<FeatureVector> > fsInserter;
  iterator begin() {return featVecList.begin();};
  iterator end() {return featVecList.end();};
  const_iterator begin() const {return featVecList.begin();};
  const_iterator end() const {return featVecList.end();};

  //////////////////////////////////////////////
  ////         modifier functions           ////
  //////////////////////////////////////////////
  unsigned int size(void) const {return featVecList.size();};
  void clear();
  void push_back(FeatureVector f); //!< add one feature vector to FeatureSet
  void push_back(const FeatureSet &other) {push_back(other.begin(), other.end());}; //!< add all features of the other feature set into this one
  template <class FeatureVectortIterator>
  void push_back(FeatureVectortIterator begin, FeatureVectortIterator end); //!< add several feature vectors to FeatureSet
  fsInserter getInserter() {return fsInserter(featVecList,featVecList.end());};

  FeatureSet clone() const; //!< deep-copies all feature-vectors
  template <class DimIndexIterator>
  FeatureSet clone(DimIndexIterator begin, DimIndexIterator) const; //!< deep-copies all feature-vectors, but only the specified dimensions (must be given in ascending order, all if begin==end)
  void splice(FeatureSet &other); //!< efficiently moves all elements from the other set into this set
  //merge(other)
  void splitUniform(double factor, FeatureSet &second); //!< splits a fraction of its features into the second set
  template <class FeatureSetIterator>
  void splitUniformCopy(FeatureSetIterator targetBegin, FeatureSetIterator targetEnd) const; //!< splits features into chunks of equal size. creates copies of features
  template <class FeatureSetIterator, class PercentageIterator>
  void splitUniformCopy(FeatureSetIterator targetBegin, FeatureSetIterator targetEnd, PercentageIterator percentageBegin, PercentageIterator percentageEnd) const; //!< splits features into several chunks of specified sizes. creates copies of features

  // TODO: implement
  //    static int         readFeatureVectorsFromFile (std::string fileName, std::string directory, std::list<FeatureVector*> &featureList, int nb=0, std::map<unsigned int,bool> *featureDimensionsIgnored = 0); //!< Reads feature vectors from a file and returns the number of features that were read
  void saveToFile(std::string fileName);
  //    static void        writeIntoStream(std::list<FeatureVector*> &featureList, std::ostream &data, std::ostream &info, int maxcount = 0);
  //    static int         readFromStream(std::istream &data, std::istream &info, std::list<FeatureVector*> &featureList, int nbToRead = 0, std::map<unsigned int,bool> *featureDimensionsIgnored = 0);


  //////////////////////////////////////////////
  ////         statistic functions          ////
  //////////////////////////////////////////////
  template <class CustomTag>
  bool isTagged() const; //!< Returns true if all FeatureVectors contain this kind of Tag

protected:
  std::list<FeatureVector> featVecList; // currently, FeatureVectors have reference counting implemented, so memory is managed automatically
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////  template implementation  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template <class CustomTag>
bool FeatureSet::isTagged() const {
  const_iterator b = begin();
  const_iterator e = end();
  while (b != e) {
    if (!b->hasTag<CustomTag>())
      return false;
    ++b;
  }
  return true;
}

template <class FeatureVectortIterator>
void FeatureSet::push_back(FeatureVectortIterator begin, FeatureVectortIterator end) {
  while (begin != end) {
    push_back(*begin);
    ++begin;
  }
}

template <class DimIndexIterator>
FeatureSet FeatureSet::clone(DimIndexIterator diBegin, DimIndexIterator diEnd) const
{
  FeatureSet other;
  for (FeatureSet::const_iterator i = begin(); i != end(); ++i) {
      other.push_back(i->clone(diBegin, diEnd));
  }
  return other;
}

template <class FeatureSetIterator>
void FeatureSet::splitUniformCopy(FeatureSetIterator outputBegin, FeatureSetIterator outputEnd) const {
  unsigned int count = 0;
  FeatureSetIterator i = outputBegin;
  while (i!=outputEnd) {
    ++count;
    ++i;
  }
  std::vector<double> percentages(count,1.0/(double)count);
  splitUniformCopy(outputBegin, outputEnd, percentages.begin(), percentages.end());
}

template <class FeatureSetIterator, class PercentageIterator>
void FeatureSet::splitUniformCopy(FeatureSetIterator output, FeatureSetIterator outputEnd, PercentageIterator percentage, PercentageIterator percentageEnd) const {
  std::list<FeatureVector>::const_iterator fvi = featVecList.begin();
  double currentIndex = 0;
  double listsize = (double)size();
  double stopAt = 0;
  while ((output != outputEnd) && (percentage != percentageEnd)) {
    stopAt += (*percentage) * listsize;
    while ((currentIndex < stopAt) && (fvi != featVecList.end())) {
      output->push_back(*fvi);
      currentIndex += 1.0;
      ++fvi;
    }
    ++output;
    ++percentage;
  }
}

#endif /* FEATURESET_H_ */
