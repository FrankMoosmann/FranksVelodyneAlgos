#ifndef FEATUREINFOTRUECLASSLABEL_H
#define FEATUREINFOTRUECLASSLABEL_H
/*!
    \file   TagTrueClassLabel.h
    \brief  Provides a class for tagging a FeatureVector with a class label
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    1.2.2010
*/

#include <set>
#include <map>
#include <sstream>
#include <boost/foreach.hpp>
#include "Tag.hpp"
#include "FeatureVector.hpp"
#include "ClassLabel.hpp"

/*!
 * \class TagTrueClassLabel
 * \brief This class holds information about classification outcome
 */
class TagTrueClassLabel : public Tag
{
  public:
    TagTrueClassLabel() : label(NULL_CLASS_LABEL) {};
    TagTrueClassLabel(const TagTrueClassLabel &source) : label(source.label) {};
    TagTrueClassLabel(ClassLabel label_) : label(label_) {};
    virtual ~TagTrueClassLabel() {};

    ClassLabel label;

    virtual Tag* clone() {return new TagTrueClassLabel(*this);};
    virtual TagType getType() {return 2;};
    static TagType Type() {return 2;};
    static std::string TypeStr() {return "TagTrueClassLabel";};

    template <class InputIterator>
    static void        attachToAll(InputIterator featureVectorBegin, InputIterator featureVectorEnd, ClassLabel label); //!< attaches the specified label to all FeatureVectors as TrueLabelTag
    static ClassLabel  getLabel(const FeatureVector &fv) {return fv.getTag<TagTrueClassLabel>().label;}; //!< return the true class label of a feature vector \throws runtime_error if no label was attached to the vector
    template <class InputIterator, class OutputIterator>
    static void        getLabels(InputIterator featureVectorBegin, InputIterator featureVectorEnd, OutputIterator classLabelInserter); //!< adds all occurring class labels into the inserter \throws run_time_error if features have no TrueClassLabel tag associated
    template <class InputIterator, class PairOutputIterator>
    static void        getLabelCounts(InputIterator featureVectorBegin, InputIterator featureVectorEnd, PairOutputIterator labelCountPairInserter); //!< get the number of occurrences of all class labels \throws run_time_error if features have no TrueClassLabel tag associated
    template <class InputIterator>
    static std::string getLabelDistribution(InputIterator featureVectorBegin, InputIterator featureVectorEnd); //!< get list of all occurring class labels as string \throws run_time_error if features have no TrueClassLabel tag associated
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////  template implementation  ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template <class InputIterator>
void TagTrueClassLabel::attachToAll (InputIterator featureVectorIterator, InputIterator FeatureVectorEnd, ClassLabel label)
{
  while (featureVectorIterator != FeatureVectorEnd) {
    FeatureVector &cFeatVec = *featureVectorIterator;
    cFeatVec.setTag(TagTrueClassLabel(label));
    ++featureVectorIterator;
  }
}

template <class InputIterator, class OutputIterator>
void TagTrueClassLabel::getLabels(InputIterator featureVectorIterator, InputIterator featureVectorEnd, OutputIterator classLabelInserter)
{
  std::set<ClassLabel> lList;
  while (featureVectorIterator != featureVectorEnd) {
    const Tag* t = featureVectorIterator->getTag(TagTrueClassLabel::Type());
    const TagTrueClassLabel* ttc = dynamic_cast<const TagTrueClassLabel*>(t);
    lList.insert(ttc->label);
    ++featureVectorIterator;
  }
  BOOST_FOREACH(ClassLabel c, lList) {
    ++classLabelInserter = c;
  }
}

template <class InputIterator, class PairOutputIterator>
void TagTrueClassLabel::getLabelCounts(InputIterator featureVectorIterator, InputIterator featureVectorEnd, PairOutputIterator labelCountPairInserter)
{
  std::map<ClassLabel, unsigned int> lMap;
  std::map<ClassLabel, unsigned int>::iterator lMapIt;
  typedef std::pair<ClassLabel, unsigned int> pairT;
  while (featureVectorIterator != featureVectorEnd) {
    const Tag* t = featureVectorIterator->getTag(TagTrueClassLabel::Type());
    if (t) {
      const TagTrueClassLabel* ttc = dynamic_cast<const TagTrueClassLabel*>(t);
      if ((lMapIt = lMap.find(ttc->label)) == lMap.end())
        lMap.insert(pairT(ttc->label,1));
      else
        lMap[ttc->label]++;
    }
    ++featureVectorIterator;
  }
  BOOST_FOREACH(pairT cp, lMap) {
    ++labelCountPairInserter = cp;
  }
}

template <class InputIterator>
std::string TagTrueClassLabel::getLabelDistribution(InputIterator featureVectorIterator, InputIterator featureVectorEnd)
{
  std::ostringstream ostr;
  std::map<ClassLabel, unsigned int> lMap;
  typedef std::pair<ClassLabel, unsigned int> pairT;
  TagTrueClassLabel::getLabelCounts(featureVectorIterator, featureVectorEnd, inserter(lMap,lMap.end()));
  BOOST_FOREACH(pairT p, lMap) ostr << " " << p.first << ":" << p.second;
  return ostr.str();
}


#endif
