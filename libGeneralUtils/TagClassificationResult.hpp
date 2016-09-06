#ifndef FEATUREINFOCLASSIFICATIONRESULT_H
#define FEATUREINFOCLASSIFICATIONRESULT_H
/*!
    \file   TagClassificationResult.h
    \brief  Provides a class for tagging a FeatureVector with a classification result
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    1.2.2010
*/

#include "Tag.hpp"
#include "ClassLabel.hpp"
#include "ClassLabelHistogram.hpp"

/*!
 * \class TagClassificationResult
 * \brief This class holds information about classification outcome
 */
class TagClassificationResult : public Tag
{
  public:
    TagClassificationResult() : label(NULL_CLASS_LABEL), histogram() {};
    TagClassificationResult(const TagClassificationResult &source) : label(source.label), histogram(source.histogram) {};
    TagClassificationResult(ClassLabel label_, ClassLabelHistogram histogram_) : label(label_), histogram(histogram_) {};
    virtual ~TagClassificationResult() {};

    ClassLabel label;
    ClassLabelHistogram histogram; //!< containes all classification-votes

    virtual Tag* clone() {return new TagClassificationResult(*this);};
    virtual TagType getType() {return 3;};
    static TagType Type() {return 3;};
    static std::string TypeStr() {return "TagClassificationResult";};
};

#endif
