#ifndef FEATUREINFO_H
#define FEATUREINFO_H
/*!
    \file   Tag
    \brief  Provides a base class for tagging FeatureVectors with different type of information
    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    17.05.2008
*/

typedef unsigned int TagType;

/*!
 * \class Tag
 * \brief Derive from this class to tag FeatureVectors with additional information
 *
 * Each derived class must have its own unique ID as TagType.
 * So far assigned:
 * 0 = Base Tag
 * 1 = TagImagePos
 * 2 = TagTrueClassLabel
 * 3 = TagClassificationResult
 * 4 = TagInvalidEntries
 * 5 = TagID
 */
class Tag
{
  public:
    //Tag() = 0;
    virtual ~Tag() {};
    
    virtual Tag* clone() {return NULL;};
    virtual TagType getType() {return 0;};
    static TagType Type() {return 0;};
    static std::string TypeStr() {return "";};
};

#endif
