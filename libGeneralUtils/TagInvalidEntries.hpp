/*
 * TagInvalidEntries.h
 *
 *  Created on: Aug 2, 2010
 *      Author: sauerland
 */

#ifndef TAGINVALIDENTRIES_H_
#define TAGINVALIDENTRIES_H_

#include "Tag.hpp"

#include <vector>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>


/*!
 * \class TagInvalidEntries
 * \brief This tag holds which feature dimensions are invalid
 */
class TagInvalidEntries : public Tag
{
  private:
    std::vector<unsigned int> entries;

  public:
    TagInvalidEntries() : entries(){};
    TagInvalidEntries(const TagInvalidEntries &source) : entries(source.entries) {};
    TagInvalidEntries(std::vector<unsigned int> invalid_entries) : entries(invalid_entries){};
    virtual ~TagInvalidEntries() {};

    bool isValid(unsigned int index) const{if(find(entries.begin(),entries.end(),index) == entries.end()) return true;return false;};
    void add(unsigned int index){entries.push_back(index);};

    virtual Tag* clone() {return new TagInvalidEntries(*this);};
    virtual TagType getType() {return 4;};
    static TagType Type() {return 4;};
    static std::string TypeStr() {return "TagInvalidEntries";};

  private:
    //! for storing and loading this class with boost:serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      (void)version; // avoids warning "unused parameter"
      using boost::serialization::make_nvp;
      ar & BOOST_SERIALIZATION_NVP(entries);
    }
};

#endif /* TAGINVALIDENTRIES_H_ */
