/*!
    \file   EnumCreation.h
    \brief  Provides a preprocessor-macro for easy creation of enums including conversion to string and boost program_option parsing
    \author Frank Moosmann (<frank.moosmann@kit.edu>)
    \date   15.11.2011

    The macro functions here can be used to define enums.
    String-conversions and a parsing function for boost::program_options
    are automatically defined at the same time

    Below is a small example code on how to use this file

    Copyright: Frank Moosmann
*/

/*****************************************
            USAGE-EXAMPLE
******************************************
#include "EnumCreation.hpp"
#define MY_ENUM_SEQ (ENUM1)(ENUM2)(ENUM3)
CREATE_ENUM(MyEnumType, MY_ENUM_SEQ)            // defines the type "MyEnumType" and functions "to_string(MyEnumType)" and "to_MyEnumType(string)"
CREATE_ENUM_VALIDATOR(MyEnumType, MY_ENUM_SEQ)  // enables simple usage within boost::program_options

MyEnumType enumVar;
boost::program_options::options_description progOpt;
progOpt.add_options()
  ("help", "produce help message")
  ("enum", po::value(&enumVar)->default_value(ENUM1), "enum test")
;

std::string enumStr = to_string(enumVar);
MyEnumType backConvert = to_MyEnumType(enumStr);
cout << "enum=" << to_string(enumVar) << " / " << backConvert << endl;
switch (enumVar) {
  case ENUM1: cout << "e1" << endl; break;
  case ENUM2: cout << "e2" << endl; break;
  case ENUM3: cout << "e3" << endl; break;
}
******************************************/

#ifndef ENUMCREATION_H_
#define ENUMCREATION_H_

#include <boost/program_options.hpp>
#include <boost/preprocessor.hpp>

/****************************************
         "private" HELPER-MACROS
         don't use them directly
*****************************************/

// the following helper-macro generates the to-string conversion for the given element
// of the form       case CM_TREE: return "TREE";
#define GEN_ENUM_TOSTR(unused_r, unused_data, elem) \
  case elem: return BOOST_PP_STRINGIZE(elem);
// the following helper-macro generates the from-string conversion for the given element
// of the form       if (str == "TREE") return CM_TREE;
#define GEN_ENUM_FROMSTR(unused_r, unused_data, elem) \
  if (str == BOOST_PP_STRINGIZE(elem)) return elem;


/****************************************
           "public" MACROS
    the only functions to be used
        from outside this file
*****************************************/

// this macro defines a string with comma-separated enum values
#define CREATE_ENUM_LISTSTR(enumtypename, enumseqence) \
    BOOST_PP_STRINGIZE( BOOST_PP_SEQ_ENUM(enumseqence) )
//    BOOST_PP_SEQ_FOR_EACH(GEN_ENUM_STR,~,enumseqence)

// this macro defines the enum and two conversion functions (to string and from string)
// the array must simply contain the enums as entries
#define CREATE_ENUM(enumtypename, enumseqence) \
  enum enumtypename { BOOST_PP_SEQ_ENUM(enumseqence) }; \
  inline std::string to_string(enumtypename e) { \
    switch (e) { \
      BOOST_PP_SEQ_FOR_EACH( GEN_ENUM_TOSTR, ~, enumseqence ) \
      default: return ""; \
    } \
  }; \
  inline enumtypename BOOST_PP_CAT(to_,enumtypename)(std::string str) { \
    BOOST_PP_SEQ_FOR_EACH( GEN_ENUM_FROMSTR, ~, enumseqence ) \
    throw std::invalid_argument("enumtypename string is invalid"); \
  }; \
  inline std::ostream &operator<<(std::ostream &ostr, const enumtypename &enumvar) { \
    ostr << to_string(enumvar); \
    return ostr; \
  };

// this macro defines a function which allows boost::program_options to parse input into an enum type
// prerequisite: the enum must have been definied by the above macro
#define CREATE_ENUM_VALIDATOR(enumtypename, array) \
  void validate(boost::any& v, const std::vector<std::string>& values, enumtypename* target_type, int) \
  { \
    (void)target_type; \
    namespace po = boost::program_options; \
    po::validators::check_first_occurrence(v); \
    const std::string& s = po::validators::get_single_string(values); \
    try { \
      v = boost::any(BOOST_PP_CAT(to_,enumtypename)(s)); \
      return; \
    } catch (std::invalid_argument &e) { \
      throw po::validation_error(po::validation_error::invalid_option_value, "invalid value"); \
    } \
  };



#endif /* ENUMCREATION_H_ */
