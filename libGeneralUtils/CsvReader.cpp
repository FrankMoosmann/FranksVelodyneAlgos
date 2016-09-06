/*
 * CsvReader.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: moosmann
 */

#include "CsvReader.hpp"

#include <fstream>
#include <stdexcept>
#include <climits>
#include <limits>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include "FMUtils.hpp"

using namespace std;
using namespace boost::algorithm;

CsvReader::CsvReader(string filename, string separators, string commentIndicators, bool extractEmtpyLines, bool extractEmptyFields)
{
  //cout << "reading " << filename << endl;
  if (commentIndicators.length()>1)
    cerr << endl << "WARNING: CsvReader can currently only handle one commentIndicator" << flush;
  ifstream desc(filename.c_str());
  if (!desc.good())
    throw invalid_argument("CsvReader: specified file "+filename+" does not exist");
  string line;
  vector<string> *splitVec = new vector<string>;
  while (!desc.eof()) {
    getline(desc,line);
//    cout << "line:" << line << endl;
    trim(line);
    if (commentIndicators != "") {
      if (line[0] == commentIndicators[0])
      //if (starts_with(line, commentIndicators, is_any_of(commentIndicators)))
        continue;
    }
    if (line.length() > 0) {
      split( *splitVec, line, is_any_of(separators)); // split line into parts, alternatively use: !(is_digit() || is_punct())
      if (!extractEmptyFields) {
        splitVec->erase(remove(splitVec->begin(), splitVec->end(), ""), splitVec->end());
      }
    }
    if ((line.length() > 0) || (extractEmtpyLines)) {
      fields.push_back(splitVec);
//      cout << "->size " << fields.size() << endl;
      splitVec = new vector<string>;
    }
  }
//  cout << "Csv:finished" << endl;
}

CsvReader::~CsvReader() {
  BOOST_FOREACH(std::vector<std::string>* vi, fields)
    delete vi;
}

unsigned int CsvReader::lineCount(void) const
{
  return fields.size();
}
unsigned int CsvReader::minFieldCount(void) const
{
  if (fields.size() == 0)
    return 0;
  unsigned int mn = numeric_limits<unsigned int>::max();
  for (BOOST_AUTO(l,fields.begin()); l!= fields.end(); ++l) {
    mn = min(mn,(unsigned int)(**l).size());
  }
//  BOOST_AUTO(l,fields.begin());
//  unsigned int mn = l->size();
//  ++l;
//  while (l!=fields.end()) {
//    mn = min(mn,l->size());
//    ++l;
//  }
  return mn;
}

unsigned int CsvReader::fieldCount(unsigned int line) const
{
  if (line >= fields.size())
    throw out_of_range("CsvReader::fieldCount: line is out of range");
  return fields[line]->size();
}

std::string CsvReader::get(unsigned int line, unsigned int field) const
{
  if (line >= fields.size())
    throw out_of_range("CsvReader::get: line is out of range");
  if (field >= fields[line]->size())
    throw out_of_range("CsvReader::get: field is out of range");
  return (*fields[line])[field];
}

