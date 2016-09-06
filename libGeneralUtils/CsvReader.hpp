/*
 * CsvReader.h
 *
 *  Created on: Feb 26, 2010
 *      Author: moosmann
 */

#ifndef CSVREADER_H_
#define CSVREADER_H_

#include <vector>
#include <string>
#include <sstream>

class CsvReader {
public:
  CsvReader(std::string filename, std::string separators = " ;", std::string commentIndicators = "", bool extractEmtpyLines = false, bool extractEmptyFields = true);
  virtual ~CsvReader();

  unsigned int lineCount(void) const; //!< returns the number of lines extracted, if case extractEmtpyLines==true it should correspond to the number of lines in the text file
  unsigned int minFieldCount(void) const; //!< returns the number of fields that all lines minimally contain
  unsigned int fieldCount(unsigned int line) const; //!< the the number of values extracted from the specified line. line must be in 0..lineCount-1
  std::string get(unsigned int line, unsigned int field) const; //!< the value at the specified position as it appeared in the text file. line must be in 0..lineCount-1 and field in 0..fieldCount(line)
  template <typename return_t>
  return_t get(unsigned int line, unsigned int field) const { //!< the value at the specified position converted to return_t. line must be in 0..lineCount-1 and field in 0..fieldCount(line)
    std::string s = get(line,field);
    std::istringstream iss(s);
    return_t t;
    iss >> t;
    return t;
  };
private:
  std::vector< std::vector<std::string>* > fields;
};

#endif /* CSVREADER_H_ */
