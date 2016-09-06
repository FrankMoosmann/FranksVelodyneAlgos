/***************************************************************************
                  FMUtils.h  -  some small helper methods
                             -------------------
    begin                : 04.11.2005
    copyright            : (C) 2005 by Frank Moosmann
    email                : frank.moosmann@inrialpes.fr
 ***************************************************************************/

#ifndef FMUTILS_H
#define FMUTILS_H

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/foreach.hpp>

/***************************************************************************
   String utility functions
 ***************************************************************************/
// to convert a string into another type
template<class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}
template<class T>
T from_string(const std::string& s, std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  T t;
  iss >> f >> t;
  return t;
}
// to convert another type into a string
template<typename T>
std::string to_string(const T& x) {
  std::ostringstream oss;
  oss << x;
  return oss.str();
}
// decompose a string into substrings
void explode(std::string s, std::string e, std::vector<std::string>* ret);
// trim all whitespaces at the end
void trim(std::string &s);

/***************************************************************************
   Functions for the list template
 ***************************************************************************/
// adds all list elements of origin into the target list (similar to splice, but does not remove originals)
template<class T>
void appendList(std::list<T> *origin, std::list<T> *target) {
  typename std::list<T>::iterator p;
  for (p=origin->begin(); p!=origin->end(); p++) {
    target->push_back(*p);
  }
}

// reads a text file into a container accessed via an output iterator
template<class OutputIterator>
void readTextFile(std::string filename, OutputIterator target) {
  std::string line;
  std::ifstream file(filename.c_str());
  if (!file.good())
    return;
  while (!file.eof())  {
    getline(file,line);
    ++target = line;
  }
}

// needed to cout vector of confidence values of the classifier
template <class vector_content_type>
std::ostream &operator<<(std::ostream &ostr, std::vector<vector_content_type> const &v)
{
  BOOST_FOREACH(vector_content_type e, v) {ostr << "\t" << e;};
  return ostr;
}


/***************************************************************************
   Function for easier output of debug information
   set DEBUG_LVL to control amount of output
 ***************************************************************************/
void debug_out(const int debug_level, const char* message);
void debug_out(const int debug_level, const char* message, double value1, double value2=0.0);
void debug_out(const int debug_level, const char* message, int value1, int value2=0);

/***************************************************************************
   Functions to get a random numbers
 ***************************************************************************/
void randomInit();
unsigned int randomInRange(const unsigned int lowest, const unsigned int highest);
double randomInRange(const double lowest, const double highest);
double randomGauss(double mean, double variance);

/***************************************************************************
   File utility functions
 ***************************************************************************/
bool fileExists(std::string filename);
int getFiles(std::list<std::string>* fileList, std::string filename, bool recursive=true, std::list<std::string>* extList = 0); // return 0 if files were found??
bool hasExtension(std::string filename, std::list<std::string>* extList);
int getSubdirs(std::list<std::string>* fileList, std::string dirName);
std::string extractFileName(std::string file, bool withExtension = true);
std::string extractPath(std::string file);

#endif /*FMUTILS_H*/
