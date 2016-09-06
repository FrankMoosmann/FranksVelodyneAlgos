/***************************************************************************
                  FMUtils.h  -  some small helper methods
                             -------------------
    begin                : 04.11.2005
    copyright            : (C) 2005 by Frank Moosmann
    email                : frank.moosmann@inrialpes.fr
 ***************************************************************************/

#include "FMUtils.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <fstream>

using namespace std;

int CurDEBUGLEVEL;

/***************************************************************************
   String functions
 ***************************************************************************/
// decompose a string into substrings
void explode(string s, string e, vector<string>* ret) {
 int iPos = s.find(e, 0);
 int iPit = e.length();
 while(iPos>-1) {
   if(iPos!=0)
     ret->push_back(s.substr(0,iPos));
   s.erase(0,iPos+iPit);
   iPos = s.find(e, 0);
 } // end while
 if(s!="")
   ret->push_back(s);
 return;
}

void trim(string &s) {
  string whitespaces (" \t\f\v\n\r");
  size_t found = s.find_last_not_of(whitespaces);
  if (found!=string::npos)
    s.erase(found+1);
}

/***************************************************************************
   Function for easier output of debug information
   set DEBUG_LVL to control amount of output
 ***************************************************************************/
void debug_out(const int debug_level, const char* message)
{
  if (debug_level <= CurDEBUGLEVEL)
  {
    for (int i=0; i<debug_level; i++) printf("  ");
    printf("%s",message);
    printf("\n");
  }
}

void debug_out(const int debug_level, const char* message, double value1, double value2)
{
  if (debug_level <= CurDEBUGLEVEL)
  {
    for (int i=0; i<debug_level; i++) printf("  ");
    printf(message, value1, value2);
    printf("\n");
  }
}

void debug_out(const int debug_level, const char* message, int value1, int value2)
{
  if (debug_level <= CurDEBUGLEVEL)
  {
    for (int i=0; i<debug_level; i++) printf("  ");
    printf(message, value1, value2);
    printf("\n");
  }
}


/***************************************************************************
   Initialisation routine for the random generator
 ***************************************************************************/
void randomInit()
{
  srand((unsigned)time(0));
}

/***************************************************************************
   Function to get a random number within a specific range
 ***************************************************************************/
unsigned int randomInRange(const unsigned int lowest, const unsigned int highest)
{
  unsigned int range = (highest-lowest)+1;
  double factor = (double)(rand()%RAND_MAX)/(double)(RAND_MAX);
  unsigned int result = lowest + (unsigned int)(factor*(double)range);
  return result;
}

double randomInRange(const double lowest, const double highest)
{
  return (double)rand() * (highest - lowest) / (double)RAND_MAX + lowest;
}

double randomGauss(double mean, double variance)
{
  /*
  double u1,u2,v;
  v = 1;
  while (v >= 1) {
    u1 = rand()/RAND_MAX; //randomInRange(0,1);
    u2 = rand()/RAND_MAX; //randomInRange(0,1);
    v = pow((2*u1 - 1),2) + pow((2*u2 - 1),2);
  }
  double x01 = (2*u1 - 1)*sqrt(-(2*log(v)/v));
  return mean + sqrt(variance)*x01;
  */
  // calclulation method of Box-Mueller:
  double y1 = randomInRange(0.0,0.99999999);
  double y2 = randomInRange(0.0,0.99999999);
  return (sqrt( -2.0 * log(y1)) * cos(2.0 * M_PI * y2)) * sqrt(variance) + mean;
}


/***************************************************************************
   check if file exists
 ***************************************************************************/
bool fileExists(string filename)
{
  bool result = false;
  fstream fin;
  fin.open(filename.c_str(),ios::in);
  if (fin.is_open())
    result=true;
  else
    fin.clear(ios::failbit);
  fin.close();
  return result;
}

/***************************************************************************
   get all files in directory and subdirectories
 ***************************************************************************/
int getFiles(list<string>* fileList, string filename, bool recursive, list<string>* extList)
{
  struct stat dir_stat;   // to retrieve the file type (file/directory)
  stat(filename.c_str(), &dir_stat);
  if (S_ISDIR(dir_stat.st_mode)) { 
    // get files in this directory and recall this function recursively
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(filename.c_str())) == NULL) {
        cerr << "Error(" << errno << ") opening " << filename << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {
      //filter . and ..
      string currfile = string(dirp->d_name);
      if ((currfile == ".") || (currfile == ".."))
        continue;
      string dirAndFile = filename + currfile;
      stat(dirAndFile.c_str(), &dir_stat);
      if (S_ISDIR(dir_stat.st_mode)) {
        if (recursive) {
          dirAndFile = dirAndFile + "/";
          getFiles(fileList, dirAndFile);
        }
      } else
      if (fileExists(dirAndFile)) { 
        if (hasExtension(dirAndFile, extList))
          fileList->push_back(dirAndFile);
      }
    }
    closedir(dp);
    return 0;
  }
  if (fileExists(filename)) { 
    if (hasExtension(filename, extList))
      fileList->push_back(filename);
    return 0;
  }
  return 1;
}

bool hasExtension(string filename, list<string>* extList) {
  bool result = true;
  if (extList != 0 ) {
    result = false;
    for (list<string>::iterator p=extList->begin(); p!=extList->end(); ++p) {
      if (filename.find(*p) != string::npos)
        result = true;
    }
  }
  return result;
}

/***************************************************************************
   get all subdirectories
 ***************************************************************************/
int getSubdirs(list<string>* fileList, string dirName)
{
  struct stat dir_stat;   // to retrieve the file type (file/directory)
  stat(dirName.c_str(), &dir_stat);
  if (S_ISDIR(dir_stat.st_mode)) { 
    // get files in this directory
    DIR *dp;
    struct dirent *dirp;
    if((dp = opendir(dirName.c_str())) == NULL) {
        cerr << "Error(" << errno << ") opening " << dirName << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {
      //filter . and ..
      string currfile = string(dirp->d_name);
      if (currfile == ".")
        continue;
      if (currfile == "..")
        continue;
      string dirAndFile = dirName + currfile;
      stat(dirAndFile.c_str(), &dir_stat);
      if (S_ISDIR(dir_stat.st_mode)) { 
        fileList->push_back(currfile);
      }
    }
    closedir(dp);
    return 0;
  }
  else  { 
    cerr << "Error opening " << dirName << endl;
    return 1;
  }
}

/***************************************************************************
   extract the file name
 ***************************************************************************/
string extractFileName(string file, bool withExtension)
{
  string result;
  size_t pos = file.rfind("/");
  if (pos == string::npos)
    result = file;
  else {
    result = file.substr(pos+1);
  }
  if (!withExtension) {
    pos = result.rfind(".");
    if (pos != string::npos)
      result = result.substr(0,pos);
  }
  return result;
}

/***************************************************************************
   extract the file name
 ***************************************************************************/
string extractPath(string file)
{
  size_t pos = file.rfind("/");
  if (pos == string::npos)
    return string("");
  else {
    return file.substr(0,pos+1); //including the "/" character
  }
}

// end
