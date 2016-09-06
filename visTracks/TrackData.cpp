/*
 * TrackData.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: moosmann
 */

#include "TrackData.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

using namespace std;

std::istream &operator>> (std::istream &istr, TrackData::Surface &s) {
  istr >> s.pt[0];
  istr >> s.pt[1];
  istr >> s.pt[2];
  istr >> s.nr[0];
  istr >> s.nr[1];
  istr >> s.nr[2];
  istr >> s.nc;
  return istr;
};

TrackData::TrackData(std::string basicLine, std::string bbLine, bool loadSurfaces) {
  istringstream istr(basicLine);
  unsigned char dummy;
  size_t nbVIP;

  if (istr.eof()) throwexcept(1);
  istr >> dummy; //"["
  istr >> uid;
  istr >> dummy; //"]"
  istr >> age;
  istr >> lastUp;
  istr >> colorR;
  istr >> colorG;
  istr >> colorB;
  for (size_t i=0; i<state.size(); ++i)
    istr >> state[i];
  for (size_t i=0; i<covar.size(); ++i)
    istr >> covar[i];
  istr >> nbVIP;
  istr >> nbPts;
  if ((loadSurfaces) && (nbPts > 0)) {
    ptsVIP.resize(nbVIP);
    ptsADD.resize(nbPts-nbVIP);
    if (istr.eof()) throwexcept(2);
    try {
      for (size_t i=0; i<ptsVIP.size(); ++i) {
        istr >> ptsVIP[i];
      }
      for (size_t i=0; i<ptsADD.size(); ++i) {
        istr >> ptsADD[i];
      }
    } catch (exception) {
      ptsVIP.resize(0);
      ptsADD.resize(0);
      throwexcept(3);
    }
  }

  try {
    istringstream istrBB(bbLine);
    istrBB >> center[0];
    istrBB >> center[1];
    istrBB >> center[2];
    istrBB >> bbAngle;
    istrBB >> bbExtensions[0];
    istrBB >> bbExtensions[1];
    istrBB >> bbExtensions[2];
    //cout << "loaded BB: " << bbExtensions[0] << "/" << bbExtensions[1] << "/" << bbExtensions[2] << flush;
  } catch (exception) {
    center.assign(0);
    bbAngle = 0;
    bbExtensions.assign(0);
  }
}

TrackData::~TrackData() {
}

size_t TrackData::getMemUsage() const {
  size_t s = 0;
  s += sizeof(uid);
  s += sizeof(age);
  s += sizeof(lastUp);
  s += sizeof(colorR);
  s += sizeof(colorG);
  s += sizeof(colorB);
  s += sizeof(state);
  s += sizeof(covar);
  s += sizeof(Surface)*ptsVIP.size();
  s += sizeof(Surface)*ptsADD.size();
  s += sizeof(center);
  return s;
}

void TrackData::throwexcept(unsigned int state) {
  (void)state;
  string msg = "failed reading TrackInfo from string, status " + boost::lexical_cast<std::string>(state);
  throw runtime_error(msg);
}
