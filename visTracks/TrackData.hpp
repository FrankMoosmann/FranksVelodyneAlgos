/*
 * TrackData.hpp
 *
 *  Created on: Mar 20, 2012
 *      Author: moosmann
 */
#ifndef TRACKDATA_H_
#define TRACKDATA_H_

#include <vector>
#include <istream>
#include <boost/array.hpp>

template <typename T, size_t N>
T norm(boost::array<T,N> arr) {
  T sum = 0;
  for (size_t i=0; i<N; ++i)
    sum += arr[i]*arr[i];
  return sqrt(sum);
}

class TrackData {
public:
  struct Surface {
    boost::array<float,3> pt;
    boost::array<float,3> nr;
    float nc;

    friend std::istream &operator>> (std::istream &istr, Surface &s);
  };

  TrackData(std::string basic, std::string boundingbox, bool loadSurfaces);
  virtual ~TrackData();

  unsigned long long        uid;
  unsigned int              age;
  unsigned int              lastUp;
  float                     colorR;
  float                     colorG;
  float                     colorB;
  boost::array<float,12>    state; //!< [rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot]
  boost::array<float,12*12> covar;
  unsigned int              nbPts;
  std::vector<Surface>      ptsVIP;
  std::vector<Surface>      ptsADD;
  boost::array<float,3>     center;
  float                     bbAngle;
  boost::array<float,3>     bbExtensions;

  size_t getMemUsage() const;

private:
  void throwexcept(unsigned int state);
};

#endif /* TRACKDATA_H_ */
