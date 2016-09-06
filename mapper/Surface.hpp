#ifndef SURFACE_H_
#define SURFACE_H_
/*
 * Surface.h
 *
 *  Created on: Nov 26, 2010
 *      Author: moosmann
 */

#include <boost/shared_ptr.hpp>
#include <MatrixDefs.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <HomogeneousTransformationMatrix.hpp>


class ScanPose {
  int64_t timestamp;
  matrixTools::DMatrix htm2w;
  matrixTools::DMatrix htm2v;

public:
  ScanPose () : timestamp(0), htm2w(matrixTools::DIdMatrix(4)), htm2v(matrixTools::DIdMatrix(4)) {};
  ScanPose (int64_t t, const matrixTools::DMatrix &h2w) : timestamp(t), htm2w(h2w), htm2v(HomogeneousTransformationMatrix::invert_HTM(h2w)) {};

  int64_t getTimestamp() const {return timestamp;};
  void setTimestamp(int64_t t) {timestamp=t;};
  const matrixTools::DVector getWorldPos() const {matrixTools::DMatrixColConst tHom(htm2w, 3); return matrixTools::ublas::vector_range< matrixTools::DMatrixColConst >(tHom, matrixTools::ublas::range(0,3));};
  const matrixTools::DMatrix& getHTM2w() const {return htm2w;};
  const matrixTools::DMatrix& getHTM2v() const {return htm2v;};
  void setHTM2w(const matrixTools::DMatrix& htm) {htm2w = htm; htm2v = HomogeneousTransformationMatrix::invert_HTM(htm);};
  void setHTM2v(const matrixTools::DMatrix& htm) {htm2v = htm; htm2w = HomogeneousTransformationMatrix::invert_HTM(htm);};

private:
  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    (void)version;
    ar & BOOST_SERIALIZATION_NVP(timestamp);
    ar & BOOST_SERIALIZATION_NVP(htm2w);
    ar & BOOST_SERIALIZATION_NVP(htm2v);
  }
};

class Surface {
public:
  typedef boost::shared_ptr< Surface > SPtr;
  typedef unsigned int PIT;
  PIT poseIdx;
  matrixTools::DVector position;
  matrixTools::DVector origPosition;
  matrixTools::DVector normal;
  unsigned char intensity;
  double normalConfidence;
  double distFromScanner;
  Surface () : position(3), origPosition(3), normal(3) {};
  Surface (PIT pi,                           matrixTools::DVector &p, matrixTools::DVector &n, unsigned char i, double nc, double dfs) : poseIdx(pi), position (p), origPosition(p),  normal(n), intensity(i), normalConfidence(nc), distFromScanner(dfs) {};
  Surface (PIT pi, matrixTools::DVector &po, matrixTools::DVector &p, matrixTools::DVector &n, unsigned char i, double nc, double dfs) : poseIdx(pi), position (p), origPosition(po), normal(n), intensity(i), normalConfidence(nc), distFromScanner(dfs) {};
  inline double operator() (int i) const {return position(i);};
private:
  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_NVP(poseIdx);
    ar & BOOST_SERIALIZATION_NVP(position);
    ar & BOOST_SERIALIZATION_NVP(origPosition);
    ar & BOOST_SERIALIZATION_NVP(normal);
    if (version == 2)
      ar & BOOST_SERIALIZATION_NVP(intensity);
    ar & BOOST_SERIALIZATION_NVP(normalConfidence);
    ar & BOOST_SERIALIZATION_NVP(distFromScanner);
  }
};

// needed for grid-like indexing of surfaces, as shared points on surfaces are used
// without shared pointers, the indexing function could be shifted directly into Surface
class SurfaceProxy {
public:
  Surface::SPtr s;
  unsigned int hitCount;
  SurfaceProxy () : s(), hitCount(0) {};
  SurfaceProxy (Surface::SPtr s_, unsigned int hitCnt=1) : s(s_), hitCount(hitCnt) {};
  inline double operator() (int i) const {return (*s)(i);};
  // is the following really needed???
  inline bool operator== (SurfaceProxy const &other) {return (s == other.s);};
  inline bool operator!= (SurfaceProxy const &other) {return (s != other.s);};
  inline Surface& operator*() {return *s;};
  inline const Surface& operator*() const {return *s;};
  inline Surface* operator->() {return s.get();};
  inline const Surface* operator->() const {return s.get();};

private:
  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    (void)version;
    ar & BOOST_SERIALIZATION_NVP(s);
  }
};


BOOST_CLASS_VERSION(Surface, 2)
BOOST_CLASS_VERSION(SurfaceProxy, 1)

#endif /* SURFACE_H_ */
