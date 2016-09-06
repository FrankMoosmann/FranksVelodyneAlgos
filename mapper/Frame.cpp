#include "Frame.hpp"

#include <cfloat>
#include <climits>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <kogmo_time.h>

using namespace std;
using namespace matrixTools;


Frame::Frame(const unsigned int horizSize, const unsigned int vertSize)
  : distance(horizSize, vertSize)
   ,intensity(horizSize, vertSize)
   ,point3D(horizSize, vertSize)
   ,distDerivativeH(horizSize-1, vertSize)
   ,distDerivativeV(horizSize, vertSize-1)
   ,connectivityH(horizSize-1, vertSize)
   ,connectivityV(horizSize, vertSize-1)
   ,normal3D(horizSize, vertSize)
   ,normalStdDevRAD(horizSize, vertSize)
   ,normalVariance3D(horizSize, vertSize)
   ,normalConfidence(horizSize, vertSize)
   ,surface(horizSize, vertSize)
   ,tmpVec3D(horizSize, vertSize)
   ,tmpFloat1(horizSize, vertSize)
{
  reset();
}

Frame::~Frame()
{
}

void Frame::reset()
{
  valid = false;
  recTime = 0;
  timeDiffToLast = 0;
  positionHTM2v = ublas::identity_matrix<double>(4);
  positionHTM2w = ublas::identity_matrix<double>(4);
  insPositionHTM2w = ublas::identity_matrix<double>(4);
  moveHTM = ublas::identity_matrix<double>(4);
  diffHTMToLast = ublas::identity_matrix<double>(4);
  distance.fill(DBL_MAX); // set externally
  intensity.fill(0); // set externally
  point3D.fill(DVector(3,DBL_MAX)); // by AlgoFrameFeatures
  distDerivativeH.fill(DBL_MAX); // by AlgoFrameFeatures
  distDerivativeV.fill(DBL_MAX); // by AlgoFrameFeatures
  connectivityH.fill(0.0f); // by AlgoFrameFeatures
  connectivityV.fill(0.0f); // by AlgoFrameFeatures
  normal3D.fill(DVector(3,0.0)); // by AlgoFrameFeatures
  normalStdDevRAD.fill(M_PI); // by AlgoFrameFeatures
  normalVariance3D.fill(DMatrix(3,3,0.0)); // by AlgoFrameFeatures
  normalConfidence.fill(0.0); // by AlgoFrameFeatures
  surface.fill(Surface::SPtr());
  tmpVec3D.fill(DVector(3,0.0)); // borders must be initialized with 3D vectors
}


std::string Frame::getTime()
{
  KogniMobil::kogmo_timestamp_string_t str;
  KogniMobil::kogmo_timestamp_to_string(recTime, str);
  return string(str);
}

