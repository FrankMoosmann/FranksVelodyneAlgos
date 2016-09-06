/*!
    \file   ICP.cpp
    \brief  Provides a base for implementations of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2008-2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/

#include "ICP.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "HomogeneousTransformationMatrix.hpp"

using namespace std;

// "upgrade" functions:
void ICP::cpSearch2cpwSearch(const ICP::DVector &searchPoint, const ICP::DVector &searchNormal, ICP::DVector &corresp, double &dist, double &weight, cpSearch cfunc, double thresh)
{
  cfunc(searchPoint, searchNormal, corresp, dist);
  weight = (dist > thresh) ? thresh/dist : 1.0;
}
void ICP::cpSearch2cpwSearch(const ICP::DVector &searchPoint, const ICP::DVector &searchNormal, ICP::DVector &corresp, double &dist, double &weight, cpSearch cfunc)
{
  cfunc(searchPoint, searchNormal, corresp, dist);
  weight = 1.0;
}
void ICP::cpnSearch2cpnwSearch(const ICP::DVector &searchPoint, const ICP::DVector &searchNormal, ICP::DVector &correspPoint, ICP::DVector &correspNormal, double &dist, double &weight, cpnSearch cfunc)
{
  cfunc(searchPoint, searchNormal, correspPoint, correspNormal, dist);
  weight = 1.0;
}
// "downgrade" functions:
void ICP::cpnSearch2cpSearch(const ICP::DVector &searchPoint, const ICP::DVector &searchNormal, ICP::DVector &correspPoint, double &dist, cpnSearch cfunc)
{
  ICP::DVector correspNormal;
  cfunc(searchPoint, searchNormal, correspPoint, correspNormal, dist);
}
void ICP::cpnwSearch2cpwSearch(const ICP::DVector &searchPoint, const ICP::DVector &searchNormal, ICP::DVector &correspPoint, double &dist, double &weight, cpnwSearch cfunc)
{
  ICP::DVector correspNormal;
  cfunc(searchPoint, searchNormal, correspPoint, correspNormal, dist, weight);
}
void ICP::cpnwSearch2cpnSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, DVector &correspNormal, double &dist, cpnwSearch cfunc)
{
  double weight;
  cfunc(searchPoint, searchNormal, correspPoint, correspNormal, dist, weight);
}
void ICP::cpwSearch2cpSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist, cpwSearch cfunc)
{
  double weight;
  cfunc(searchPoint, searchNormal, correspPoint, dist, weight);
}


ICP::ICPBase::ICPBase()
{
  R = ublas::identity_matrix<double>(3);
  t = ublas::zero_vector<double>(3);
  avgDist = DBL_MAX;
  wAvgDist = DBL_MAX;
  avgErr = DBL_MAX;
  wAvgErr = DBL_MAX;
}

void ICP::ICPBase::reset()
{
  reset(DIdMatrix(3), DZeroVector(3)); // virtual -> call of derived class function (if implemented)
}

void ICP::ICPBase::reset(const DVector &x)
{
  if (x.size() != 6)  throw std::out_of_range("ICP::reset: parameter vector must be 6-dimensional");
  DMatrix Rot; DVector trans;
  HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt(x[2],x[1],x[0],x[3],x[4],x[5],Rot,trans);
  reset(Rot, trans); // virtual -> call of derived class function (if implemented)
}

void ICP::ICPBase::reset(const DMatrix &htm)
{
  if (htm.size1() != 4)  throw std::out_of_range("ICP::reset: HTM is not 4x4");
  if (htm.size2() != 4)  throw std::out_of_range("ICP::reset: HTM is not 4x4");
  DMatrix Rot; DVector trans;
  HomogeneousTransformationMatrix::HTM_2_Rt(htm,Rot,trans);
  reset(Rot, trans); // virtual -> call of derived class function (if implemented)
}

void ICP::ICPBase::reset(const DMatrix &rot, const DVector &trans)
{
  if (rot.size1() != 3)  throw std::out_of_range("ICP::reset: Rotation matrix is not 3x3");
  if (rot.size2() != 3)  throw std::out_of_range("ICP::reset: Rotation matrix is not 3x3");
  if (trans.size() != 3)  throw std::out_of_range("ICP::reset: translation vector is not 3-dimensional");
  R = rot;
  t = trans;
  avgDist = DBL_MAX;
  wAvgDist = DBL_MAX;
  avgErr = DBL_MAX;
  wAvgErr = DBL_MAX;
}

bool ICP::ICPBase::iterate(unsigned int maxIter, double minDist)
{
  bool estimateValid = false;
  for (unsigned int i=0; i<maxIter; ++i) {
    estimateValid = iterate();
    if (avgDist <= minDist) break;
  }
  return estimateValid;
}

bool ICP::ICPBase::iterateUntilConvergence(unsigned int minIter, unsigned int maxIter, double minDist)
{
  const double maxRelAvgDistChange(0.01);
  const double maxRelErrorInc(1.2); // allow maximum 20% error increase
  const unsigned int buffSize(3);
  double avgDistChange[buffSize] = {1+maxRelAvgDistChange,1+maxRelAvgDistChange,1+maxRelAvgDistChange};
  unsigned int avgDistIdx = 0;
  double firstAvgDist = getLastAvgWeightDist(); // this error is independent of the energy used
  double oldAvgDist = firstAvgDist;
  double currAvgDist = 0.0;
  bool icpSucceeded = false;
  //cout << " fd=" << firstAvgDist << flush;
  for (unsigned int i=0; i<maxIter; ++i) {
    cout << "." << flush;
    icpSucceeded = iterate();
    currAvgDist = getLastAvgWeightDist();
    if (std::isnan(currAvgDist)) cerr << endl << "ICP:getLastAvgWeightDist()=NAN" << flush;
    avgDistChange[avgDistIdx] = (currAvgDist == 0.0) ? 0.0 : fabs(oldAvgDist-currAvgDist)/(oldAvgDist+currAvgDist);
    //cout << " d=" << currAvgDist << "(delta:" << avgDistChange[avgDistIdx] << ") " << flush;
    if (i+1 >= minIter) {
      if (   (avgDistChange[avgDistIdx] < maxRelAvgDistChange) // if avgDist doesn't change much
          && (avgDistChange[(avgDistIdx-1+buffSize)%buffSize] < maxRelAvgDistChange)
          && (avgDistChange[(avgDistIdx-2+buffSize)%buffSize] < maxRelAvgDistChange)
          && (icpSucceeded)) {
        cout << "successC(i=" << i << ",d="<<currAvgDist<<")" << flush;
        return true;
      }
      if (currAvgDist <= minDist) {
        cout << "successD(i=" << i << ",d="<<currAvgDist<<")" << flush;
        return true;
      }
      if (maxRelErrorInc*firstAvgDist < currAvgDist) { // at most 20% error increase
        cout << "failed(i=" << i << ",d="<<currAvgDist<<")" << flush;
        return false;
      }
    }
    oldAvgDist = currAvgDist;
    avgDistIdx = (avgDistIdx+1)%buffSize;
  }
  return icpSucceeded;
}

void ICP::ICPBase::getEstimate(DMatrix &htm) const
{
  htm = ublas::zero_matrix<double>(4,4);
  ublas::matrix_range<DMatrix> rot(htm, ublas::range (0,3), ublas::range (0,3));
  rot = R;
  htm(0,3) = t(0);
  htm(1,3) = t(1);
  htm(2,3) = t(2);
  htm(3,3) = 1.;
}

matrixTools::DMatrix ICP::ICPBase::getHTMEstimate() const
{
  DMatrix htm = ublas::zero_matrix<double>(4,4);
  ublas::matrix_range<DMatrix> rot(htm, ublas::range (0,3), ublas::range (0,3));
  rot = R;
  htm(0,3) = t(0);
  htm(1,3) = t(1);
  htm(2,3) = t(2);
  htm(3,3) = 1.;
  return htm;
}

void ICP::ICPBase::getEstimate(DMatrix &Rot, DVector &trans) const
{
  Rot = R;
  trans = t;
}

void ICP::ICPBase::getEstimate(DVector &transParam_) const
{
  transParam_.resize(6);
  HomogeneousTransformationMatrix::Rt_2_YawPitchRollXYZ(R, t, transParam_(2), transParam_(1), transParam_(0), transParam_(3), transParam_(4), transParam_(5));
}

matrixTools::DVector ICP::ICPBase::getParamEstimate() const
{
  DVector transParam_(6);
  HomogeneousTransformationMatrix::Rt_2_YawPitchRollXYZ(R, t, transParam_(2), transParam_(1), transParam_(0), transParam_(3), transParam_(4), transParam_(5));
  return transParam_;
}

