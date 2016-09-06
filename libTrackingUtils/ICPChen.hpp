/*!
    \file   ICPChen.h
    \brief  Provides an implementation of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_CHEN_H_
#define ICP_CHEN_H_

#include <iostream>

#include "ICP.hpp"

namespace ICP {

namespace ublas = boost::numeric::ublas;
using namespace matrixTools;

/*!
 * \class ICPChen
 *
 * \brief Iterative Closest Point algorithm after [Chen91]
 *
 * Minimizes Point-Plane distance
 * In case the correspondence search throws and weighting is disabled, the registration result will be poor.
 *
 * Any correspondence function can be used for the ICP:
 *   KdTreeNeighborSearch kdTree(scene);
 *   ICPChen::CorrespondanceFunction corFunc = boost::bind(&KdTreeNeighborSearch::findClosestNeighbor, &kdTree, _1, _2, _3, _4);
 *   ICPChen icp(corWeightFunc,...);
 *   icp.iterate(5,0.001);
 */
template <typename DVectorInputIterator, typename SearchResultsT>
class ICPChen : public ICPBase
{
public:
  typedef typename ICP::NeighborSearch<SearchResultsT>::ncpSearch corrFuncT;

  ICPChen(DVectorInputIterator ptBegin, DVectorInputIterator ptEnd, DVectorInputIterator nBegin, DVectorInputIterator nEnd, corrFuncT corrFunc, unsigned int nCorr, double regCoeff);
  //ICPChen(cpnwSearch corrFunc, DVectorInputIterator ptBegin, DVectorInputIterator ptEnd, DVectorInputIterator nBegin, DVectorInputIterator nEnd);
  ~ICPChen();

  virtual void reset(const DMatrix &R, const DVector &t); //!< reset ICP to a given initial transformation
  virtual bool iterate(); //!< compute next iteration of ICP estimate
  void getEstimatedVariance(DMatrix &covar) const; //!< get the covariance of the current 6-DOF parameter vector
  DMatrix getEstimatedVariance() const; //!< get the covariance of the current 6-DOF parameter vector

//  static void test();

private:

  void init();
  inline double sign(double a) {if (a<0) return -1; if (a>0) return 1; return 0;};

  corrFuncT searchCorrespondence;
  unsigned int nbCorresp;
  const double regCoeff;
  bool useWeights;
  const double COVAR_MAX;

  // input variables
  DVectorInputIterator ptBegin;
  DVectorInputIterator ptEnd;
  DVectorInputIterator nBegin;
  DVectorInputIterator nEnd;
  unsigned int nbModelPts; // = n

  // additional output variables to ICPBase:
  DMatrix covar; // covariance (6x6) of 6-DOF transformation parameters (roll,yaw,pitch,tx,ty,tz)

  // internal variables (allocate only once)
  DVector  center; // size (3 x 1), point that is used as virtual center, so calculations are numerically more robust
  DCMatrix transModelPts; // size (3 x n)
  DCMatrix transModelNrm; // size (3 x n)
  DVector  weights; // size n <-- using a boost-band-matrix and prod() is really slow, thus use a vector
  DCMatrix A;  // data matrix of linear energy model
  DCMatrix At;// transposed data matrix of linear energy model
  DVector  b; // right hand side of linear energy model
  DVector  varN; // additional residual variance due to normal estimation error (used for covariance calculation)
  DCMatrix AtA;  // multiplied data matrix
  DVector  Atb; // multiplied b vector
  DMatrix  Ri; // rotation estimation of current iteration (3x3)
  DVector  ti; // translation estimation of current iteration (3)
};

} // namespace ICP

#include "ICPChen.tcc"

#endif /*ICP_CHEN_H_*/
