/*!
    \file   ICPPlanePlane.h
    \brief  Provides an implementation of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_PLANEPLANE_H_
#define ICP_PLANEPLANE_H_

#include <iostream>

#include "ICP.hpp"

namespace ICP {

namespace ublas = boost::numeric::ublas;
using namespace matrixTools;

/*!
 * \class ICPPlanePlane
 *
 * \brief A variant of the Iterative Closest Point algorithm
 *
 * Minimizes Plane-Plane distances
 *
 * Any correspondence function can be used for the ICP:
 *   KdTreeNeighborSearch kdTree(scene);
 *   ICPPlanePlane::CorrespondanceFunction corFunc = boost::bind(&KdTreeNeighborSearch::findClosestNeighbor, &kdTree, _1, _2, _3, _4);
 *   ICPPlanePlane icp(corWeightFunc,...);
 *   icp.iterate(5,0.001);
 */
template <typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
class ICPPlanePlane : public ICPBase
{
public:
  typedef typename ICP::NeighborSearch<SearchResultsT>::ncpSearch corrFuncT;

  ICPPlanePlane(DVectorInputIterator ptBegin, DVectorInputIterator ptEnd,
      DVectorInputIterator nBegin, DVectorInputIterator nEnd,
      DInputIterator nConfBegin, DInputIterator nConfEnd, // can use (double*)NULL if not available
      DInputIterator nStdDevBegin, DInputIterator nStdDevEnd,
      corrFuncT corrFunc, unsigned int nCorr, double regCoeff);
  ~ICPPlanePlane();

  virtual bool iterate(); //!< compute next iteration of ICP estimate

private:
  // input variables
  DVectorInputIterator ptBegin;
  DVectorInputIterator ptEnd;
  DVectorInputIterator nBegin;
  DVectorInputIterator nEnd;
  DInputIterator nConfBegin;
  DInputIterator nConfEnd;
  DInputIterator nStdDevBegin;
  DInputIterator nStdDevEnd;
  unsigned int nbModelPts; // = n

  corrFuncT searchCorrespondence;
  unsigned int nbCorresp;
  const double regCoeff;
  bool useWeights;

  // additional output variables to ICPBase:
  DMatrix covar; // covariance (6x6) of 6-DOF transformation parameters (roll,yaw,pitch,tx,ty,tz)

  // internal variables (allocate only once)
  DCMatrix transModelPts; // size (3 x n)
  DCMatrix transModelNrm; // size (3 x n)
  DCMatrix CorrPts; // size (3 x n*nbCorresp)
  DCMatrix CorrNrm; // size (3 x n*nbCorresp)
  DVector  CorrWeights; // size n
  DVector  CorrNConf; // size n
  DVector  CorrNStdDevRAD; // size n
  DCMatrix H; // data matrix for rotation estimation
  DCMatrix U,Vt; // for SVD decomposition of H
  DVector  Dv; // for SVD decomposition of H
  DCMatrix A;  // data matrix for translation estimation
  DCMatrix At;// transposed data matrix
  DVector  b; // right hand side for translation estimation
  DCMatrix AtA;  // multiplied data matrix
  DVector  Atb; // multiplied b vector
  DMatrix  Ri; // rotation estimation of current iteration (3x3)
  DVector  ti; // translation estimation of current iteration (3)
  double   cumWeight; // summed weights of correspondences
};

} // namespace ICP

#include "ICPPlanePlane.tcc"

#endif /*ICP_PLANEPLANE_H_*/
