/*!
    \file   ICPBesl.h
    \brief  Provides an implementation of the popular ICP algorithm
    \author  Oliver Pink (<pink@mrt.uka.de>),
            Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2008-2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_BESL_H_
#define ICP_BESL_H_

#include <iostream>

#include "ICP.hpp"

namespace ICP {

namespace ublas = boost::numeric::ublas;
using namespace matrixTools;

/*!
 * \class ICPBesl
 *
 * \brief Iterative Closest Point algorithm after [Besl92]
 *
 * Minimizes point-to-point distances
 *
 * Any correspondence function can be used for the ICP:
 *   KdTreeNeighborSearch kdTree(scene);
 *   ICPBesl::CorrespondanceFunction corFunc = boost::bind(&KdTreeNeighborSearch::findClosestNeighbor, &kdTree, _1, _2, _3, _4);
 *   ICPBesl icp(corFunc);
 *   icp.iterate(5,0.001);
 */
template <typename DVectorInputIterator, typename SearchResultsT>
class ICPBesl : public ICPBase
{
public:
  typedef typename ICP::NeighborSearch<SearchResultsT>::ncpSearch corrFuncT;

  ICPBesl(DVectorInputIterator ptBegin, DVectorInputIterator ptEnd, DVectorInputIterator nBegin, DVectorInputIterator nEnd, corrFuncT corrFunc, unsigned int nCorr);
  ~ICPBesl();

  virtual bool iterate(); //!< compute next iteration of ICP estimate
//  virtual void reset(const DMatrix &R, const DVector &t); //!< reset ICP to a given initial transformation

private:

  void init();

  bool FindCorrespondences();
  bool ComputeTransformation();

  corrFuncT searchCorrespondence;
  unsigned int nbCorresp;
  bool useWeights;

  // Input variables:
  DVectorInputIterator ptBegin;
  DVectorInputIterator ptEnd;
  DVectorInputIterator nBegin;
  DVectorInputIterator nEnd;
  unsigned int ptCount;

  // Output variables in ICPBase

  // Inner variables:
  DCMatrix  transformedModel;  // R+t applied to the model (size 3 x m)
  DCMatrix  worldCorresp;  // Corresponding Scene Points (size 3 x m*c)
  DVector  weights; // size m*c <-- using a boost-band-matrix and prod() is really slow, thus use a vector
  DVector   modelCenter; // Weighted center of Model
  DVector   worldCenter; // Weighted center of matched world points
  DCMatrix  centeredWorld;  // Centered matched Scene Points (size 3 x m*c)
  DCMatrix  centeredModel;  // Centered template (Image) (size 3 x m*c)
  double    cumWeight; // Accumulative weight count

  DCMatrix  tmp; // temporarily used
  DCMatrix  K; // temporarily used: weighted covariance matrix
  DVector   A; // temporarily used for SVD-decomposition of K
  DCMatrix  U; // temporarily used for SVD-decomposition of K
  DCMatrix  Vt; // temporarily used for SVD-decomposition of K
};

} // namespace ICP

#include "ICPBesl.tcc" // template implementation

#endif /*ICP_BESL_H_*/

