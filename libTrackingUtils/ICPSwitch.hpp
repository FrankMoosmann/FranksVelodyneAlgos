/*!
    \file   ICPSwitch.h
    \brief  Provides an implementation of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_SWITCH_H_
#define ICP_SWITCH_H_

#include <iostream>

#include "ICP.hpp"
#include "ICPLinearized.hpp"

namespace ICP {

  /*!
   * \class ICPSwitch
   *
   * \brief Iterative Closest Point algorithm for linearized energy functions
   *
   * Uses any linearized energy functions on points+normals
   * Provides uncertainty measures for retrieved estimate
   * Switches between two energy functions after a fixed number of iterations
   * On switch, a user-given function is called
   *
   * Any correspondence function can be used for the ICP:
   *   KdTreeNeighborSearch kdTree(scene);
   *   ICPSwitch::CorrespondanceFunction corFunc = boost::bind(&KdTreeNeighborSearch::findClosestNeighbor, &kdTree, _1, _2, _3, _4);
   *   ICPSwitch icp(corWeightFunc,...);
   *   icp.iterate(5,0.001);
   */
  template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
  class ICPSwitch : public ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>
  {
  public:
    typedef ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT> ICPLin;
    typedef typename ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::EstimationMode EstimationMode;
    typedef typename ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::corrFuncT corrFuncT;
    typedef boost::function<void(bool)> switchFunc;

    ICPSwitch(DVectorInputIterator ptBegin, DVectorInputIterator ptEnd,
                  DVectorInputIterator nBegin, DVectorInputIterator nEnd,
                  DInputIterator wBegin, DInputIterator wEnd, // weights, can use (double*)NULL if not available
                  DInputIterator nConfBegin, DInputIterator nConfEnd, // can use (double*)NULL if not available
                  DInputIterator nStdDevBegin, DInputIterator nStdDevEnd,
                  DMatrixInputIterator pCovarBegin, DMatrixInputIterator pCovarEnd,
                  DMatrixInputIterator nCovarBegin, DMatrixInputIterator nCovarEnd,
                  corrFuncT corrFunc, EnergyFunction &energy1, EnergyFunction &energy2,
                  unsigned int nbIterSwitch, switchFunc sFunc, bool sFuncInitVal,
                  unsigned int nbCorresp, double regularizationCoeff, EstimationMode estMode);
    ~ICPSwitch();

    virtual void reset(const DMatrix &R, const DVector &t); //!< reset ICP to a given initial transformation
    virtual bool iterate(); //!< compute next iteration of ICP estimate

    // internal variables (allocate only once)
    EnergyFunction *energy1;
    EnergyFunction *energy2;
    unsigned int nbIterSwitch;
    switchFunc sFunc;
    bool sFuncInitVal;
    unsigned int currIter;
  };

} // namespace ICP

#include "ICPSwitch.tcc"

#endif /*ICP_SWITCH_H_*/
