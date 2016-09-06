/*!
    \file   ICP.h
    \brief  Provides a base for implementations of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2008-2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_H_
#define ICP_H_

#include <cfloat>
#include <boost/function.hpp>
#include "MatrixDefs.hpp"

namespace ICP {

  namespace ublas = boost::numeric::ublas;
  using namespace matrixTools;

  // The following function-types are used to find correspondences: for a given point+normal, the function returns the point(+normal)-correspondence and its distance (and weight)
  // All functions must either return valid results or throw an exception
  typedef boost::function<void (const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist)>
    cpSearch; //!< returns corresponding point and distance (in case returned dist == DBL_MAX correspondence was rejected)
  typedef boost::function<void (const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist, double &weight)>
    cpwSearch; //!< returns corresponding point, weight, and distance (in case returned dist == DBL_MAX or weight == 0.0 correspondence was rejected)
  typedef boost::function<void (const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, DVector &correspNormal, double &dist)>
    cpnSearch; //!< returns corresponding point, normal, and distance (in case returned dist == DBL_MAX correspondence was rejected)
  typedef boost::function<void (const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, DVector &correspNormal, double &dist, double &weight)>
    cpnwSearch; //!< returns corresponding point, normal, weight, and distance (in case returned dist == DBL_MAX  or  weight == 0.0 correspondence was rejected)
  template<class SearchResultsType>
  struct NeighborSearch // typedef-template works only by putting it into a template-struct
  {
    typedef boost::function<SearchResultsType (const DVector &searchPoint, const DVector &searchNormal, double normalConfidence, unsigned int nnCount)>
      ncpSearch;
  };

  // Converting Functions: Usage (e.g.):
  // ICP::cpnwSearch cpnwFunc = boost::bind(&ProjectedNeighborSearch::findClosestProjectedNeighbor, &pnSearch, _1, _2, _3, _4, _5, _6);
  // ICP::cpwSearch corWeightFunc = boost::bind(&ICP::cpnwSearch2cpwSearch, _1, _2, _3, _4, _5, cpnwFunc);
  // - "upgrade":
  void cpSearch2cpwSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &corresp, double &dist, double &weight, cpSearch cfunc, double thresh); //!< set weight=1 if dist<thresh, weight=thresh/dist otherwise
  void cpSearch2cpwSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &corresp, double &dist, double &weight, cpSearch cfunc); //!< set weight = 1
  void cpnSearch2cpnwSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, DVector &correspNormal, double &dist, double &weight, cpnSearch cfunc); //!< set weight = 1
  // - "downgrade":
  void cpnSearch2cpSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist, cpnSearch cfunc);
  void cpnwSearch2cpwSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist, double &weight, cpnwSearch cfunc);
  void cpnwSearch2cpnSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, DVector &correspNormal, double &dist, cpnwSearch cfunc);
  void cpwSearch2cpSearch(const DVector &searchPoint, const DVector &searchNormal, DVector &correspPoint, double &dist, cpwSearch cfunc);


  /*!
   * \class ICP
   *
   * \brief Base class for the Iterative Closest Point algorithm
   *
   * Aligns two point clouds by optimizing rotation and translation
   * At each iteration the following steps are carried out:
   * 1) Current transformation is applied to model-points
   * 2) For each model point correspondences (points) are retrieved from world
   * 3) Transformation is calculated by minimizing some error function
   * 4) Calculated transformation is then combined with current transformation
   */
  class ICPBase {
  public:
    ICPBase();
    virtual ~ICPBase() {};

    void reset(); //!< reset ICP to a zero transformation (will call the R/t-virtual function below)
    void reset(const DVector &transParam); //!< reset ICP to a given initial transformation (rx-roll,ry-pitch,rz-yaw,tx,ty,tz) (will call the R/t-virtual function below)
    void reset(const DMatrix &htm); //!< reset ICP to a given initial transformation given as Homogeneous Matrix Transform
    virtual void reset(const DMatrix &R, const DVector &t); //!< reset ICP to a given initial transformation
    virtual bool iterate() = 0; //!< compute next iteration of ICP estimate. Returns true if iteration was successful. To be implemented by derived classes
    bool iterate(unsigned int maxIter, double minDist); //!< iterate until either maxIter is exceeded or minDist is undershot
    bool iterateUntilConvergence(unsigned int minIter = 3, unsigned int maxIter = UINT_MAX, double minDist = 0.0); //!< iterate until convergence or until either maxIter is exceeded or minDist is undershot

    void getEstimate(DMatrix &R, DVector &t) const; //!< get the current estimate as R+t (3x3 rotation matrix + 3x1 translation vector) a model points p is transformed via this R,t by R*p+t to obtain the new transformed point
    void getEstimate(DMatrix &htm) const; //!< get the current estimate as HTM (4x4 homogeneous transformation matrix). a model points p is transformed via this HTM by HTM*p to obtain the new transformed point
    DMatrix getHTMEstimate() const; //!< get the current estimate as HTM (4x4 homogeneous transformation matrix). a model points p is transformed via this HTM by HTM*p to obtain the new transformed point
    void getEstimate(DVector &transParam) const; //!< get the current 6-DOF parameter vector (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
    DVector getParamEstimate() const; //!< get the current 6-DOF parameter vector (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
    double getLastAvgDist() const {return avgDist;}; //!< get averaged distance of correspondences (before, not after refinement!!)
    double getLastAvgWeightDist() const {return wAvgDist;}; //!< get weighted averaged distance of correspondences (before, not after refinement!!)
    double getLastAvgError() const {return avgErr;}; //!< get averaged error (before, not after refinement!!)
    double getLastAvgWeightError() const {return wAvgErr;}; //!< get weighted averaged error (before, not after refinement!!)

  protected:
    // output variables, to be set by derived classes:
    DMatrix R; // rotation estimation (3x3)
    DVector t; // translation estimation (3)
    double  avgDist; // average distance, not weighted
    double  wAvgDist; // weighted average distance
    double  avgErr; // average error, not weighted
    double  wAvgErr; // weighted average error
  };


} // namespace ICP

#endif /* ICP_H_ */
