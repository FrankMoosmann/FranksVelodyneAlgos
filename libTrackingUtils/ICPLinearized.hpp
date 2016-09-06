/*!
    \file   ICPLinearized.h
    \brief  Provides an implementation of the popular ICP algorithm
    \author  Frank Moosmann (<moosmann@mrt.uka.de>)
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef ICP_LINEARIZED_H_
#define ICP_LINEARIZED_H_

#include <iostream>

#include "ICP.hpp"

namespace ICP {

  /*!
   * \class EnergyFunction
   * \brief Derive from this base class to implement an energy function for linearized ICP below
   */
  class EnergyFunction {
  public:
    struct Surface {
      DVector p;
      DVector n;
//      double pStdDevM;
      double nConfidence;
      double nStdDevRAD;
      DMatrix pCovar3D;
      DMatrix nCovar3D;
      Surface()
        : p(DVector(3)), n(DVector(3)), /*pStdDevM(0.0),*/ nConfidence(0.0), nStdDevRAD(M_PI/4), pCovar3D(DZeroMatrix(3)), nCovar3D(DZeroMatrix(3)) {};
      Surface(const DVector &p_, const DVector &n_, /*const double pStdDevM_,*/ const double nConfidence_,  const double nStdDevRAD_, const DMatrix &pCovar3D_,  const DMatrix &nCovar3D_)
        : p(p_), n(n_), /*pStdDevM(pStdDevM_),*/ nConfidence(nConfidence_), nStdDevRAD(nStdDevRAD_), pCovar3D(pCovar3D_), nCovar3D(nCovar3D_) {};
    };
    EnergyFunction() {};
    virtual ~EnergyFunction() {};
    //! Has to calculate the (non-linear) error of the given surface-correspondence
    virtual double calculate(const Surface &s1, const Surface &s2) const = 0;
    //! Has to calculate the 6x6 hessian matrix and 6d gradient vector for a given surface-correspondence
    virtual void approximateAt(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DVector &gradient) const = 0;
    //! Has to calculate the 6x6 covariance matrix of the 6d transformation parameter for a given surface-correspondence
    virtual void approximateCovar(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DSMatrix &TransMeasCovar) const {(void)s1; (void)s2; Hessian = ublas::zero_matrix<double>(6,6); TransMeasCovar = ublas::zero_matrix<double>(6,6);};
  };

  /*!
   * \class ICPLinearized
   *
   * \brief Iterative Closest Point algorithm for linearized energy functions
   *
   * Uses any linearized energy functions on points+normals
   * Provides uncertainty measures for retrieved estimate
   *
   * Any correspondence function can be used for the ICP:
   *   KdTreeNeighborSearch kdTree(scene);
   *   ICPLinearized::CorrespondanceFunction corFunc = boost::bind(&KdTreeNeighborSearch::findClosestNeighbor, &kdTree, _1, _2, _3, _4);
   *   ICPLinearized icp(corWeightFunc,...);
   *   icp.iterate(5,0.001);
   */
  template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
  class ICPLinearized : public ICPBase
  {
  public:
    typedef typename ICP::NeighborSearch<SearchResultsT>::ncpSearch corrFuncT;
    enum EstimationMode {Mode3R3T, Mode1R3T};

    ICPLinearized(DVectorInputIterator ptBegin, DVectorInputIterator ptEnd,
                  DVectorInputIterator nBegin, DVectorInputIterator nEnd,
                  DInputIterator wBegin, DInputIterator wEnd, // weights, can use (double*)NULL if not available
                  DInputIterator nConfBegin, DInputIterator nConfEnd, // can use (double*)NULL if not available
                  DInputIterator nStdDevBegin, DInputIterator nStdDevEnd,
                  DMatrixInputIterator pCovarBegin, DMatrixInputIterator pCovarEnd,
                  DMatrixInputIterator nCovarBegin, DMatrixInputIterator nCovarEnd,
                  corrFuncT corrFunc, EnergyFunction &energy,
                  unsigned int nbCorresp, double regularizationCoeff, EstimationMode estMode);
    ~ICPLinearized();

    virtual void reset(const DMatrix &R, const DVector &t); //!< reset ICP to a given initial transformation
    virtual bool iterate(); //!< compute next iteration of ICP estimate
    bool establishCorrespondences(); //!< compute only the correspondences
    bool minimizeError(); //!< compute only the new transformation
    void getErrorStatistics(double &avgWError, DSMatrix &measCovar, DSMatrix &hessian) const; //!< get the analytical covariance of the current 6-DOF parameter vector
    DSMatrix getEstimatedVariance() const; //!< get the covariance as inverse hessian of the current 6-DOF parameter vector
    DSMatrix getAnalyticalVariance() const; //!< get the covariance from implicit function theorem of the current 6-DOF parameter vector
    DSMatrix getSampledVariance(bool correspondenceSearch = false); //!< get the sampled covariance of the current 6-DOF parameter vector
    void testSampledVariance(bool correspondenceSearch = false); //!< test the sampled covariance of the current 6-DOF parameter vector
    /*! \brief get the error value for the current ICP estimate
     * \return the summed, weighted error over all correspondences.
     */
    double getError() const;
    double getError(const DMatrix &rot, const DVector &trans) const; //!< return the error value for the given transformation of the original point cloud, error value is as with getError()
    double getAvgWeightError() const {return getError()/cumWeight;}; //!< returns the average error per correspondence for the current ICP estimate
    double getAvgWeightError(const DMatrix &rot, const DVector &trans) const {return getError(rot,trans)/cumWeight;}; //!< returns the average error per correspondence for the given transformation of the original point cloud, error value is as with getError()
    double getAvgWeight() const; //!< returns the average weight
    double getSignificantWeightRatio() const; //!< returns the percentage of correspondences with weight>0.5
    double getOutlierRatio() const {return outlierRatio;}; //!< returns the percentage of correspondences with weight>0.5

    DVector getModelCenter() const {return modelCenter;}; //!< only valid if an iteration was already carried out
    DVector getTransModelCenter() const {return transModelCenter;}; //!< only valid if an iteration was already carried out

  //  static void test();

  //protected: // TODO (9): make protected again!!

    inline double sign(double a) {if (a<0) return -1; if (a>0) return 1; return 0;};

    /*!
     * \class EmpCovarEstimator This inner class allows for empirical covariance estimation
     *
     * Therefore, the desired average error per correspondence must be given
     * (that can be expected due to measurement noise, wrong correspondences, etc)
     * The error function is then sampled, and the error value together with the given expected
     * error used to determine a scaled parameter point that would produce the expected error
     * (based on a quadratic function model). All these scaled parameter points "x" are then
     *  used to calculate a covariance matrix in a normal fashion: Covar = 1/n * SUM x*xT
     */
    class EmpCovarEstimator {
    public:
      EmpCovarEstimator(ICPLinearized &icp, unsigned int nbSamples, bool resetCorresp);
      ~EmpCovarEstimator();
      void generateSample(const DVector &relTrans);
      DSMatrix getCovar(bool normalizeVolume = false);
      double getLastError();
      DVector getLastSample();
    private:
      ICPLinearized &icp;
      unsigned int nb; // number of samples that will be generated
      bool correspondenceSearch; // whether to re-find correspondences for each sample
      DVector x0; // current icp transformation as vector (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
      double b0; // error at current transformation
      double expAvgError; // expected error per correspondence
      double errMin;
      DVector xt; // sampled transformation (total) = (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
      DMatrix Rot; DVector trans; // sampled transformation as Rot-Matrix and trans-Vector
      DMatrix A; // sampled points centered around current transformation, will be scaled by function values in order to create covariance matrix
      unsigned int idx;
      //  DMatrix A(nb,21); // 21 parameters to estimate (6+4+5+3+2+1) of the diagonal matrix
      //  DVector b(nb);
      // temporary variables:
      //double nx; // norm of the sampled transformation (relative to current transformation)
      double b; // error at sampled transformation
      std::vector<double> bs;
      double scale;
    };

    EstimationMode estMode;
    corrFuncT searchCorrespondence;
    EnergyFunction *energy;
    unsigned int nbCorresp;
    const double regCoeff;
    bool useWeights;
    const double COVAR_MAX;
    const double ERROR_MIN;

    // input variables
    DVectorInputIterator ptBegin;
    DVectorInputIterator ptEnd;
    DVectorInputIterator nBegin;
    DVectorInputIterator nEnd;
    DMatrixInputIterator pCovarBegin;
    DMatrixInputIterator pCovarEnd;
    DMatrixInputIterator nCovarBegin;
    DMatrixInputIterator nCovarEnd;
    DInputIterator wBegin;
    DInputIterator wEnd;
    DInputIterator nConfBegin;
    DInputIterator nConfEnd;
    DInputIterator nStdDevBegin;
    DInputIterator nStdDevEnd;
    unsigned int nbModelPts; // = n

    // additional output variables to ICPBase:
    DSMatrix covar; // covariance (6x6) of 6-DOF transformation parameters (roll,yaw,pitch,tx,ty,tz)

    // internal variables (allocate only once)
    unsigned int iteration;
    DVector  center; // size (3 x 1), point that is used as virtual center, so calculations are numerically more robust
    DSMatrix HS;  // sum of hessians of the linearized model, size (6x6)
    DSMatrix Hp;  // hessian of the current point-correspondence, size (6x6)
    DVector  gS; // sum of the gradients of the linearized model, size (6)
    DVector  gp; // gradient of the current point-correspondence, size (6)
    //DVector  varN; // additional residual variance due to normal estimation error (used for covariance calculation)
    DMatrix  Ri; // rotation estimation of current iteration (3x3)
    DVector  ti; // translation estimation of current iteration (3)
    bool     lastIterTransOnly; // flag indicates whether last iteration calculated translation only

    // temporary: store correspondences
    DCMatrix transModelPts; // size (3 x n)
    DCMatrix transModelNrm; // size (3 x n)
    std::vector<DSMatrix> transPCovar; // size n
    std::vector<DSMatrix> transNCovar; // size n
    DVector  modelCenter; // size (3 x 1)
    DVector  transModelCenter; // size (3 x 1)
    DVector  correspCenter; // size (3 x 1)
    DCMatrix CorrPts; // size (3 x n*nbCorresp)
    DCMatrix CorrNrm; // size (3 x n*nbCorresp)
    DVector  CorrDist; // size (3 x n*nbCorresp) correspondence distance
    DVector  CorrErr; // size (3 x n*nbCorresp) correspondence energy errors
    DVector  CorrWeights; // size n
    DVector  CorrNConf; // size n
    DVector  CorrNStdDevRAD; // size n
    std::vector<DSMatrix> CorrPCovar; // size n*nbCorresp
    std::vector<DSMatrix> CorrNCovar; // size n*nbCorresp
    double   cumWeight; // summed weights of correspondences
    double   outlierRatio;
  };

  /*!
   * \class NormalConfLinCombEnergy
   * \brief Energy function that combines two energies depending on the normal confidence
   */
  class NormalConfLinCombEnergy : public ICP::EnergyFunction {
  private:
    const ICP::EnergyFunction &enNormal;
    const ICP::EnergyFunction &enPoint;
    mutable DSMatrix tmpHessian1, tmpHessian2, tmpCovar1, tmpCovar2;
    mutable DVector tmpGradient1, tmpGradient2;
  public:
    NormalConfLinCombEnergy(const EnergyFunction &normalEn, const EnergyFunction &ptEn)
      : EnergyFunction(), enNormal(normalEn), enPoint(ptEn), tmpHessian1(6,6), tmpHessian2(6,6), tmpCovar1(6,6), tmpCovar2(6,6), tmpGradient1(6), tmpGradient2(6) {};
    virtual double calculate(const Surface &s1, const Surface &s2) const {
//      std::cout << "LC: " << s1.p << "  " << s2.p << std::flush;
//      std::cout << " --> En=" << enNormal.calculate(s1,s2) << std::flush;
//      std::cout << "  Ep=" << enPoint.calculate(s1,s2) << std::flush;
//      std::cout << "  w=" << s2.nConfidence << std::endl;
      return s2.nConfidence * enNormal.calculate(s1,s2) + (1.0-s2.nConfidence) * enPoint.calculate(s1,s2);
    }
    virtual void approximateAt(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DVector &gradient) const {
      enNormal.approximateAt(s1, s2, tmpHessian1, tmpGradient1);
      enPoint.approximateAt(s1, s2, tmpHessian2, tmpGradient2);
      gradient = s2.nConfidence * tmpGradient1 + (1.0-s2.nConfidence) * tmpGradient2;
      Hessian = s2.nConfidence * tmpHessian1 + (1.0-s2.nConfidence) * tmpHessian2;
    }
    virtual void approximateCovar(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DSMatrix &TransMeasCovar) const {
      enNormal.approximateCovar(s1, s2, tmpHessian1, tmpCovar1);
      enPoint.approximateCovar(s1, s2, tmpHessian2, tmpCovar2);
      Hessian = s2.nConfidence * tmpHessian1 + (1.0-s2.nConfidence) * tmpHessian2;
      TransMeasCovar = s2.nConfidence * tmpCovar1 + (1.0-s2.nConfidence) * tmpCovar2;
    }
  };

  /*!
   * \class Point2PointEnergy
   * \brief Linearized energy function of the Point-2-Point ICP after Besl
   *
   * Minimizes the error (dT*d), i.e. the square distance vector d=(R*p + t - q)
   */
  class Point2PointEnergy : public ICP::EnergyFunction {
  private:
    mutable DVector d;
    mutable double pxs,pys,pzs;
  public:
    Point2PointEnergy() : EnergyFunction(), d(3) {};
    virtual double calculate(const Surface &s1, const Surface &s2) const {
      //std::cout << "PtPt: " << s1.p << "  " << s2.p << std::flush;
      d = s1.p-s2.p;
      return ublas::inner_prod(d, d);
    }
    virtual void approximateAt(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DVector &gradient) const {
      Hessian.resize(6);
      gradient.resize(6);
      const DVector &p(s1.p);
      const DVector &q(s2.p);
      d = p-q;
      DVectorRange(gradient,ublas::range(3,6)) = d;
      DVectorRange c(gradient,ublas::range(0,3));
      cross_product(p, d, c);
      gradient *= 2;
      pxs = p(0)*p(0);
      pys = p(1)*p(1);
      pzs = p(2)*p(2);
      Hessian(0,0) = 2 * (pys+pzs);
      Hessian(0,1) = -2 * p(0)*p(1);
      Hessian(0,2) = -2 * p(0)*p(2);
      Hessian(0,3) = 0;
      Hessian(0,4) = -2 * p(2);
      Hessian(0,5) = 2 * p(1);
      Hessian(1,1) = 2 * (pxs+pzs);
      Hessian(1,2) = -2 * p(1)*p(2);
      Hessian(1,3) = 2 * p(2);
      Hessian(1,4) = 0;
      Hessian(1,5) = -2 * p(0);
      Hessian(2,2) = 2 * (pxs+pys);
      Hessian(2,3) = -2 * p(1);
      Hessian(2,4) = 2 * p(0);
      Hessian(2,5) = 0;
      Hessian(3,3) = 2;
      Hessian(3,4) = 0;
      Hessian(3,5) = 0;
      Hessian(4,4) = 2;
      Hessian(4,5) = 0;
      Hessian(5,5) = 2;
    }
    virtual void approximateCovar(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DSMatrix &TransMeasCovar) const {
      const DVector &p(s1.p);
      const DVector &q(s2.p);
      double px = p(0);
      double py = p(1);
      double pz = p(2);
      double qx = q(0);
      double qy = q(1);
      double qz = q(2);

      DVector unused_gradient;
      approximateAt(s1, s2, Hessian, unused_gradient);

      DMatrix Ezx = ublas::zero_matrix<double>(6,6); // corresponds to d^2 Err / dz dx with x as below and z as above
      Ezx(1,0) = 2 * qz;
      Ezx(2,0) = -2 * qy;
      Ezx(3,0) = 2;
      Ezx(0,1) = -2 * qz;
      Ezx(2,1) = 2 * qx;
      Ezx(4,1) = 2;
      Ezx(0,2) = 2 * qy;
      Ezx(1,2) = -2 * qx;
      Ezx(5,2) = 2;
      Ezx(1,3) = -2 * pz;
      Ezx(2,3) = 2 * py;
      Ezx(3,3) = -2;
      Ezx(0,4) = 2 * pz;
      Ezx(2,4) = -2 * px;
      Ezx(4,4) = -2;
      Ezx(0,5) = -2 * py;
      Ezx(1,5) = 2 * px;
      Ezx(5,5) = -2;
      DMatrix EzxT = ublas::trans(Ezx);

      // [px, py, pz, nx, ny, nz, qx, qy, qz]
      DSMatrix CovZ = ublas::zero_matrix<double>(6,6); // covariance of the measurements "z", here: [px, py, pz, qx, qy, qz]
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(0,3), ublas::range(0,3)) = s1.pCovar3D;
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(3,6), ublas::range(3,6)) = s2.pCovar3D;

//      std::cout << "§" << Ezx.size1() << "_" << Ezx.size2() << std::flush;
//      std::cout << ";" << CovZ.size1() << "_" << CovZ.size2() << std::flush;
      DMatrix tmp = ublas::prod(Ezx,CovZ);
//      std::cout << "#" << tmp.size1() << "_" << tmp.size2() << std::flush;
//      std::cout << "%" << EzxT.size1() << "_" << EzxT.size2() << std::flush;
      TransMeasCovar = ublas::prod(tmp,EzxT);
    }
  };

  /*!
   * \class Point2PlaneEnergy
   * \brief Energy function of the Point-2-Plane ICP after Chen
   *
   * Minimizes the error (dT*m)^2, i.e. the distance vector d=(R*p + t - q) projected onto the q-surface with normal m
   * Thus, "measurements" are px,py,pz (s1.p) and qx,qy,qz with mx,my,mz (s2.p and s2.n)
   * Gradient (=Ez) corresponds to d Err / dz with z=[px, py, pz, mx, my, mz, qx, qy, qz]
   * Hessian (=Ezz) corresponds to d^2 Err / dz^2 with z=[px, py, pz, mx, my, mz, qx, qy, qz]
   * Ezx corresponds to d^2 Err / dz dx with z=[px, py, pz, mx, my, mz, qx, qy, qz] and x=[rx, ry, rz, tx, ty, tz]
   * TransMeasCovar is thus the measurement covariances transformed into x=[rx, ry, rz, tx, ty, tz]-space
   */
  class Point2PlaneEnergy : public ICP::EnergyFunction {
  private:
    mutable DVector a;
    mutable double b;
  public:
    Point2PlaneEnergy() : EnergyFunction(), a(6) {};
    virtual double calculate(const Surface &s1, const Surface &s2) const {
      return pow(ublas::inner_prod(s1.p-s2.p, s2.n),2);
    }
    virtual void approximateAt(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DVector &gradient) const {
      DVectorRange c(a,ublas::range(0,3));
      cross_product(s1.p, s2.n, c);
      DVectorRange(a,ublas::range(3,6)) = s2.n;
      b = ublas::inner_prod((s1.p-s2.p),s2.n);
//      std::cout << std::endl << "a:" << a << "  b:" << b << std::flush;
      gradient = 2 * b * a;
      Hessian = 2 * ublas::outer_prod(a,a);
    }
    virtual void approximateCovar(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DSMatrix &TransMeasCovar) const {
      const DVector &p(s1.p);
      //const DVector &n(s1.n);
      const DVector &q(s2.p);
      const DVector &m(s2.n);
      double px = p(0);
      double py = p(1);
      double pz = p(2);
      double qx = q(0);
      double qy = q(1);
      double qz = q(2);
      double mx = m(0);
      double my = m(1);
      double mz = m(2);

      DVector unused_gradient;
      approximateAt(s1, s2, Hessian, unused_gradient);
      DMatrix Ezx(6,9);
      Ezx(0,0) = 2 * mx * (-pz * my + py * mz);
      Ezx(1,0) = 2 * mx * (pz * mx - px * mz) - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * mz;
      Ezx(2,0) = 2 * mx * (-py * mx + px * my) + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * my;
      Ezx(3,0) = 2 * mx * mx;
      Ezx(4,0) = 2 * my * mx;
      Ezx(5,0) = 2 * mz * mx;
      Ezx(0,1) = 2 * my * (-pz * my + py * mz) + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * mz;
      Ezx(1,1) = 2 * (pz * mx - px * mz) * my;
      Ezx(2,1) = 2 * (-py * mx + px * my) * my - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * mx;
      Ezx(3,1) = 2 * my * mx;
      Ezx(4,1) = 2 * my * my;
      Ezx(5,1) = 2 * mz * my;
      Ezx(0,2) = 2 * mz * (-pz * my + py * mz) - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * my;
      Ezx(1,2) = 2 * (pz * mx - px * mz) * mz + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * mx;
      Ezx(2,2) = 2 * (-py * mx + px * my) * mz;
      Ezx(3,2) = 2 * mz * mx;
      Ezx(4,2) = 2 * mz * my;
      Ezx(5,2) = 2 * mz * mz;
      Ezx(0,3) = 2 * (px - qx) * (-pz * my + py * mz);
      Ezx(1,3) = 2 * (pz * mx - px * mz) * (px - qx) + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pz;
      Ezx(2,3) = 2 * (-py * mx + px * my) * (px - qx) - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * py;
      Ezx(3,3) = 4 * (px - qx) * mx + 2 * (py - qy) * my + 2 * (pz - qz) * mz;
      Ezx(4,3) = 2 * my * (px - qx);
      Ezx(5,3) = 2 * mz * (px - qx);
      Ezx(0,4) = 2 * (py - qy) * (-pz * my + py * mz) - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pz;
      Ezx(1,4) = 2 * (pz * mx - px * mz) * (py - qy);
      Ezx(2,4) = 2 * (-py * mx + px * my) * (py - qy) + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * px;
      Ezx(3,4) = 2 * (py - qy) * mx;
      Ezx(4,4) = 4 * (py - qy) * my + 2 * (px - qx) * mx + 2 * (pz - qz) * mz;
      Ezx(5,4) = 2 * mz * (py - qy);
      Ezx(0,5) = 2 * (pz - qz) * (-pz * my + py * mz) + 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * py;
      Ezx(1,5) = 2 * (pz * mx - px * mz) * (pz - qz) - 2 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * px;
      Ezx(2,5) = 2 * (-py * mx + px * my) * (pz - qz);
      Ezx(3,5) = 2 * (pz - qz) * mx;
      Ezx(4,5) = 2 * my * (pz - qz);
      Ezx(5,5) = 4 * (pz - qz) * mz + 2 * (px - qx) * mx + 2 * (py - qy) * my;
      Ezx(0,6) = -2 * mx * (-pz * my + py * mz);
      Ezx(1,6) = -2 * mx * (pz * mx - px * mz);
      Ezx(2,6) = -2 * mx * (-py * mx + px * my);
      Ezx(3,6) = -2 * mx * mx;
      Ezx(4,6) = -2 * my * mx;
      Ezx(5,6) = -2 * mz * mx;
      Ezx(0,7) = -2 * my * (-pz * my + py * mz);
      Ezx(1,7) = -2 * (pz * mx - px * mz) * my;
      Ezx(2,7) = -2 * (-py * mx + px * my) * my;
      Ezx(3,7) = -2 * my * mx;
      Ezx(4,7) = -2 * my * my;
      Ezx(5,7) = -2 * mz * my;
      Ezx(0,8) = -2 * mz * (-pz * my + py * mz);
      Ezx(1,8) = -2 * (pz * mx - px * mz) * mz;
      Ezx(2,8) = -2 * (-py * mx + px * my) * mz;
      Ezx(3,8) = -2 * mz * mx;
      Ezx(4,8) = -2 * mz * my;
      Ezx(5,8) = -2 * mz * mz;
      DMatrix EzxT = ublas::trans(Ezx);

      DSMatrix CovZ = ublas::zero_matrix<double>(9,9); // covariance of the measurements "z", here: [px, py, pz, mx, my, mz, qx, qy, qz]
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(0,3), ublas::range(0,3)) = s1.pCovar3D;
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(3,6), ublas::range(3,6)) = s2.nCovar3D;
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(6,9), ublas::range(6,9)) = s2.pCovar3D;
//      std::cout << std::endl << "covar_meas(i):" << std::endl << aligned_write(CovZ) << std::endl;
//      std::cout << std::endl << "Ezx(i):" << std::endl << aligned_write(Ezx) << std::endl;

//      std::cout << "§" << Ezx.size1() << "_" << Ezx.size2() << std::flush;
//      std::cout << ";" << CovZ.size1() << "_" << CovZ.size2() << std::flush;
      DMatrix tmp = ublas::prod(Ezx,CovZ);
//      std::cout << "#" << tmp.size1() << "_" << tmp.size2() << std::flush;
//      std::cout << "%" << EzxT.size1() << "_" << EzxT.size2() << std::flush;
      TransMeasCovar = ublas::prod(tmp,EzxT);
    }
  };

  /*!
   * \class Point2NoisyPlaneEnergy
   * \brief Noisified Point-2-Plane ICP
   *
   * Minimizes the error (dT*n)^2/(dT*d) which corresponds to the angle
   */
  class Point2NoisyPlaneEnergy : public ICP::EnergyFunction {
  private:
    const double expon;
    mutable DVector d,c,k;
    mutable double pxs,pys,pzs,ds,dl,d3,dtm,dtms,tA,tB;
  public:
    Point2NoisyPlaneEnergy(double exponent = 1.0) : EnergyFunction(), expon(exponent), d(3), c(3), k(3) {};
    virtual double calculate(const Surface &s1, const Surface &s2) const {
      d = s1.p-s2.p;
      double dsquare = ublas::inner_prod(d, d);
      if (dsquare < 0.0001) // avoid division by zero. define error of overlapping points as zero
        return 0.0;
//      std::cout << "  " << ublas::inner_prod(d, s2.n) << " / sqrt(" << dsquare << ") / " << s2.nStdDevRAD << std::flush;
//      std::cout << "  -> Enp=" << pow(asin(ublas::inner_prod(d, s2.n)/sqrt(dsquare)),2)/std::max(0.01,s2.nStdDevRAD) << std::endl;
//      return pow(asin(ublas::inner_prod(d, s2.n)/sqrt(dsquare)),2)/std::max(0.01,s2.nStdDevRAD);
      return pow(ublas::inner_prod(d, s2.n),2)/pow(dsquare,expon)/std::max(0.01,s2.nStdDevRAD);
//      return pow(ublas::inner_prod(d, s2.n),2)/dsquare/s2.nStdDevRAD;
    }
    virtual void approximateAt(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DVector &gradient) const {
      Hessian.resize(6);
      gradient.resize(6);
      const DVector &p(s1.p);
      //const DVector &n(s1.n);
      const DVector &q(s2.p);
      const DVector &m(s2.n);

      //   (dLinT, n)^2 / ((dLinT, dLin), sigalph)):
      double alpha = 1.0/std::max(s2.nStdDevRAD,0.01); // minimum 0,5°
      double px = p(0);
      double py = p(1);
      double pz = p(2);
      double qx = q(0);
      double qy = q(1);
      double qz = q(2);
      double mx = m(0);
      double my = m(1);
      double mz = m(2);
      double dspq = (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1));
      double dpq = dspq*dspq;
      gradient(0) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * (-my * pz + mz * py) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py);
      gradient(1) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * (mx * pz - mz * px) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px);
      gradient(2) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * (-mx * py + my * px) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px);
      gradient(3) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * mx - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * px - 0.2e1 * qx);
      gradient(4) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * my - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * py - 0.2e1 * qy);
      gradient(5) = 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) / dspq * mz - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * pz - 0.2e1 * qz);
      Hessian(0,0) = 0.2e1 * alpha * pow(-my * pz + mz * py, 0.2e1) / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py, 0.2e1) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * pz * pz + 0.2e1 * py * py);
      Hessian(0,1) = 0.2e1 * alpha * (mx * pz - mz * px) / dspq * (-my * pz + mz * py) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (mx * pz - mz * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * px * py;
      Hessian(0,2) = 0.2e1 * alpha * (-mx * py + my * px) / dspq * (-my * pz + mz * py) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (-mx * py + my * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * px * pz;
      Hessian(0,3) = 0.2e1 * alpha * mx / dspq * (-my * pz + mz * py) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * mx + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * px - 0.2e1 * qx);
      Hessian(0,4) = 0.2e1 * alpha * my / dspq * (-my * pz + mz * py) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * my + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * pz;
      Hessian(0,5) = 0.2e1 * alpha * mz / dspq * (-my * pz + mz * py) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-my * pz + mz * py) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * mz + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * py;
//      Hessian(1,0) = 0.2e1 * alpha * (nx * pz - nz * px) / dspq * (-ny * pz + nz * py) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-ny * pz + nz * py) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (nx * pz - nz * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * px * py;
      Hessian(1,1) = 0.2e1 * alpha * pow(mx * pz - mz * px, 0.2e1) / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (mx * pz - mz * px) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px, 0.2e1) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * pz * pz + 0.2e1 * px * px);
      Hessian(1,2) = 0.2e1 * alpha * (-mx * py + my * px) / dspq * (mx * pz - mz * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (mx * pz - mz * px) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (-mx * py + my * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * py * pz;
      Hessian(1,3) = 0.2e1 * alpha * mx / dspq * (mx * pz - mz * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (mx * pz - mz * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * mx + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * pz;
      Hessian(1,4) = 0.2e1 * alpha * my / dspq * (mx * pz - mz * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (mx * pz - mz * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * my + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * py - 0.2e1 * qy);
      Hessian(1,5) = 0.2e1 * alpha * mz / dspq * (mx * pz - mz * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (mx * pz - mz * px) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * mz + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * px;
//      Hessian(2,0) = 0.2e1 * alpha * (-nx * py + ny * px) / dspq * (-ny * pz + nz * py) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-ny * pz + nz * py) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (-nx * py + ny * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * px * pz;
//      Hessian(2,1) = 0.2e1 * alpha * (-nx * py + ny * px) / dspq * (nx * pz - nz * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (nx * pz - nz * px) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (-nx * py + ny * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * py * pz;
      Hessian(2,2) = 0.2e1 * alpha * pow(-mx * py + my * px, 0.2e1) / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-mx * py + my * px) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px, 0.2e1) - alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * (0.2e1 * py * py + 0.2e1 * px * px);
      Hessian(2,3) = 0.2e1 * alpha * mx / dspq * (-mx * py + my * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-mx * py + my * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * mx + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * py;
      Hessian(2,4) = 0.2e1 * alpha * my / dspq * (-mx * py + my * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-mx * py + my * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * my + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq * px;
      Hessian(2,5) = 0.2e1 * alpha * mz / dspq * (-mx * py + my * px) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-mx * py + my * px) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * mz + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * pz - 0.2e1 * qz);
//      Hessian(3,0) = 0.2e1 * alpha * nx / dspq * (-ny * pz + nz * py) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-ny * pz + nz * py) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * nx + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * px - 0.2e1 * qx);
//      Hessian(3,1) = 0.2e1 * alpha * nx / dspq * (nx * pz - nz * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (nx * pz - nz * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * nx + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * pz;
//      Hessian(3,2) = 0.2e1 * alpha * nx / dspq * (-nx * py + ny * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-nx * py + ny * px) * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * nx + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * py;
      Hessian(3,3) = 0.2e1 * alpha * mx * mx / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * mx * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(0.2e1 * px - 0.2e1 * qx, 0.2e1) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq;
      Hessian(3,4) = 0.2e1 * alpha * my / dspq * mx - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * mx * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * px - 0.2e1 * qx) * my + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * py - 0.2e1 * qy);
      Hessian(3,5) = 0.2e1 * alpha * mz / dspq * mx - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * mx * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * px - 0.2e1 * qx) * mz + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * pz - 0.2e1 * qz);
//      Hessian(4,0) = 0.2e1 * alpha * ny / dspq * (-ny * pz + nz * py) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-ny * pz + nz * py) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * ny + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * pz;
//      Hessian(4,1) = 0.2e1 * alpha * ny / dspq * (nx * pz - nz * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (nx * pz - nz * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * ny + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * py - 0.2e1 * qy);
//      Hessian(4,2) = 0.2e1 * alpha * ny / dspq * (-nx * py + ny * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-nx * py + ny * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * ny + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * px;
//      Hessian(4,3) = 0.2e1 * alpha * ny / dspq * nx - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * nx * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * px - 0.2e1 * qx) * ny + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * py - 0.2e1 * qy);
      Hessian(4,4) = 0.2e1 * alpha * my * my / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * my * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(0.2e1 * py - 0.2e1 * qy, 0.2e1) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq;
      Hessian(4,5) = 0.2e1 * alpha * mz / dspq * my - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * my * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * (0.2e1 * py - 0.2e1 * qy) * mz + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * py - 0.2e1 * qy) * (0.2e1 * pz - 0.2e1 * qz);
//      Hessian(5,0) = 0.2e1 * alpha * nz / dspq * (-ny * pz + nz * py) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-ny * pz + nz * py) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * nz + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * py;
//      Hessian(5,1) = 0.2e1 * alpha * nz / dspq * (nx * pz - nz * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (nx * pz - nz * px) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * nz + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * dpq * px;
//      Hessian(5,2) = 0.2e1 * alpha * nz / dspq * (-nx * py + ny * px) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-nx * py + ny * px) * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * nz + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) * (0.2e1 * pz - 0.2e1 * qz);
//      Hessian(5,3) = 0.2e1 * alpha * nz / dspq * nx - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * nx * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * px - 0.2e1 * qx) * nz + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * pz - 0.2e1 * qz);
//      Hessian(5,4) = 0.2e1 * alpha * nz / dspq * ny - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * ny * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * alpha * (nx * (px - qx) + ny * (py - qy) + nz * (pz - qz)) * dpq * (0.2e1 * py - 0.2e1 * qy) * nz + 0.2e1 * alpha * pow(nx * (px - qx) + ny * (py - qy) + nz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * (0.2e1 * py - 0.2e1 * qy) * (0.2e1 * pz - 0.2e1 * qz);
      Hessian(5,5) = 0.2e1 * alpha * mz * mz / dspq - 0.4e1 * alpha * (mx * (px - qx) + my * (py - qy) + mz * (pz - qz)) * dpq * mz * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * pow(dspq, -0.3e1) * pow(0.2e1 * pz - 0.2e1 * qz, 0.2e1) - 0.2e1 * alpha * pow(mx * (px - qx) + my * (py - qy) + mz * (pz - qz), 0.2e1) * dpq;
    }
    virtual void approximateCovar(const Surface &s1, const Surface &s2, DSMatrix &Hessian, DSMatrix &TransMeasCovar) const {
      const DVector &p(s1.p);
      //const DVector &n(s1.n);
      const DVector &q(s2.p);
      const DVector &m(s2.n);
      double px = p(0);
      double py = p(1);
      double pz = p(2);
      double qx = q(0);
      double qy = q(1);
      double qz = q(2);
      double mx = m(0);
      double my = m(1);
      double mz = m(2);
      double sigalph = std::max(s2.nStdDevRAD,0.01); // minimum 0,5°

      DVector unused_gradient;
      approximateAt(s1, s2, Hessian, unused_gradient);

      DMatrix Ezx(6,9); // corresponds to d^2 Err / dz dx with x as below and z as above
      Ezx(0,0) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py);
      Ezx(1,0) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qz;
      Ezx(2,0) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qy;
      Ezx(3,0) = 0.2e1 * mx * mx / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph - 0.4e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * pow(0.2e1 * px - 0.2e1 * qx, 0.2e1) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      Ezx(4,0) = 0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * my + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(5,0) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,1) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qz;
      Ezx(1,1) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px);
      Ezx(2,1) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qx;
      Ezx(3,1) = 0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * my + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(4,1) = 0.2e1 * my * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph - 0.4e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * pow(0.2e1 * py - 0.2e1 * qy, 0.2e1) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      Ezx(5,1) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,2) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qy;
      Ezx(1,2) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * qx;
      Ezx(2,2) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * pz - 0.2e1 * qz) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px);
      Ezx(3,2) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(4,2) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(5,2) = 0.2e1 * mz * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph - 0.4e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * pow(0.2e1 * pz - 0.2e1 * qz, 0.2e1) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      Ezx(0,3) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py);
      Ezx(1,3) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * pz;
      Ezx(2,3) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * py;
      Ezx(3,3) = 0.2e1 * mx / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph;
      Ezx(4,3) = 0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(5,3) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (px - qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (px - qx) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,4) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * pz;
      Ezx(1,4) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px);
      Ezx(2,4) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * px;
      Ezx(3,4) = 0.2e1 * mx / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (0.2e1 * px - 0.2e1 * qx);
      Ezx(4,4) = 0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph;
      Ezx(5,4) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (py - qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (py - qy) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,5) = 0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * py;
      Ezx(1,5) = 0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * px;
      Ezx(2,5) = 0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px);
      Ezx(3,5) = 0.2e1 * mx / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (0.2e1 * px - 0.2e1 * qx);
      Ezx(4,5) = 0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(5,5) = 0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * (pz - qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (pz - qz) * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph;
      Ezx(0,6) = -0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py);
      Ezx(1,6) = -0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * pz;
      Ezx(2,6) = -0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * py;
      Ezx(3,6) = -0.2e1 * mx * mx / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * px - 0.2e1 * qx) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * mx + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (0.2e1 * px - 0.2e1 * qx) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      Ezx(4,6) = -0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * my + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(5,6) = -0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mx * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * px + 0.2e1 * qx) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,7) = -0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * pz;
      Ezx(1,7) = -0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px);
      Ezx(2,7) = -0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * px;
      Ezx(3,7) = -0.2e1 * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * my - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * mx + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (0.2e1 * px - 0.2e1 * qx);
      Ezx(4,7) = -0.2e1 * my * my / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * py - 0.2e1 * qy) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * my + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (0.2e1 * py - 0.2e1 * qy) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      Ezx(5,7) = -0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * my * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * py + 0.2e1 * qy) * (0.2e1 * pz - 0.2e1 * qz);
      Ezx(0,8) = -0.2e1 * (-pz * my + py * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (-pz * my + py * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (-0.2e1 * (py - qy) * pz + 0.2e1 * (pz - qz) * py) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * py;
      Ezx(1,8) = -0.2e1 * (pz * mx - px * mz) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (pz * mx - px * mz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (0.2e1 * (px - qx) * pz - 0.2e1 * (pz - qz) * px) - 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * px;
      Ezx(2,8) = -0.2e1 * (-py * mx + px * my) / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mz + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (-py * mx + px * my) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (-0.2e1 * (px - qx) * py + 0.2e1 * (py - qy) * px);
      Ezx(3,8) = -0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * mx + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * px - 0.2e1 * qx) * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * mx + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (0.2e1 * px - 0.2e1 * qx);
      Ezx(4,8) = -0.2e1 * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph * my + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (0.2e1 * py - 0.2e1 * qy) * mz - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * my + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (0.2e1 * py - 0.2e1 * qy);
      Ezx(5,8) = -0.2e1 * mz * mz / (pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1)) / sigalph + 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * mz * (0.2e1 * pz - 0.2e1 * qz) - 0.2e1 * ((px - qx) * mx + (py - qy) * my + (pz - qz) * mz) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * mz + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.3e1) / sigalph * (-0.2e1 * pz + 0.2e1 * qz) * (0.2e1 * pz - 0.2e1 * qz) + 0.2e1 * pow((px - qx) * mx + (py - qy) * my + (pz - qz) * mz, 0.2e1) * pow(pow(px - qx, 0.2e1) + pow(py - qy, 0.2e1) + pow(pz - qz, 0.2e1), -0.2e1) / sigalph;
      DMatrix EzxT = ublas::trans(Ezx);

      DSMatrix CovZ = ublas::zero_matrix<double>(9,9); // covariance of the measurements "z", here: [px, py, pz, mx, my, mz, qx, qy, qz]
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(0,3), ublas::range(0,3)) = s1.pCovar3D;
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(3,6), ublas::range(3,6)) = s2.nCovar3D;
      ublas::matrix_range< DSMatrix >(CovZ, ublas::range(6,9), ublas::range(6,9)) = s2.pCovar3D;

//      std::cout << "§" << Ezx.size1() << "_" << Ezx.size2() << std::flush;
//      std::cout << ";" << CovZ.size1() << "_" << CovZ.size2() << std::flush;
      DMatrix tmp = ublas::prod(Ezx,CovZ);
//      std::cout << "#" << tmp.size1() << "_" << tmp.size2() << std::flush;
//      std::cout << "%" << EzxT.size1() << "_" << EzxT.size2() << std::flush;
      TransMeasCovar = ublas::prod(tmp,EzxT);
    }
  };

} // namespace ICP

#include "ICPLinearized.tcc"

#endif /*ICP_LINEARIZED_H_*/
