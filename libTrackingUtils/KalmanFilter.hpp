#ifndef KALMANFILTER_H_
#define KALMANFILTER_H_

#include <stdexcept>

#include "MatrixDefs.hpp"

namespace mdefs = matrixTools;

/*!
 * \class KalmanFilter
 * \brief Implements a basic linear Kalman Filter for sequential state estimation
 *
 * This class allows several filters to share their prediction-covariance-matrix, system matrix, and measurement matrix.
 * This reduces memory requirements if many tracks exist.
 */
class KalmanFilter
{
public:
  KalmanFilter(const mdefs::DVector &initState, const mdefs::DSMatrix &initVariance);
  KalmanFilter(const mdefs::DVector &initState, const mdefs::DSMatrix &initVariance, const mdefs::DSMatrix &predictVar, const mdefs::DMatrix &sysMat, const mdefs::DMatrix &measMat); //!< constructor which gives the filter its own copy of the specified matrices
  KalmanFilter(const mdefs::DVector &initState, const mdefs::DSMatrix &initVariance, const mdefs::DSMatrix &predictVar, const mdefs::DMatrix &sysMat, const mdefs::DMatrix &sysMatT, const mdefs::DMatrix &measMat, const mdefs::DMatrix &measMatT); //!< constructor that makes the filter reference the given external matrices
  KalmanFilter(KalmanFilter &other);
  virtual ~KalmanFilter();
  
  void setPredictionVariance(const mdefs::DSMatrix &matrix);
  void setSystemModel(const mdefs::DMatrix &matrix);
  void setMeasurementModel(const mdefs::DMatrix &matrix);
  void setPredictionVarianceRef(const mdefs::DSMatrix &matrix); //!< makes the KalmanFilter use the specified matrix as prediction variance. reference will be stored, so no data is copied => the specified matrix must be valid as long as the KalmanFilter is running
  void setSystemModelRef(const mdefs::DMatrix &matrix, const mdefs::DMatrix &matrixTransp); //!< makes the KalmanFilter use the specified matrix as system model. reference will be stored, so no data is copied => the specified matrix must be valid as long as the KalmanFilter is running
  void setMeasurementModelRef(const mdefs::DMatrix &matrix, const mdefs::DMatrix &matrixTransp); //!< makes the KalmanFilter use the specified matrix as measurement model. reference will be stored, so no data is copied => the specified matrix must be valid as long as the KalmanFilter is running

  void predict();
  void update(mdefs::DVector &z, mdefs::DMatrix &R); //!< update: measurement vector and covariance matrix
  void getState(mdefs::DVector &state) const {state = x;}; //!< returns state vector
  void getCovar(mdefs::DMatrix &covar) const {covar = P;}; //!< returns covariance matrix of current state
  
private:
  KalmanFilter() :Q(&fullQ),F(&fullF),Ft(&fullFt),H(&fullH),Ht(&fullHt){}; //!< make private so it cannot be used

  // Kalman Filter variables:
  mdefs::DVector            x; //!< current state (n x 1)
  mdefs::DSMatrix           P; //!< covariance of the state (n x n)
  // The following allows to either assign each KalmanFilter object its own matrix or to use one common matrix for many Filters, so memory requirements are minimal
  mdefs::DSMatrix            fullQ; //!< prediction-covariance of the state (size n x n), might not be used if reference below refers to external matrix
  mdefs::DMatrix            fullF; //!< system matrix for predicting state (size n x n), might not be used if reference below refers to external matrix
  mdefs::DMatrix            fullFt; //!< transposed system matrix, might not be used if reference below refers to external matrix
  mdefs::DMatrix            fullH; //!< measurement matrix (size n x m), m = size of measurement vector, might not be used if reference below refers to external matrix
  mdefs::DMatrix            fullHt; //!< transposed measurement matrix, might not be used if reference below refers to external matrix
  const mdefs::DSMatrix      *Q; //!< pointer to the prediction-covariance of the state (size n x n)
  const mdefs::DMatrix      *F; //!< pointer to the system matrix for predicting state (size n x n)
  const mdefs::DMatrix      *Ft; //!< pointer to the transposed system matrix
  const mdefs::DMatrix      *H; //!< pointer to the measurement matrix (size n x m), m = size of measurement vector
  const mdefs::DMatrix      *Ht; //!< pointer to the transposed measurement matrix
};

#endif /*KALMANFILTER_H_*/
