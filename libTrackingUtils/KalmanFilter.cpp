#include "KalmanFilter.hpp"

#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>
#include <climits>
#include <set>
#include "MatrixDefs.hpp"

using namespace std;
using namespace matrixTools;

///////////////////////////////////////////////////////
//////////////  Impl. KalmanFilter   ///////////////
///////////////////////////////////////////////////////
KalmanFilter::KalmanFilter(const DVector &initState, const DSMatrix &initVariance) :
   x(initState)
  ,P(initVariance)
  ,fullQ(ublas::zero_matrix<double>(1,1))
  ,fullF(ublas::zero_matrix<double>(1,1))
  ,fullFt(ublas::zero_matrix<double>(1,1))
  ,fullH(ublas::zero_matrix<double>(1,1))
  ,fullHt(ublas::zero_matrix<double>(1,1))
  ,Q(&fullQ)
  ,F(&fullF)
  ,Ft(&fullFt)
  ,H(&fullH)
  ,Ht(&fullHt)
{
}

KalmanFilter::KalmanFilter(const DVector &initState, const DSMatrix &initVariance, const mdefs::DSMatrix &predictVar, const mdefs::DMatrix &sysMat, const mdefs::DMatrix &measMat) :
   x(initState)
  ,P(initVariance)
  ,Q(&fullQ)
  ,F(&fullF)
  ,Ft(&fullFt)
  ,H(&fullH)
  ,Ht(&fullHt)
{
  setPredictionVariance(predictVar);
  setSystemModel(sysMat);
  setMeasurementModel(measMat);
}

KalmanFilter::KalmanFilter(const mdefs::DVector &initState, const mdefs::DSMatrix &initVariance, const mdefs::DSMatrix &predictVar, const mdefs::DMatrix &sysMat, const mdefs::DMatrix &sysMatT, const mdefs::DMatrix &measMat, const mdefs::DMatrix &measMatT) :
   x(initState)
  ,P(initVariance)
  ,Q(&fullQ)
  ,F(&fullF)
  ,Ft(&fullFt)
  ,H(&fullH)
  ,Ht(&fullHt)
{
  setPredictionVarianceRef(predictVar);
  setSystemModelRef(sysMat, sysMatT);
  setMeasurementModelRef(measMat, measMatT);
}

KalmanFilter::KalmanFilter(KalmanFilter &other) :
   x(other.x)
  ,P(other.P)
  ,Q(other.Q)
  ,F(other.F)
  ,Ft(other.Ft)
  ,H(other.H)
  ,Ht(other.Ht)
{
  if (other.Q == &other.fullQ) { // if pointer refers to private matrix, copy it
    fullQ = other.fullQ;
    Q = &fullQ; // re-set pointer
  }
  if (other.F == &other.fullF) { // if pointer refers to private matrix, copy it
    fullF = other.fullF;
    F = &fullF; // re-set pointer
  }
  if (other.Ft == &other.fullFt) { // if pointer refers to private matrix, copy it
    fullFt = other.fullFt;
    Ft = &fullFt; // re-set pointer
  }
  if (other.H == &other.fullH) { // if pointer refers to private matrix, copy it
    fullH = other.fullH;
    H = &fullH; // re-set pointer
  }
  if (other.Ht == &other.fullHt) { // if pointer refers to private matrix, copy it
    fullHt = other.fullHt;
    Ht = &fullHt; // re-set pointer
  }
}

KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::setPredictionVariance(const mdefs::DSMatrix &matrix)
{
  fullQ = matrix;
  Q = &fullQ;
}

void KalmanFilter::setSystemModel(const mdefs::DMatrix &matrix)
{
  fullF = matrix;
  fullFt = ublas::trans(fullF);
  F = &fullF;
  Ft = &fullFt;
}

void KalmanFilter::setMeasurementModel(const mdefs::DMatrix &matrix)
{
  fullH = matrix;
  fullHt = ublas::trans(fullH);
  H = &fullH;
  Ht = &fullHt;
}

void KalmanFilter::setPredictionVarianceRef(const mdefs::DSMatrix &matrix)
{
  Q = &matrix;
}

void KalmanFilter::setSystemModelRef(const mdefs::DMatrix &matrix, const mdefs::DMatrix &matrixTransp)
{
  F = &matrix;
  Ft = &matrixTransp;
}

void KalmanFilter::setMeasurementModelRef(const mdefs::DMatrix &matrix, const mdefs::DMatrix &matrixTransp)
{
  H = &matrix;
  Ht = &matrixTransp;
}

void KalmanFilter::predict()
{
  x = prod(*F,x); // state prediction: x=F*x
  DMatrix Qnew = prod(P,*Ft);
  Qnew = prod(*F,Qnew);
  Qnew = Qnew + *Q;
  ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SAQnew(Qnew);
  P = SAQnew; // covariance prediction P=F*P*F^t+Q
  for (unsigned int i=0; i<12; ++i) {
    assert(!std::isnan(x(i)) && "KalmanFilter::predict: state x is NaN");
    for (unsigned int j=0; j<12; ++j) {
      assert(!std::isnan(P(i,j)) && "KalmanFilter::predict: covariance P is NaN");
    }
  }
}

void KalmanFilter::update(DVector &z, DMatrix &R)
{
  namespace lapack = boost::numeric::bindings::lapack;

  // check input variables:
  assert((R.size1() == z.size()) && "KalmanFilter::update: measurement covariance has wrong dimensionality");
  assert((R.size2() == z.size()) && "KalmanFilter::update: measurement covariance has wrong dimensionality");
  for (unsigned int i=0; i<z.size(); ++i) {
    assert(!std::isnan(z(i)) && "KalmanFilter::update: measurement is NaN");
    for (unsigned int j=0; j<z.size(); ++j) {
      assert(!std::isnan(R(i,j)) && "KalmanFilter::update: covariance is NaN");
    }
  }

  // calculate expected measurement
  DVector zExp = prod(*H,x);
  
  // calculate the uncertainty of the state transformed to the measurement space (S=H*P*Ht+R / Si=inverse(S))
  DMatrix S = prod(P,*Ht);
  S = prod(*H,S);
  S = S + R;
  for (unsigned int i=0; i<S.size1(); ++i) {
    for (unsigned int j=0; j<S.size2(); ++j) {
      assert(!std::isnan(S(i,j)) && "KalmanFilter::update: S is NaN");
    }
  }
  // now, S should be a symmetric, positive semi-definite matrix (i.e. a valid covariance matrix)
  assert(S.size1() == S.size2());

  DCMatrix Si = invSym(S); // calculate inverse of S
//  cout << "  Si'=" << Si << flush;
  for (unsigned int i=0; i<Si.size1(); ++i) {
    for (unsigned int j=0; j<Si.size2(); ++j) {
      assert(!std::isnan(Si(i,j)) && "KalmanFilter::update: Si is NaN");
    }
  }
  
  // calculate Kalman Gain W=P*Ht*Si
  DMatrix W = prod(*Ht,Si);
//  cout << "  W=" << W << flush;
  W = prod(P,W);
//  cout << "  W'=" << W << flush;
  for (unsigned int i=0; i<W.size1(); ++i) {
    for (unsigned int j=0; j<W.size2(); ++j) {
      assert(!std::isnan(W(i,j)) && "KalmanFilter::update: W is NaN");
    }
  }
  
  // update state and covariance x=x+W(z-zExp), P=P-W*S*Wt
  x = x + prod(W,z-zExp);
  DMatrix Pnew = prod(S,trans(W));
  for (unsigned int i=0; i<Pnew.size1(); ++i) {
    for (unsigned int j=0; j<Pnew.size2(); ++j) {
      assert(!std::isnan(Pnew(i,j)) && "PCTrack::update: Pnew is NaN");
    }
  }
  Pnew = prod(W,Pnew);
//  cout << "  Pnew'=" << Pnew << flush;
  for (unsigned int i=0; i<Pnew.size1(); ++i) {
    for (unsigned int j=0; j<Pnew.size2(); ++j) {
      assert(!std::isnan(Pnew(i,j)) && "PCTrack::update: Pnew' is NaN");
    }
  }
  ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SAPnew(Pnew);
  P = P-SAPnew;

//  cout << "VAR:" << P;
//  cout << "  P=" << P << flush;
  for (unsigned int i=0; i<12; ++i) {
    assert(!std::isnan(x(i)) && "PCTrack::update: new state is NaN");
    for (unsigned int j=0; j<12; ++j) {
      assert(!std::isnan(P(i,j)) && "PCTrack::update: new covariance is NaN");
    }
  }
}
