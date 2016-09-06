#include "PointCloudTrack.hpp"

#define BOOST_FILESYSTEM_VERSION 2

#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <climits>
#include <set>
#include "LidarFrame.hpp"
#include "MatrixDefs.hpp"
#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;
namespace fs = boost::filesystem;

///////////////////////////////////////////////////////
////////////////   Static Functions   /////////////////
///////////////////////////////////////////////////////
PointCloudTrack::IdT PointCloudTrack::getNewUid()
{
  static IdT latestUid = 0; // will only be initialized on 1st call!!
  latestUid++;
  return latestUid;
}

DMatrix PointCloudTrack::F = ublas::identity_matrix<double>(12,12);
DMatrix PointCloudTrack::Ft = ublas::identity_matrix<double>(12,12);
DMatrix PointCloudTrack::H = ublas::zero_matrix<double>(6,12);
DMatrix PointCloudTrack::Ht = ublas::zero_matrix<double>(12,6);
DSMatrix PointCloudTrack::Q = ublas::zero_matrix<double>(12,12);

void PointCloudTrack::initializeFH(const DSMatrix &predictionVariance)
{
  static bool FHinitialized = false; // will only be initialized on 1st call!!
  if (!FHinitialized) {
    FHinitialized = true;
    ParameterHeap* params = ParameterHeap::get();
    for (int i=0; i<6; ++i) {F(i,i+6)=params->defaultDeltaT;} // assume constant
//    for (int i=0; i<6; ++i) {H(i,i+6)=1;} // velocity as measurement
    for (int i=0; i<6; ++i) {H(i,i)=1;} // position as measurement
    Ft = trans(F);
    Ht = trans(H);
    Q = predictionVariance;
    cout << endl << " Track::SystemMatrix=" << F;
    cout << endl << " Track::MeasMatrix=" << H;
    cout << endl << " Track::PredictCovar=" << Q;
  }
}

///////////////////////////////////////////////////////
//////////////  Impl. PointCloudTrack   ///////////////
///////////////////////////////////////////////////////
PointCloudTrack::PointCloudTrack(const DVector &initState, const DSMatrix &initVariance, const DSMatrix &predictionVariance, float colC_, float rowC_) :
   mergeDecision(MDKeep)
  ,uid(getNewUid())
  ,colorR(0.2 + 0.8*(double)rand()/(double)RAND_MAX)
  ,colorG(0.2 + 0.8*(double)rand()/(double)RAND_MAX)
  ,colorB(0.2 + 0.8*(double)rand()/(double)RAND_MAX)
//  ,idx(UINT_MAX)
  ,age(0)
  ,lastUpdateCounter(0)
  ,isMovingObj(false)
  ,colC(colC_)
  ,rowC(rowC_)
  ,x(initState)
  ,P(initVariance)
  ,expectingPrediction(true)
  ,pointCloudPGrid(ParameterHeap::get()->trkPGridWidth) // grid-width in meter
  ,pointCloudNGrid(ParameterHeap::get()->trkNGridWidth)
  ,pointsCleared(false)
  ,wasCloned(false)
//  ,projPixCount(0)
  ,histXPredict(ublas::zero_vector<double>(1))
  ,histPPredict(ublas::zero_matrix<double>(1,1))
  ,histZExp(ublas::zero_vector<double>(1))
  ,histZ(ublas::zero_vector<double>(1))
  ,histR(ublas::zero_matrix<double>(1,1))
  ,histKGain(ublas::zero_matrix<double>(1,1))
  ,histRegUsedFeatureMatching(false)
{
  if (ParameterHeap::get()->outDir == "") {
    worldmapFilename =  "";
  } else {
    worldmapFilename =  ParameterHeap::get()->outDir+"/world.p3d";
    fs::path wmfp(worldmapFilename);
    if ((uid == 1) && (fs::exists(wmfp))) // remove existing world track
      fs::remove(wmfp);
  }
  assert((x.size() == 12) && "PointCloudTrack: Init state is not 12-dimensional");
  assert((P.size1() == 12) && "PointCloudTrack: Init covar is not 12x12-dimensional");
  assert((P.size2() == 12) && "PointCloudTrack: Init covar is not 12x12-dimensional");
  initializeFH(predictionVariance);
  StateHistory hist;
  hist.x=x;
  hist.P=P;
  hist.wasMeasured=true;
  getRt2EgoCS(hist.R, hist.t);
  histAllStates.push_back(hist);
}

PointCloudTrack::PointCloudTrack(PointCloudTrack &other) :
   mergeDecision(other.mergeDecision)
  ,uid(other.uid)
  ,colorR(other.colorR)
  ,colorG(other.colorG)
  ,colorB(other.colorB)
//  ,idx(other.idx)
  ,age(other.age)
  ,lastUpdateCounter(other.lastUpdateCounter)
  ,isMovingObj(other.isMovingObj)
  ,colC(other.colC)
  ,rowC(other.rowC)
  ,x(other.x)
  ,P(other.P)
  ,expectingPrediction(other.expectingPrediction)
  ,pointCloudPGrid(other.pointCloudPGrid)
  ,pointCloudNGrid(other.pointCloudNGrid)
  ,pointsCleared(other.pointsCleared)
  ,worldmapFilename(other.worldmapFilename)
  ,wasCloned(false)
//  ,projPixCount(other.projPixCount)
//  ,segments(other.segments)
//  ,segLastFrame(other.segLastFrame)
  ,histAllStates(other.histAllStates)
  ,histXPredict(other.histXPredict)
  ,histPPredict(other.histPPredict)
  ,histZExp(other.histZExp)
  ,histZ(other.histZ)
  ,histR(other.histR)
  ,histKGain(other.histKGain)
  ,histRegICPError(other.histRegICPError)
  ,histRegUsedFeatureMatching(other.histRegUsedFeatureMatching)  
{
  other.wasCloned = true;
  // as 1 object already exists ("other") that constructor has already initialized the static members F/H
}

PointCloudTrack::~PointCloudTrack()
{
  if ((uid == 1) && (wasCloned = false) && (worldmapFilename != "")) { // destroy world track and no copy exists
    ofstream mapFile(worldmapFilename.c_str(), ios_base::app); // append point cloud to file
    BOOST_AUTO(piEnd,pointCloudNGrid.end());
    for (BOOST_AUTO(pi,pointCloudNGrid.begin()); pi!=piEnd; ++pi) {
      const mdefs::DVector &pos = pi->s->position;
      mapFile << pos[0] << " " << pos[1] << " " << pos[2] << " ";
    }
  }
}

void PointCloudTrack::merge(PointCloudTrack *other, bool clearCurrentAppearance, bool addAppearance, bool checkNeighbors, unsigned int historyIdx)
{
  // check prerequisite: both merging tracks should be in the same pred/update-state
  if (expectingPrediction != other->expectingPrediction) {
    cerr << endl << "merging tracks with different expectedPredition-flag" << flush;
  }
  //age = max(age, other->age); // this should indicate the history-length
  lastUpdateCounter = min(lastUpdateCounter, other->lastUpdateCounter);
  isMovingObj = isMovingObj || other->isMovingObj;
  // keep uid, colorR/G/B and idx, as they represent this track object
  // keep position (was arbitrarily initialized, thus we can select this one)
//  x(0..5)
  // mix only velocities:
//  double mix1Fac = (double)(getPointCount()*age)/(double)(getPointCount()*age + other->getPointCount()*other->age);
//  double mix2Fac = 1.0 - mix1Fac;
//  x(6) = mix1Fac*x(6) + mix2Fac*other->x(6);
//  x(7) = mix1Fac*x(7) + mix2Fac*other->x(7);
//  x(8) = mix1Fac*x(8) + mix2Fac*other->x(8);
//  x(9) = mix1Fac*x(9) + mix2Fac*other->x(9);
//  x(10) = mix1Fac*x(10) + mix2Fac*other->x(10);
//  x(11) = mix1Fac*x(11) + mix2Fac*other->x(11);
//  P = mix1Fac*P + mix2Fac*other->P;
  // to map to current Track-CS, 1st execute R2/t2, then R1/t1 --> R'=R1*R2, t'=R1*t2+t1
  if (addAppearance && clearCurrentAppearance && !pointsCleared) { // only clear on first merge so that multiple tracks can be merged into this one
    pointCloudPGrid.clear();
    pointCloudNGrid.clear();
    pointsCleared = true;
  }
  if (addAppearance) {
    DMatrix R1, R1T, pC, nC;
    DVector t1;
    other->getRt2OtherTrkCS(R1, t1, *this, historyIdx, historyIdx);
    R1T = trans(R1);
    PointIterator os, osend; other->getPoints(os, osend);
    //cout << endl << " [" << getUID() << "]" << flush;
    while (os != osend){
      Surface::SPtr s = os->s;
      pC = prod(s->ptCovar,R1T);
      pC = prod(R1,pC);
      nC = prod(s->nrCovar,R1T);
      nC = prod(R1,nC);
      addPointTrkCS(prod(R1,s->position) + t1, prod(R1,s->normal), s->normalConfidence, s->normalStdDevRAD, checkNeighbors, pC, nC);
      ++os;
    }
  }
  other->pointCloudPGrid.clear();
  other->pointCloudNGrid.clear();
  // keep state history and all other history-variables
}

void PointCloudTrack::predict()
{
  if (expectingPrediction == false) // if no measurements are available at least an empty update() should have been called inbetween two prediction steps
    cerr << "PointCloudTrack::predict(): expected update, but prediction was called" << endl;
  expectingPrediction = false;
//  cout << " PREDICT:" << x;
  x = prod(F,x); // state prediction: x=F*x
//  cout << "->" << x;
//  Matrix Qnew = prod(F,prod(P,Ft))+Q;
  DMatrix Qnew = prod(P,Ft);
  Qnew = prod(F,Qnew);
  Qnew = Qnew + Q;
  ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SAQnew(Qnew);
  P = SAQnew; // covariance prediction P=F*P*F^t+Q
//  cout << "VAR:(" << P(0,0) << "," << P(1,1) << "," << P(2,2) << "," << P(3,3) << "," << P(4,4) << "," << P(5,5)
//       << "," << P(6,6) << "," << P(7,7) << "," << P(8,8) << "," << P(9,9) << "," << P(10,10) << "," << P(11,11) << ")";
//  segments.swap(segLastFrame); // move the currently associated segments into segLastFrame
//  segments.clear();
//  projPixCount = 0;
  // save prediction state in history variables:
  histXPredict = x;
  histPPredict = P;
  pointsCleared = false;

  for (unsigned int i=0; i<12; ++i) {
    assert(!std::isnan(x(i)) && "PCTrack::predict: state x is NaN");
    for (unsigned int j=0; j<12; ++j) {
      assert(!std::isnan(P(i,j)) && "PCTrack::predict: covariance P is NaN");
    }
  }
}

void PointCloudTrack::update(DVector &z, DMatrix &R)
{
  namespace lapack = boost::numeric::bindings::lapack;

  if (expectingPrediction == true) { // update was already called at least once
    cerr << "PointCloudTrack::predict(): expected prediction, but update was called" << endl;
    histAllStates.pop_back(); // remove state added at last update-call
    --age;
  }
  lastUpdateCounter = 0;

  // check input variables:
//  cout << endl << " Update: predicted x=" << x;
//  cout << "P=" << P;
//  cout << endl << ", z=" << z;
//  cout << "R=" << R << flush;
  assert((z.size() == 6) && "PCTrack::update: measurement has wrong dimensionality");
  assert((R.size1() == 6) && "PCTrack::update: covariance has wrong dimensionality");
  assert((R.size2() == 6) && "PCTrack::update: covariance has wrong dimensionality");
  for (unsigned int i=0; i<6; ++i) {
    assert(!std::isnan(z(i)) && "PCTrack::update: measurement is NaN");
    for (unsigned int j=0; j<6; ++j) {
      assert(!std::isnan(R(i,j)) && "PCTrack::update: covariance is NaN");
    }
  }

  // TODO(9): reset coordinate system (and thus whole track state) if it moved too far away

  // calculate expected measurement
  DVector zExp = prod(H,x);
  histZExp = zExp;
  histZ = z;
  histR = R;
//  cout << endl << " ->zExp=" << zExp;
  
  // calculate the uncertainty of the state transformed to the measurement space (S=H*P*Ht+R / Si=inverse(S))
  DMatrix S = prod(P,Ht);
  S = prod(H,S);
  S = S + R;
//  cout << "  S=" << S << flush;
  for (unsigned int i=0; i<S.size1(); ++i) {
    for (unsigned int j=0; j<S.size2(); ++j) {
      assert(!std::isnan(S(i,j)) && "PCTrack::update: S is NaN");
    }
  }
  // now, S should be a symmetric, positive semi-definite matrix (i.e. a valid covariance matrix)
  assert(S.size1() == S.size2());

  DCMatrix Si = invSym(S); // calculate inverse of S

  //  cout << "  Si'=" << Si << flush;
  for (unsigned int i=0; i<Si.size1(); ++i) {
    for (unsigned int j=0; j<Si.size2(); ++j) {
      assert(!std::isnan(Si(i,j)) && "PCTrack::update: Si is NaN");
    }
  }
  
  // calculate Kalman Gain W=P*Ht*Si
  DMatrix W = prod(Ht,Si);
//  cout << "  W=" << W << flush;
  W = prod(P,W);
  histKGain = W;
//  cout << "  W'=" << W << flush;
  for (unsigned int i=0; i<W.size1(); ++i) {
    for (unsigned int j=0; j<W.size2(); ++j) {
      assert(!std::isnan(W(i,j)) && "PCTrack::update: W is NaN");
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
//  cout << "  new x=" << x << flush;
//  cout << "  new P=" << P << flush;

//  cout << "VAR:" << P;
  StateHistory hist;
  hist.x=x;
  hist.P=P;
  hist.wasMeasured=true;
  getRt2EgoCS(hist.R, hist.t);
  histAllStates.push_back(hist);
  ++age;

//  cout << "  P=" << P << flush;
  for (unsigned int i=0; i<12; ++i) {
    assert(!std::isnan(x(i)) && "PCTrack::update: new state is NaN");
    for (unsigned int j=0; j<12; ++j) {
      assert(!std::isnan(P(i,j)) && "PCTrack::update: new covariance is NaN");
    }
  }
  expectingPrediction = true;

//  cout << endl << "state-history:";
//  BOOST_FOREACH(StateHistory &hist, histAllStates)
//    cout << "  (" << hist.x << ")";
//  cout << flush;
}

void PointCloudTrack::update()
{
  if (expectingPrediction == true) { // update was already called at least once -> do nothing
    cerr << "PointCloudTrack::predict(): expected prediction, but (empty) update was called" << endl;
  } else { // update is called for the first time (after predict)
    ++age;
    ++lastUpdateCounter;
    histZExp = ublas::zero_vector<double>(1);
    histZ = ublas::zero_vector<double>(1);
    histR = ublas::zero_matrix<double>(1,1);
    histKGain = ublas::zero_matrix<double>(1,1);
    StateHistory hist;
    hist.x=x;
    hist.P=P;
    hist.wasMeasured=false;
    getRt2EgoCS(hist.R, hist.t);
    histAllStates.push_back(hist);
    // no measurement, so just keep prediction
  }
  expectingPrediction = true;
}

void PointCloudTrack::addPoint(const DVector &pEgoCS, const DVector &nEgoCS, const double normalConfidence, const double normalStdDevRAD, bool checkNeighbors, const mdefs::DMatrix &ptc, const mdefs::DMatrix &nrc)
{
  assert(!std::isnan(pEgoCS(0)) && "pEgoCS(0) is NaN");
  assert(!std::isnan(pEgoCS(1)) && "pEgoCS(1) is NaN");
  assert(!std::isnan(pEgoCS(2)) && "pEgoCS(2) is NaN");
  assert(!std::isnan(nEgoCS(0)) && "nEgoCS(0) is NaN");
  assert(!std::isnan(nEgoCS(1)) && "nEgoCS(1) is NaN");
  assert(!std::isnan(nEgoCS(2)) && "nEgoCS(2) is NaN");
  /* Transformation of a point with uncertainty:
   * naming: q:=pTrackCS, p:=pEgoCS
   * 1) mean:
   *   q = R * p + t, where R,t from YawPitchRollXYZ_2_Rt_inv
   *   q = R^T * (p - t), where R,t from YawPitchRollXYZ_2_Rt, which we use in the following as t=(tx,ty,tz) and R are independently obtained straight from the state vector!
   *   ^^^^^ this is a non-linear function of yaw,pitch,roll,tx,ty,tz (track-state) and px,py,pz
   * 2) covariance:
   *   J = dTransform / d(roll,pitch,yaw,tx,ty,tz,px,py,pz) is the Jacobian of the non-linear transformation
   *   Covar_pTrackCS = J * Covar_Merge * J^T, with Covar_Merge being a 9x9 matrix: (C_state6x6    C_state/point)
   *                                                                                (C_state/point C_point3x3)
   *   With the assumption of point-measurement/track-state being independent we get two separate Jacobians:
   *   J_state = dTransform / d(roll,pitch,yaw,tx,ty,tz) = (M3x3 -R) being composed of a complicated M-matrix and the negative rotation matrix
   *   J_point = dTransform / d(px,py,pz) = R
   *   --> Covar_pTrackCS = J_state * Covar_state * J_state^T + R * Covar_point * R^T
   *                      = M*C1*M^T - 2*R*C2*M^T + R*C3*R^T  + R*Cp*R^T  (E1), with the Covar_state matrix being decomposed into (C1 C2, C2 C3) and Cp being "Covar_point"
   *
   * 3) simplification of covariance calculus:
   *   To keep calculation and memory requirement minimal, we will only store 1 covariance value for each point, assuming that its distribution is sigma*I3
   *   We thus simplify the different terms of E1 as follows:
   *   The first term gives the influence of the rotation uncertainty, we search the maximum uncertainty Cm1 = max{C11,C22,C33} of C1 and calculate ||p-t||^2*(sin(sqrt(Cm1)))^2
   *   We ignore the second term, assuming that angle and positional uncertainty are non-correlated
   *   We search the maximum uncertainty Cm3 = max{C11,C22,C33} of C3, replacing C3 by Cm3*I3. This makes R*R^T redundant, leading to the value Cm only
   *   As the Covar_point is given by one single covariance value only, sigma, the R*R^T collapses, leaving only this sigma
   *   Thus our final uncertainy value is:
   *   ==> Covar_pTrackCS = ||p-t||^2*(sin(sqrt(Cm1)))^2 + Cm3 + Covar_pEgoCS
   */
  DMatrix R; DVector t; getRt2EgoCS(R, t);
  DMatrix Rt = trans(R);
  DVector pRel = pEgoCS - t;
  DVector pTrackCS = prod(Rt,pRel);
  DVector nTrackCS = prod(Rt,nEgoCS);
  DMatrix ptCTrackCS = prod(Rt,ptc); ptCTrackCS = prod(ptCTrackCS,R);
  DMatrix nrCTrackCS = prod(Rt,nrc); nrCTrackCS = prod(nrCTrackCS,R);
  addPointTrkCS(pTrackCS, nTrackCS, normalConfidence, normalStdDevRAD, checkNeighbors, ptCTrackCS, nrCTrackCS);
}

void PointCloudTrack::addPointTrkCS(const DVector &pTrkCS, const DVector &nTrkCS, const double normalConfidence, const double normalStdDevRAD, bool checkNeighbors, const mdefs::DMatrix &ptCTrackCS, const mdefs::DMatrix &nrCTrackCS)
{
  assert(!std::isnan(pTrkCS(0)) && "pTrackCS(0) is NaN");
  assert(!std::isnan(pTrkCS(1)) && "pTrackCS(1) is NaN");
  assert(!std::isnan(pTrkCS(2)) && "pTrackCS(2) is NaN");
  assert(!std::isnan(nTrkCS(0)) && "nTrackCS(0) is NaN");
  assert(!std::isnan(nTrkCS(1)) && "nTrackCS(1) is NaN");
  assert(!std::isnan(nTrkCS(2)) && "nTrackCS(2) is NaN");

  /*
  double Cm1 = max(P(0,0),max(P(1,1),P(2,2))); // maximum covariance of angle uncertainty
  double Cm3 = max(P(3,3),max(P(4,4),P(5,5))); // maximum covariance of position uncertainty
  double Cpe = stdDevPerMeter*stdDevPerMeter*norm_2_square(pEgoCS);
  double varTrackCS = norm_2_square(pRel)*pow(sin(sqrt(Cm1)),2) + Cm3 + Cpe;
//  cout << "\t " << pEgoCS << " " << Cpe << " " << Cm1 << " " << Cm3 << " " << varTrackCS << flush;
  assert(!std::isnan(varTrackCS) && "varTrackCS is NaN");
  */

  Surface::SPtr s(new Surface(pTrkCS, nTrkCS, /*varTrackCS,*/ normalConfidence, normalStdDevRAD, ptCTrackCS, nrCTrackCS));
  if (age == 0) s->isVIP = true; // points added at track-creation are Very Important Points

  // decide whether to add point and which existing point to delete
  //const unsigned int adeptNbNSearched(6);
  ParameterHeap *params = ParameterHeap::get();
  const unsigned int nbNMinFound(3);
  const double filterMinWDSum(1.0); // dist=0.5 => wd=exp(-sqrt(dist))=0.5, VelodyneSLAM: 1.5
  const double minNormalConfidence(0.6); // VelodyneSLAM: 0.9
  const double minNormalCoincidence(cos(15./180.*M_PI)); // =0.966, VelodyneSLAM: 10.
  const double wFacVIP(10); // VIP points are this time more important
  const unsigned int searchCellRad(params->trkPGridSearchRad);
  const double nMaxDist(sqrt(3.)*2*(double)searchCellRad * pointCloudPGrid.getCellExtent());
  // idea for adapting point:
  //   if its normal confidence is low, don't change anything
  //   if its normal confidence is high enough:
  //     search for the next neighbors within distance threshold
  //     #if neighbors are not in all directions, continue (don't adept on the edge of the scan)
  //     if its normal does not match the neighboring normals, set normal confidence low and continue
  //     adept this point that it best represents a local plane together with the neighbors
  //     therefore move point "p" along its normal vector "n" by "c" meters in order to minimize the squared distances to the neighbors' planes (q_i/m_i):
  //     E=sum{m_i^T*d_i}^2, with d_i=(p + c*n)-q_i --> E=sum{(p-q_i)^T m_i + c* n^T m_i}^2=: sum{b_i + c*a_i}^2
  //     thus, the least squares solution for c is (AtA)^-1*Atb = sum{a_i*b_i}/sum{a_i*a_i}

  if ((checkNeighbors) && (age > 0) && (normalConfidence > minNormalConfidence)) {
    // get neighbors from pointCloudPGrid
    list<SurfacePProxy> neighbors;
    back_insert_iterator< list<SurfacePProxy> > nInserter(neighbors);
    pointCloudPGrid.get_cubic_neighbors(nInserter, SurfacePProxy(s), searchCellRad);
    if (neighbors.size() < nbNMinFound) {
      neighbors.clear();
      nInserter = back_insert_iterator< list<SurfacePProxy> >(neighbors);
      pointCloudPGrid.get_cubic_neighbors(nInserter, SurfacePProxy(s), searchCellRad+1);
    }
    if (neighbors.size() >= nbNMinFound) {
      unsigned int nbNValid = 0;
      double wSum = 0.0;
      double wdSum = 0.0;
      double aSum = 0.0; // LS-solution for adaptation
      double bSum = 0.0; // LS-solution for adaptation
      double avgNCI = 0.0;
      double avgNC = 0.0;
      BOOST_FOREACH(SurfacePProxy nn, neighbors) { // looping over neighbors, sorted by distance
        DVector di = s->position - nn->position;
        double dist = norm_2(di);
        if ((dist > nMaxDist) || (dist == 0.0)) // avoid div 0 below
          continue;
        double a = ublas::inner_prod(s->normal, nn->normal); // =1 if normals are in same direction, -1 if opposite direction
        double b = ublas::inner_prod(di, nn->normal); // =dist if other point lies in normal direction
        double wd = exp(-sqrt(dist));
        double nci = max(0.0, a); // normal coincidence
        double w = wd * (0.5 + 0.5 * fabs(b / dist)) * nci*nci;
        //               ^^ full weight along normal   ^^ normal coincidence
        if (nn->isVIP)
          w *= wFacVIP;
        wdSum += wd;
        wSum += w;
        aSum += a * w * a;
        bSum += a * w * b;
        avgNC += wd * nn->normalConfidence;
        avgNCI += wd * nci;
        nbNValid++;
      }
      avgNC /= max(1e-8,wdSum);
      avgNCI /= max(1e-8,wdSum);
      double c = -bSum / max(1e-8,aSum); // LS-adaption-factor, avoid DIV0
//      cout << " a" << flush;
      if ((fabs(c) < 0.5) && (wSum >= 0.1) && (wdSum >= filterMinWDSum) && (nbNValid >= nbNMinFound)
          && (avgNC >= minNormalConfidence) && (avgNCI >= minNormalCoincidence)) {
        s->position += s->normal * c;
//        cout << "_by" << c << flush;
      } else {
//        cout << "_SKIP" << flush;
//        if (!(fabs(c) < 0.5)) cout << "_c" << c << flush;
//        if (!(wSum >= 0.1)) cout << "_ws" << wSum << flush;
//        if (!(wdSum >= filterMinWDSum)) cout << "_wds" << wdSum << flush;
//        if (!(nbNValid >= nbNMinFound)) cout << "_nn" << nbNValid << flush;
//        if (!(avgNC >= minNormalConfidence)) cout << "_anc" << avgNC << flush;
//        if (!(avgNCI >= minNormalCoincidence)) cout << "_anci" << acos(avgNCI)*180./M_PI << "°" << flush;
      }
    } // end: if |neighbors|>0
  } // end: if adept

  // if a point in the corresponding cell already exists, don't add new point
  SurfacePProxy sp(s);
  if (pointCloudPGrid.begin(sp) == pointCloudPGrid.end(sp)) { // iterators on points within same cell
    // add new point to appearance point cloud
    assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::addPointTrkCS_1: different point count");
    pointCloudPGrid.insert(SurfacePProxy(s));
    pointCloudNGrid.insert(SurfaceNProxy(s));
    assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::addPointTrkCS_2: different point count");
  }
}

void PointCloudTrack::removeFarPoints(double egoDist)
{
  const double maxDS = egoDist*egoDist;
  DMatrix R;
  DVector t,p;
  getRt2EgoCS(R,t);
  //TODO(9): increase performance by using grid-structure (delete cells that are too far)
  std::vector<Surface::SPtr> neighbors2Delete; // collect points to be deleted, so for-loop will not be confused
  neighbors2Delete.reserve(pointCloudNGrid.size());
  BOOST_AUTO(piEnd,pointCloudNGrid.end());
  for (BOOST_AUTO(pi,pointCloudNGrid.begin()); pi!=piEnd; ++pi) {
    p = ublas::prod(R,(pi->s.get())->position) + t;
    if (norm_2_square(p) > maxDS)
      neighbors2Delete.push_back(pi->s);
  }
  // really delete points
  assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::removeFarPoints_1: different point count");
  ofstream *mapFile = NULL;
  if ((uid == 1) && (worldmapFilename != "")) mapFile = new ofstream(worldmapFilename.c_str(), ios_base::app); // append
  BOOST_FOREACH(Surface::SPtr s, neighbors2Delete) {
    if (mapFile) (*mapFile) << s->position[0] << " " << s->position[1] << " " << s->position[2] << " ";
    pointCloudPGrid.remove(SurfacePProxy(s));
    pointCloudNGrid.remove(SurfaceNProxy(s));
  }
  delete mapFile;
  assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::removeFarPoints_2: different point count");
}

void PointCloudTrack::removeNeighborhood(DVector pointEgoCS, int cellRange)
{
  DMatrix R;
  DVector t,pTrackCS;
  getRt2TrkCS(R,t);
  pTrackCS = ublas::prod(R,pointEgoCS) + t;
  removeNeighborhoodTrkCS(pTrackCS, cellRange);
}

void PointCloudTrack::removeNeighborhoodTrkCS(DVector pTrackCS, int range)
{
  Surface::SPtr s(new Surface(pTrackCS, DZeroVector(3)));
  SurfacePProxy spp(s);
  const PointCloudPGrid::index_t idx = pointCloudPGrid.getIndex(spp);
  for (int i1 = -range; i1<=range; ++i1) {
    for (int i2 = -range; i2<=range; ++i2) {
      for (int i3 = -range; i3<=range; ++i3) {
        PointCloudPGrid::index_t idxN = idx;
        idxN[0] += i1;
        idxN[1] += i2;
        idxN[2] += i3;
        PointCloudPGrid::iterator i = pointCloudPGrid.begin(idxN);
        PointCloudPGrid::iterator iE = pointCloudPGrid.end(idxN);
        if (i==iE) continue; // no points in this cell
        PointCloudPGrid::grid_iterator gi = pointCloudPGrid.fromIterator(i); // must be BEFORE ++i in loop (will jump to next cell)
        // first remove all corresponding surfaces from normal-grid
        assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::removeNeighborhoodTrkCS_1: different point count");
        for (; i!=iE; ++i) {
          Surface::SPtr s = i->s;
          pointCloudNGrid.remove(SurfaceNProxy(s));
        }
        // then delete the cell in point-grid
        assert((gi != pointCloudPGrid.gridEnd()) && "PointCloudTrack::removeNeighborhoodTrkCS: gi points to end()");
        pointCloudPGrid.erase(gi);
        assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::removeNeighborhoodTrkCS_2: different point count");
      }
    }
  }
}

PointCloudTrack::PointIterator PointCloudTrack::erasePoint(PointCloudTrack::PointIterator it)
{
  if (it == pointCloudNGrid.end())
    return it;
  assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::erasePoint: different point count");
  Surface::SPtr s = it->s;
  PointCloudPGrid::iterator pit = pointCloudPGrid.find(SurfacePProxy(s));
  assert((pit != pointCloudPGrid.end()) && "PointCloudTrack::erasePoint: point not found in P-Grid");
  pointCloudPGrid.erase(pit); // should work since it!=end()
  it = pointCloudNGrid.erase(it);
  assert((pointCloudPGrid.size() == pointCloudNGrid.size()) && "PointCloudTrack::erasePoint_2: different point count");
  return it;
}

bool PointCloudTrack::getRt2OtherTrkCS(mdefs::DMatrix &R, mdefs::DVector &t, const PointCloudTrack &other, unsigned int historyIdxThis, unsigned int historyIdxOther) const
{
  assert( (this->getAge()+1 == this->histAllStates.size()) && "PointCloudTrack::getRt2OtherTrkCS: history-size doesn't match age");
  if (this->getAge() < historyIdxThis)
    throw runtime_error("PointCloudTrack::getRt2OtherTrkCS: age / history size of this too small");
  if (other.getAge() < historyIdxOther)
    throw runtime_error("PointCloudTrack::getRt2OtherTrkCS: age / history size of other too small");
  
  // historyIdxThis/Other is 0 for current state, but state is at end of vector
  const StateHistory& hi_t =       histAllStates[      histAllStates.size()-1-historyIdxThis];
  const StateHistory& hi_o = other.histAllStates[other.histAllStates.size()-1-historyIdxOther];
  bool allStatesValid = hi_t.wasMeasured && hi_o.wasMeasured;
  
  // 2 transformation: 1) this to EgoCS   2) EgoCS to other (=inverse R/t)
  DVector t_t, t_o;
  DMatrix R_t, R_o;
  t_t = hi_t.t;
  R_t = hi_t.R;
  R_o = ublas::trans(hi_o.R); // inverse rotation
  t_o = -ublas::prod(R_o,hi_o.t); // inverse translation using already inverted rotation
  R = ublas::prod(R_o, R_t);
  t = ublas::prod(R_o, t_t) + t_o;
//  cout << endl << "state T: " << hi_t.t << " state O: " << hi_o.t << " -> diff: " << t << flush;
  return allStatesValid;
}

void PointCloudTrack::getStatistics(unsigned int &pointCnt, double &minDist, double &maxDist) const
{
  pointCnt = getPointCount();
  minDist = DBL_MAX;
  maxDist = 0;
  DMatrix R;
  DVector t,v;
  getRt2EgoCS(R,t);
  PointConstIterator pi, piend;
  getPoints(pi, piend);
  while (pi != piend) {
    DVector v = pi->s->position; // realative to Track-CS
    v = prod(R,v) + t; // make relative to Ego-CS
    double d = norm_2(v);
    if (d > maxDist) maxDist = d;
    if (d < minDist) minDist = d;
    ++pi;
  }
}

double PointCloudTrack::getAvgNormalConfidence() const {
  unsigned int nbPts = 0;
  double nConf = 0;
  PointConstIterator pi, piend;
  getPoints(pi, piend);
  while (pi != piend) {
    nbPts++;
    nConf += pi->s->normalConfidence;
    ++pi;
  }
  return nConf / (double)nbPts;
}

double PointCloudTrack::pointNormalDistributionRatio()
{
  double maxCells = pow(2*M_PI/pointCloudNGrid.getCellExtent(), pointCloudNGrid.dimensionality());
  return (double)(pointCloudNGrid.getCellCount()) / maxCells;
}

