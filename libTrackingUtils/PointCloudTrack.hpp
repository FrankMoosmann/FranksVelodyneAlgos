#ifndef POINTCLOUDTRACK_H_
#define POINTCLOUDTRACK_H_

#include <stdexcept>
#include <list>
#include <set>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>

#include "HomogeneousTransformationMatrix.hpp"
#include "MatrixDefs.hpp"
#include "GridND.hpp"
class LidarFrame;

//namespace ublas = boost::numeric::ublas;
namespace mdefs = matrixTools;
namespace htm = HomogeneousTransformationMatrix;

class PointCloudTrack
{
  friend class LidarFrame;

public:
  typedef boost::shared_ptr< PointCloudTrack >      SPtr;
  typedef unsigned long long                        IdT;
  // TODO (5): hash surfaces in only one 6-D grid (point+normal) with adaptive size
  // this will probably make subsampling better
  struct  Surface {
    typedef boost::shared_ptr< Surface > SPtr;
    bool isVIP; // Very Important Point = point added at track creation
    mdefs::DVector position;
    mdefs::DVector normal;
//    double variance; // each point has uncertainty relative to current position. assumed that covariance matrix = covariance * I3 (3x3 identity matrix)
    double normalConfidence;
    double normalStdDevRAD;
    mdefs::DMatrix ptCovar;
    mdefs::DMatrix nrCovar;
    Surface (const mdefs::DVector &p, const mdefs::DVector &n) //!< dummy-constructor
    : isVIP(false), position (p), normal(n) {};
    Surface (const mdefs::DVector &p, const mdefs::DVector &n, /*double v,*/ double nc, double nsd, const mdefs::DMatrix &ptc, const mdefs::DMatrix &nrc)
    : isVIP(false), position (p), normal(n), /*variance(v),*/ normalConfidence(nc), normalStdDevRAD(nsd), ptCovar(ptc), nrCovar(nrc) {};
    inline double operator() (int i) {return position(i);};
  };
  struct SurfacePProxy { // needed for grid-like indexing of 3D point coordinates
    Surface::SPtr s;
    SurfacePProxy () {};
    SurfacePProxy (Surface::SPtr s_) : s(s_) {};
    inline double operator() (int i) const {return (s->position)(i);};
    inline bool operator== (SurfacePProxy const &other) const {return (s == other.s);};
    inline Surface& operator*() {return *s;};
    inline Surface* operator->() {return s.get();};
  };
  struct SurfaceNProxy { // needed for grid-like indexing of 3D point coordinates
    Surface::SPtr s;
    SurfaceNProxy () {};
    SurfaceNProxy (Surface::SPtr s_) : s(s_) {};
    inline double operator() (int i) const {
      mdefs::DVector &n = s->normal;
      switch (i) {
        case 0: return atan2(n(1),n(0)); break; // yaw-angle
        case 1: return atan2(n(2),sqrt(pow(n(0),2)+pow(n(1),2))); break; // pitch-angle
        default: throw std::out_of_range("");
      }
    };
    inline bool operator== (SurfaceNProxy const &other) const {return (s == other.s);};
    inline Surface& operator*() {return *s;};
    inline Surface* operator->() {return s.get();};
  };
  struct  StateHistory {
    mdefs::DVector x; // state-vector
    mdefs::DMatrix P; // covariance matrix
    mdefs::DMatrix R; // rotation matrix pTrack => pEgo
    mdefs::DVector t; // translation vector pTrack => pEgo
    bool wasMeasured;
  };
//  typedef GridND<3,SurfacePProxy>   PointCloudPGrid;
  typedef GridND<2,SurfaceNProxy>   PointCloudNGrid;
  typedef GridND<3,SurfacePProxy,SingleElementCell<SurfacePProxy> >   PointCloudPGrid;
  typedef PointCloudNGrid::iterator         PointIterator;
  typedef PointCloudNGrid::const_iterator   PointConstIterator;
  typedef std::vector< StateHistory >::const_iterator StateIterator;
  typedef std::vector< StateHistory >::const_reverse_iterator RStateIterator;

  PointCloudTrack(const mdefs::DVector &initState, const mdefs::DSMatrix &initVariance, const mdefs::DSMatrix &predictionVariance, float colC, float rowC);
  PointCloudTrack(PointCloudTrack &other);
  virtual ~PointCloudTrack();
  
  void merge(PointCloudTrack *other, bool clearCurrentAppearance = false, bool addAppearance = true, bool checkNeighbors = true, unsigned int historyIdx = 0); //!< move all contents of other into this

  // methods concerning general track information
  IdT getUID() const {return uid;};
  double getColorR() const {return colorR;};
  double getColorG() const {return colorG;};
  double getColorB() const {return colorB;};
  void setColor(double r, double g, double b) {colorR=r; colorG=g; colorB=b;};
  //unsigned int getIndex() const {return idx;};
  //void setIndex(unsigned int i) {idx = i;};
  unsigned int getAge() const {return age;};
  unsigned int getLastUpdateCounter() const {return lastUpdateCounter;};
  bool isMovingObject() const {return isMovingObj;};
  void isMovingObject(bool b) {isMovingObj = b;};

  // methods concerning the track-state (rotations are 1)rz=yaw 2)ry=pitch 3)rx=roll):
  void                         predict();
  void                         update(mdefs::DVector &z, mdefs::DMatrix &R); //!< update: measurement vector [rx,ry,rz,tx,ty,tz] and covariance matrix
  void                         update(); //!< empty update, i.e. no measurement was made
  const mdefs::DVector&        getState() const {return x;}; //!< returns state vector [rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot]
  void                         getState(mdefs::DVector &state) const {state = x;}; //!< returns state vector [rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot]
  const mdefs::DSMatrix&       getCovar() const {return P;}; //!< returns covariance matrix of current state
  void                         getCovar(mdefs::DMatrix &covar) const {covar = P;}; //!< returns covariance matrix of current state
  void getPosition(mdefs::DVector &p) const {p.resize(3); p(0)=x(3); p(1)=x(4); p(2)=x(5);}; // returns the current position [x,y,z] (m)
  void getVelocity(mdefs::DVector &v) const {v.resize(3); v(0)=x(9); v(1)=x(10); v(2)=x(11);}; // returns the current velocity vector [xDot,yDot,zDot] (m/s)
  void getRt2EgoCS(mdefs::DMatrix &R, mdefs::DVector &t) const {htm::YawPitchRollXYZ_2_Rt(x(2), x(1), x(0), x(3), x(4), x(5), R, t);}; //!< get Rot-Matrix and Transl-Vector, so a Point in Track-CS can be transformed to Ego-CS: p'=R*p+t
  void getRt2TrkCS(mdefs::DMatrix &R, mdefs::DVector &t) const {htm::YawPitchRollXYZ_2_Rt_inv(x(2), x(1), x(0), x(3), x(4), x(5), R, t);}; //!< get Rot-Matrix and Transl-Vector, so a Point in Ego-CS can be transformed to Track-CS: p'=R*p+t
  void getHTM2EgoCS(mdefs::DMatrix &HTM) const                 {HTM = htm::YawPitchRollXYZ_2_HTM(x(2), x(1), x(0), x(3), x(4), x(5));}; //!< get HTM, so a homogeneous Point in Track-CS can be transformed to Ego-CS: ph'=HTM*ph
  void getHTM2TrkCS(mdefs::DMatrix &HTM) const                 {HTM = htm::YawPitchRollXYZ_2_HTM_inv(x(2), x(1), x(0), x(3), x(4), x(5));}; //!< get HTM, so a homogeneous Point in Ego-CS can be transformed to Track-CS: ph'=HTM*p
  /*! \brief get Rot-Matrix and Transl-Vector, so a Point in this Track-CS can be transformed to other tracks' CS: p'=R*p+t
   *  historyIndex [0..histCount-1] thereby refers to the relative state in the past.
   *  \return Returns true if both states were measured correctly and not predicted. */
  bool getRt2OtherTrkCS(mdefs::DMatrix &R, mdefs::DVector &t, const PointCloudTrack &other, unsigned int historyIdxThis = 0, unsigned int historyIdxOther = 0) const;
  static void state2Rt_T2E(const mdefs::DVector &state, mdefs::DMatrix &R, mdefs::DVector &t)        {htm::YawPitchRollXYZ_2_Rt(state(2), state(1), state(0), state(3), state(4), state(5), R, t);}; //!< get Rot-Matrix and Transl-Vector, so a Point in Track-CS can be transformed to Ego-CS: p'=R*p+t
  static void state2Rt_E2T(const mdefs::DVector &state, mdefs::DMatrix &R, mdefs::DVector &t)        {htm::YawPitchRollXYZ_2_Rt_inv(state(2), state(1), state(0), state(3), state(4), state(5), R, t);}; //!< get Rot-Matrix and Transl-Vector, so a Point in Ego-CS can be transformed to Track-CS: p'=R*p+t
  static void state2HTM_T2E(const mdefs::DVector &state, mdefs::DMatrix &HTM)                        {HTM = htm::YawPitchRollXYZ_2_HTM(state(2), state(1), state(0), state(3), state(4), state(5));}; //!< get HTM, so a homogeneous Point in Track-CS can be transformed to Ego-CS: ph'=HTM*ph
  static void Rt2state_T2E(const mdefs::DMatrix &R, const mdefs::DVector &t, mdefs::DVector &state)  {state.resize(6); htm::Rt_2_YawPitchRollXYZ(R, t, state(2), state(1), state(0), state(3), state(4), state(5));}; //!< turn a Rot-Matrix and Transl-Vector into a state vector
  static void HTM2state_T2E(const mdefs::DMatrix &HTM, mdefs::DVector &state)                        {state.resize(6); htm::HTM_2_YawPitchRollXYZ(HTM, state(2), state(1), state(0), state(3), state(4), state(5));}; //!< turn a HTM into a state vector
  
  // methods concerning the appearance point cloud:
  void getCenter(float &colC_, float &rowC_) const {colC_ = colC; rowC_ = rowC;};
  void addPoint(const mdefs::DVector &pEgoCS, const mdefs::DVector &nEgoCS, const double normalConfidence, const double normalStdDevRAD, bool checkNeighbors, const mdefs::DMatrix &ptc, const mdefs::DMatrix &nrc); //!< add a surface specified relative to Ego-CS
  void addPointTrkCS(const mdefs::DVector &pTrkCS, const mdefs::DVector &nTrkCS, const double normalConfidence, const double normalStdDevRAD, bool checkNeighbors, const mdefs::DMatrix &ptCTrkCS, const mdefs::DMatrix &nrCTrkCS); //!< add a surface specified relative to Track-CS
  void removeFarPoints(double egoDist); //!< removes all stored points which have a distance to the ego vehicle higher than the specified threshold
  void removeNeighborhood(mdefs::DVector pointEgoCS, int cellRange = 0); //!< removes all points in a cubic cell-neighborhood around the given point
  void removeNeighborhoodTrkCS(mdefs::DVector pointTrkCS, int cellRange = 0); //!< removes all points in the same cell as the given point
  PointIterator erasePoint(PointIterator it); //!< returns iterator on next element. only valid if inner grid-containers' erase method has this property (ok for list and SingleElementCell)
  void clearPoints() {pointCloudPGrid.clear(); pointCloudNGrid.clear();}; //!< removes all stored points
  size_t getPointCount() const {return pointCloudPGrid.size();};
  void getPoints(PointIterator &begin, PointIterator &end) {begin = pointCloudNGrid.begin();end = pointCloudNGrid.end();}; //!< get points, specified relative to Track-CS!
  void getPoints(PointConstIterator &begin, PointConstIterator &end) const {begin = pointCloudNGrid.begin();end = pointCloudNGrid.end();}; //!< get points, specified relative to Track-CS!
  void getStatistics(unsigned int &pointCount, double &minDist, double &maxDist) const;
  double getAvgNormalConfidence() const;
  double pointNormalDistributionRatio(); //!< returns 1 if normals in all directions are included and lower values otherwise. with a scan from one viewpoint a maximum value of 0.5 can be achieved
  template <size_t N>
  boost::array<unsigned int, N> getMotion2NormalHistogram(const mdefs::DVector &state1, const mdefs::DVector &state2, double *moveMainNormalDir = NULL) const; //!< calculate how many points move along their normal vector [N-1] or perpendicular to it [0] when track moves from state1 to state2 (6D-Vector)
  template <typename InsertIterator>
  unsigned int samplePointsUniformly(const unsigned int desiredCount, InsertIterator container) const; //!< sample uniformly over point location distribution. returns effective count sampled (<=desiredCount)
  template <typename InsertIterator>
  unsigned int sampleNormalsUniformly(const unsigned int desiredCount, InsertIterator container) const; //!< sample uniformly over normal vector distribution. returns effective count sampled (<=desiredCount)

  // TODO (7): remove this variable as soon as track-management works
//  switch (track->mergeDecision) {
//  case PointCloudTrack::MDIgnore :
//    break;
//  case PointCloudTrack::MDKeep :
//    break;
//  case PointCloudTrack::MDMerge :
//    break;
//  }
  enum MergeDecision {MDIgnore, MDKeep, MDMerge};
  MergeDecision mergeDecision;
  PointCloudTrack::SPtr trackMaxMergeScore; //!< track this one should be merged into (buffer for decision)
  //tmMergeDecision
  //if (track->getMergeTrackIcpSuccess()) {
  bool getMergeTrackIcpSuccess() {return getMergeTrackIcpSuccess(trackMaxMergeScore.get());}; //!< only works if trackMaxMergeScore contains a valid pointer!
  bool getMergeTrackIcpSuccess(PointCloudTrack *linktrack) {return /*(linktrack->getAge() > getAge()+1) &&*/ (linktrack->histAllStates[linktrack->histAllStates.size()-1-getAge()].wasMeasured);};


private:
  PointCloudTrack(): uid(0) {}; //!< make private so that it cannot be used
  static IdT getNewUid();
  static void initializeFH(const mdefs::DSMatrix &predictionVariance); //!< specifies system and measurement matrix

  // General variables:
  const IdT                 uid; //!< unique id indicating the track
  double                    colorR;
  double                    colorG;
  double                    colorB;
  //unsigned int              idx; //!< index = position in trackList, so all current tracks have index 0..trackcount-1
  unsigned int              age; //!< indicates the age of the track in number of frames
  unsigned int              lastUpdateCounter; //!< in number of frames (reset to 0 when updated)
  bool                      isMovingObj; //!< true if object had at least once a speed > 0 (relative to world!!!)
  float                     colC, rowC; //!< column/row of center at time of creation

  // Kalman Filter variables:
  mdefs::DVector            x; //!< (Kalman-Filter-)State of the current track relative to Ego-CS [rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot]
  mdefs::DSMatrix           P; //!< covariance of the state (size 12 x 12)
  static mdefs::DSMatrix    Q; //!< (constant) prediction-covariance of the state (size 12 x 12)
  static mdefs::DMatrix     F; //!< System matrix for predicting state (size 12 x 12)
  static mdefs::DMatrix     Ft; //!< Transposed system matrix
  static mdefs::DMatrix     H; //!< Measurement matrix (size 12 x 6), 6 = size of measurement vector
  static mdefs::DMatrix     Ht; //!< transposed measurement matrix
  bool                      expectingPrediction; //!< true after update() was called to indicate that prediction should be next

  // Appearance variables (relative to the current state coordinate system):
  // surfaces are stored in both grids via shared-pointers, thus adding/removing points MUST be accounted for in both grids!
  PointCloudPGrid    pointCloudPGrid; //!< 3D point index on the point cloud, uses pointers to the real surface points
  PointCloudNGrid    pointCloudNGrid; //!< 2D normal index on the point cloud via angles, uses pointers to the real surface points
  bool               pointsCleared; //!< flag that is set true when appearance is cleared during merge and set false on prediction
  double weightPProxy(const SurfacePProxy &p) const {
    return p.s->isVIP ? 5.0 : 1.0;
  }
  double weightNProxy(const SurfaceNProxy &n) const {
    return n.s->isVIP ? 5.0 : 1.0;
  }
  std::string               worldmapFilename; //!< filename of p3d-mapfile to write world track into
  bool                      wasCloned; //!< flag whether this track was already cloned from using the copy constructor

  // Frame-based temporary variables:
  //unsigned int               projPixCount; //!< pixels assigned to this track in the current frame by projection (not segmentation)
  //std::set< LidarSegment* > segments; //!< stores which segments this track overlays
  //std::set< LidarSegment* > segLastFrame; //!< stores which segments this track overlayed in the last frame, needed for feature-matching
  
public:
  // Variables used for debugging (can be removed as soon as everything works and as calculateEgoMotion uses the current velocity instead of this history):
  std::vector< StateHistory > histAllStates; //!< contains all the previous states of the track (relative to Ego-CS). current state is at back(). history-size = age+1 (history include init-state)
  mdefs::DVector            histXPredict; // stores the predicted state (which is overwritten when updating)
  mdefs::DMatrix            histPPredict; // Covariance matrix of state
  mdefs::DVector            histZExp; // Expected measurement
  mdefs::DVector            histZ; // Measurement
  mdefs::DMatrix            histR; // Covariance of measurement
  mdefs::DMatrix            histKGain; // Kalman-Gain
  std::list< double >       histRegICPError; // stores the errors of each iteration. first error is used to decide upon FeatureMatching
  bool                      histRegUsedFeatureMatching;
};


template <size_t N>
boost::array<unsigned int, N> PointCloudTrack::getMotion2NormalHistogram(const mdefs::DVector &state1, const mdefs::DVector &state2, double *moveMainNormalDir) const
{
  using namespace matrixTools;
  assert(state1.size() >= 6 && "PointCloudTrack::getMotion2NormalHistogram: state is not >= 6-dimensional");
  assert(state2.size() >= 6 && "PointCloudTrack::getMotion2NormalHistogram: state is not >= 6-dimensional");
  boost::array<unsigned int, N> result; // number of normals per direction (relative to move), [0]=perpendicular
  boost::array<double, N>       avgMove; // projected move onto this direction
  for (size_t idx=0; idx<N; ++idx) { // initialize
    result[idx] = 0;
    avgMove[idx] = 0.0;
  }
  DMatrix R1, R2, R; DVector t1, t2, t;
  state2Rt_T2E(state2, R2, t2); // at state 2, transform point into Ego-CS
  state2Rt_E2T(state1, R1, t1); // at state 1, transform back to Track-CS
  R = ublas::prod(R1, R2); // combine the two transformations
  t = ublas::prod(R1, t2) + t1;
  for (PointConstIterator pi = pointCloudNGrid.begin(); pi != pointCloudNGrid.end(); ++pi) {
    const DVector &point = pi->s->position;
    DVector move = ublas::prod(R, point) + t - point; // move-vector, relative to Track-CS
    double norm = ublas::norm_2(move);
    if (norm > 0) {
      move /= norm; // scale to length
      double cp = ublas::inner_prod(move, pi->s->normal); // cp is now in range -1 to 1, 0 meaning perpendicular
      size_t bin = fabs(cp)*(double)(N-1); // will give a bin in 0..N-1
      result[bin]++;
      avgMove[bin] += cp*norm; // move projected on normal vector
    }
  }
  int maxIdx = 0;
  unsigned int maxVal = 0;
  for (size_t idx=0; idx<N; ++idx) { // make moves average
    avgMove[idx] /= (double)(max(1u,result[idx]));
    if (result[idx] > maxVal) {
      maxVal = result[idx];
      maxIdx = idx;
    }
  }
  if (moveMainNormalDir)
    *moveMainNormalDir = avgMove[maxIdx];
  return result;
};
template <typename InsertIterator>
unsigned int PointCloudTrack::samplePointsUniformly(const unsigned int desiredCount, InsertIterator container) const
{
  std::list<SurfacePProxy> sampledPoints;
  std::back_insert_iterator< std::list<SurfacePProxy> > backInsIter(sampledPoints);
  PointCloudPGrid::weight_func wf = boost::bind(&PointCloudTrack::weightPProxy, this, _1);
  unsigned int nb = pointCloudPGrid.uniform_sample(desiredCount, backInsIter, wf);
  BOOST_FOREACH(const SurfacePProxy &surface, sampledPoints)
    *container++ = surface.s;
  return nb;
};
template <typename InsertIterator>
unsigned int PointCloudTrack::sampleNormalsUniformly(const unsigned int desiredCount, InsertIterator container) const
{
  std::list<SurfaceNProxy> sampledPoints;
  std::back_insert_iterator< std::list<SurfaceNProxy> > backInsIter(sampledPoints);
  PointCloudNGrid::weight_func wf = boost::bind(&PointCloudTrack::weightNProxy, this, _1);
  unsigned int nb = pointCloudNGrid.uniform_sample(desiredCount, backInsIter, wf);
  BOOST_FOREACH(const SurfaceNProxy &surface, sampledPoints)
    *container++ = surface.s;
  return nb;
};

#endif /*POINTCLOUDTRACK_H_*/
