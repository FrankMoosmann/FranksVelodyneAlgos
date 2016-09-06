#ifndef ALGOREGISTRATION_H_
#define ALGOREGISTRATION_H_

#include <set>

#include <ANN/ANN.h> // (approximate) nearest neighbor - for feature matching
#include "kdtree2.hpp"
#include "LidarFrame.hpp"
#include "PointCloudTrack.hpp"
#include "LidarImageProjector.hpp"
#include "ICPLinearized.hpp"

namespace mdefs = matrixTools;

/*!
 * \class SearchResultsBase
 * \brief A Base class used for ICP
 */
class SearchResultsBase {
private:
  SearchResultsBase *sr; // this base class must only be instantiated for refering to a derived object
public:
  SearchResultsBase(SearchResultsBase* sr_); // takes memory responsibility
  virtual ~SearchResultsBase();
  virtual bool next(); //!< increments iterator and returns false in case end is reached
  virtual mdefs::DVector getCorrespPoint(); //!< \throws std::exception if invalid correspondence
  virtual mdefs::DVector getCorrespNormal(); //!< \throws std::exception if invalid correspondence
  virtual double getNormConf(); //!< \throws std::exception if invalid correspondence
  virtual double getNormStdDevRAD(); //!< \throws std::exception if invalid correspondence
  virtual mdefs::DSMatrix getPCovar(); //!< \throws std::exception if invalid correspondence
  virtual mdefs::DSMatrix getNCovar(); //!< \throws std::exception if invalid correspondence
  virtual double getDist(); //!< \throws std::exception if invalid correspondence
  virtual double getWeight(); //!< \throws std::exception if invalid correspondence
};

/*!
 * \class ProjectedNeighborSearch
 * \brief allows to search for closest point, including projection onto local surface
 *
 * This class is used by the registration procedure below.
 * It allows for efficient neighbor search among the 3D points of the provided frame.
 * The first method searches the closest neighbor,
 * the second projects the search point onto the local surface of the closest neighbor
 * and returns the projected point.
 * Both additionally check whether the neighbor was already assigned by using
 * the binary mask "regAvailPixels" of the frame provided in the constructor
 * (throwing an exception in case it is not available)
 */
class ProjectedNeighborSearch {
public:
  ProjectedNeighborSearch(FrameSPtr frame, const LidarImageProjector *projector = 0);
  virtual ~ProjectedNeighborSearch();

  void useNormalWeighting(bool use) {useNormalWeight = use;}; //init: params->regIcpUseNormalWeight
  void useOcclusionWeighting(bool use) {useOcclWeight = use;}; //init: params->regIcpUseOcclWeight
  void useDistDiffWeighting(bool use) {useDistDiffWeight = use;}; //init: params->regIcpUseDistDiffWeight
  void useNormalSearch(bool use) {kdTree = use ? kdTree6D : kdTree3D;}; //init: kdTree6D
  void useHorizSearch(bool use) {kdTree = use ? kdTreeStretch : kdTree3D;};
  void usePixMaskImage(bool use) {usePixMask = use;};

  class SearchResults : public SearchResultsBase {
  public:
    enum ProjMode {PROJECT, DONTPROJECT};
    SearchResults(const ProjectedNeighborSearch &pns, ProjMode pmode_) : SearchResultsBase(NULL), ns(pns), pmode(pmode_), index(0), correspValid(false) {};
    virtual ~SearchResults() {};

    virtual bool next();
    virtual mdefs::DVector getCorrespPoint();
    virtual mdefs::DVector getCorrespNormal();
    virtual double getNormConf();
    virtual double getNormStdDevRAD();
    virtual mdefs::DSMatrix getPCovar();
    virtual mdefs::DSMatrix getNCovar();
    virtual double getDist();
    virtual double getWeight();
  private:
    const ProjectedNeighborSearch &ns;
    ProjMode pmode;
    unsigned int index;
    void assertCorrespValid(); //!< calculates correspondence point+normal, in case they are not yet valid.  \throws std::exception if invalid correspondence
    bool correspValid;
    mdefs::DVector correspPoint;
    mdefs::DVector correspNormal;
    mdefs::DSMatrix correspPCovar;
    mdefs::DSMatrix correspNCovar;
    double correspNormConf;
    double correspNormStdDev;
  };
  friend class SearchResults; // allow SearchResults access to private variables of ProjectedNeighborSearch

  //! \returns SearchResults object containing nnCount neighbors for given search point + normal
  SearchResultsBase findClosestNeighbor(const mdefs::DVector &searchPoint, const mdefs::DVector &searchNormal, double normalConfidence, unsigned int nnCount) const;
  //! \returns SearchResults object containing nnCount neighbors for given search point + normal
  SearchResultsBase findClosestProjectedNeighbor(const mdefs::DVector &searchPoint, const mdefs::DVector &spNormal, double normalConfidence, unsigned int nnCount) const;
  //! public buffer to store all correspondences that are searched for
  mutable std::list<mdefs::DVector> correspondenceBuffer; // for DEBUG purposes
private:
  // status variables:
  FrameSPtr             frame;
  const LidarImageProjector *projector;
  double                maxDist;
  double                weightZ;
  double                 projDistTolerance;
  bool                  useNormalWeight;
  bool                  useDistDiffWeight;
  bool                  useOcclWeight;
  bool                  usePixMask;
  double                occlWeight;
  // kd-Tree storage:
  boost::multi_array<float,2> dataPts;
  boost::multi_array<float,2> dataPtsStretch;
  int                   *idxMap;
  kdtree2::KDTree       *kdTree;
  kdtree2::KDTree       *kdTree6D;
  kdtree2::KDTree       *kdTree3D;
  kdtree2::KDTree       *kdTreeStretch;
  // query variables, might even be changed in "const" retrieval-functions:
  mutable mdefs::DVector searchPoint; // buffered query point
  mutable mdefs::DVector searchNormal; // buffered query normal
  mutable double pixPlaneDist; // distance of query point to scene, positive if in front of scene point/plane
  mutable kdtree2::KDTreeResultVector result; // query-result: index
  mutable std::vector<float> queryPt; // query point, transformed in kd-tree format
  // private helper functions:
  void findNeighbors(const mdefs::DVector &searchPoint, const mdefs::DVector &searchNormal, double normalConfidence, unsigned int nnCount) const; //!< calls the kd-tree search function, storing results in the kd-tree query arrays
  void getColRow(unsigned int index, int &col, int &row) const; //!< access the kd-tree query result at the specified index and convert into col/row
};

/*!
 * \class TrackNeighborSearch
 * \brief allows to search for closest points in PointCloudTracks
 *
 * This class is used by the registration procedure below.
 */
class TrackNeighborSearch {
public:
  TrackNeighborSearch(PointCloudTrack::SPtr track);
  virtual ~TrackNeighborSearch();

  class SearchResults : public SearchResultsBase {
  public:
    SearchResults(const TrackNeighborSearch &tns) : SearchResultsBase(NULL), ns(tns), currResIndex(0) {};
    virtual ~SearchResults() {};

    virtual bool next();
    virtual mdefs::DVector getCorrespPoint();
    virtual mdefs::DVector getCorrespNormal();
    virtual double getNormConf();
    virtual double getNormStdDevRAD();
    virtual mdefs::DSMatrix getPCovar();
    virtual mdefs::DSMatrix getNCovar();
    virtual double getDist();
    virtual double getWeight();
  private:
    const TrackNeighborSearch &ns;
    mutable unsigned int currResIndex;
  };
  friend class SearchResults; // allow SearchResults access to private variables of TrackNeighborSearch

  //! \returns this as SearchResults object containing nnCount neighbors for given search point + normal
  SearchResultsBase findNeighbors(const mdefs::DVector &searchPoint, const mdefs::DVector &searchNormal, double normalConfidence, unsigned int nnCount) const;

private:
  // input variables:
  double                                 maxDist;
  std::vector<PointCloudTrack::Surface*> surfaces;
  // kd-Tree storage:
  boost::multi_array<float,2> dataPts;
  kdtree2::KDTree             *kdTree;
  // query variables, might even be changed in "const" retrieval-functions:
  mutable mdefs::DVector searchPoint; // buffered query point
  mutable mdefs::DVector searchNormal; // buffered query normal
  mutable std::vector<float> queryPt; // query point, transformed in kd-tree format
  mutable kdtree2::KDTreeResultVector result; // query-result: index
};

/*!
 * \class FeatureMatching
 * \brief allows to search for points closest in feature space
 */
class FeatureMatching {
public:
  FeatureMatching(FrameSPtr currFrame, FrameSPtr lastFrame);
  virtual ~FeatureMatching();
  //! returns col,row,dist of found correspondence \throws exception if no correspondence is found
  void findClosestNeighbor(int colL, int rowL, int &colC, int &rowC, double &dist) const;
private:
  // status variables:
  FrameSPtr     cf;
  FrameSPtr     lf;
  double         featFac;
  double         vertFac;
  double         maxDistanceSqr;
  // kd-Tree storage:
  ANNpointArray dataPts;
  int           *idxMap;
  ANNkd_tree    *kdTree;
  // query variables, might even be changed in const-functions:
  mutable ANNpoint      queryPt; // query point
  mutable ANNidxArray   nnIdx; // query-result: index
  mutable ANNdistArray  dists; // query-result: distance of closest neighbor
};

/*!
 * \class TrackRegistration
 * \brief can be used to register all tracks of one frame with the next frame
 */
class TrackRegistration {
public:
  typedef boost::shared_ptr< TrackRegistration > SPtr;
  TrackRegistration(FrameSPtr currentF, FrameSPtr lastF, const LidarImageProjector &projector);
  virtual ~TrackRegistration();

  bool nextTrack(); // returns false if nothing to be done
  FrameSPtr getFrame() {return cf;};
  std::pair<double,double> calculateError(PointCloudTrack::SPtr track, const mdefs::DVector &state); // point-plane and point-point-error
  void usePixMask(bool use = true) {pnSearch.usePixMaskImage(use);};
  void flushTrackProjListAndExtend(); // method to project all remaining tracks, not projected yet (unsuccessful ICP registrations)
  void extendProjTrackIDs(); // interpolate missing pixels and increase borders

//private:
// for debugging purposes made public
  struct TrackSubsampled {
    PointCloudTrack::SPtr track;
    int trackBufferIdx;
    int trackIdx;
    unsigned int trackPtCnt; 
    mdefs::DCMatrix pointsTrackCS;
    mdefs::DCMatrix pointNormalsTrackCS;
    std::vector<mdefs::DMatrix> pointCovarCS;
    std::vector<mdefs::DMatrix> normalCovarCS;
    mdefs::DVector initialWeights;
    mdefs::DVector normalConfidence;
    mdefs::DVector normalStdDevRAD;
    mdefs::DMatrix RInit;
    mdefs::DVector tInit;
    // used for debugging ICP: one correspondence-list for each iteration
    // each inner list contains searchPoint-correspPoint pairs
    std::list< std::list<mdefs::DVector> > icpCorrespondenceBuffer; // TODO (9): remove this debug-code
  };
  typedef boost::shared_ptr< TrackSubsampled > SubTrackSPtr;
  struct SubTrackSPtrComp {
    bool operator()(const SubTrackSPtr &p1, const SubTrackSPtr &p2) const {
      if (p1->trackBufferIdx == p2->trackBufferIdx)
        return (p1->trackPtCnt > p2->trackPtCnt);
      else
        return (p1->trackBufferIdx > p2->trackBufferIdx);
    };
  };
  typedef std::multiset<SubTrackSPtr,SubTrackSPtrComp> SubTrackSet;
  typedef std::list<SubTrackSPtr> SubTrackLst;

  SubTrackLst subsampledTracks; //!< holds the subsampled tracks where the registration is performed on
  SubTrackLst::iterator nextTrackToRegister;
  SubTrackLst::iterator currTrackRegistered;
  int lastTrackBufferIdx;
  unsigned int nbRegistered;

private:
  static ICP::EnergyFunction::Surface convert(const PointCloudTrack::Surface &s);
  static ICP::EnergyFunction::Surface convert(SearchResultsBase &sr);
  static SubTrackSPtr generate(PointCloudTrack::SPtr track, int bufferIndex, int trackIndex, double minSamplingPtCount);
  static SubTrackSPtr generateFromWorld(FrameSPtr cf, PointCloudTrack::SPtr track, int bufferIndex, int trackIndex, double maxDist);
  void projectToImage(SubTrackSPtr subtrack, bool icpSuccessful);
  std::list<SubTrackSPtr> trackProjList; // tracks with unsuccessful ICP registration, will be projected at last
  FrameSPtr cf;
  FrameSPtr lf;
  const LidarImageProjector &projector;
  ProjectedNeighborSearch pnSearch; // search in current frame
  //TrackNeighborSearch wtSearch; // search in world track
  FeatureMatching *fMatcher; // pointer as it might not be instantiated
  ICP::Point2PlaneEnergy enPtPl;
  ICP::Point2PointEnergy enPtPt;
  ICP::NormalConfLinCombEnergy enComb;
  const double matchFeatFac;
  const double matchVertFac;
  const double matchMaxDistance;
  const bool useFMatch;
  const double fMatchAvgDistThresh;
  const double icp4DPtCountThresh;
  const double icpAvgDistThresh;
  const double icpRegCoeff;
  const unsigned int icpNbCorresp;
  const double icpErrOk;
  const double icpMahalDistState;
  const double projDistTolerance;
  const double lidarMeasStdDevPerMeter;
  const double maxSampleDist;
};


/*!
 * \brief method to match features of two frames. Assumes that both frames have valid features associated
 * 
 * \param cf  Current frame which is matched against next frame
 * \param nf  Relative frame serving as matching basis
 */
void calculateMatches(FrameSPtr cf, FrameSPtr nf);

/*!
 * \brief method to calculate Registration for all segments of a frame
 * 
 * \param cf  Current frame which is matched against next frame
 * \param nf  Relative frame serving as matching basis
 * \param vertFac  set this > 1 to weight z-coordinate, thus prefering neighbors with similar height 
 */
//void calculateRegistration(FrameSPtr cf, FrameSPtr nf, const double vertFac = 2.0f, const unsigned int minSamplingPtCount = 20, const double avgDistThresh = 0.05, const double projDistTolerance = 0.1);



#endif /*ALGOREGISTRATION_H_*/
