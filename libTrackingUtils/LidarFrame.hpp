#ifndef LIDARFRAME_H_
#define LIDARFRAME_H_

#include <vector>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include "MatrixDefs.hpp"
#include "LidarImage.hpp"
#include "LidarSegment.hpp"
#include "PointCloudTrack.hpp"
#include "TwoLayerGraph.hpp"

class LidarFrame;
typedef boost::shared_ptr< LidarFrame > FrameSPtr;

namespace mdefs = matrixTools;
/*! 
 * \class LidarImageFrame
 * \brief This is the core class containing all data for signal processing on one Lidar Frame
 * 
 * An object of this class represents one Lidar frame and all associated data. More specifically it holds
 * - The associated position of the car in Cartesian absolute coordinates (in form of a HTM)
 * - The Lidar data itself as (1) distance image, (3) absolute Cartesian images for x,y,z respectively
 * - The first derivative in horiz and vert direction of the distance values
 * - A "pixel connectivity" image
 * - Associated normal vectors
 * - And more....
 */
class LidarFrame
{
public:
  LidarFrame(const unsigned int horizSize, const unsigned int vertSize);
  virtual ~LidarFrame();
  size_t memorySize() const;

public:
  // TODO(9): use frame-mutex everywhere!!!
  boost::mutex                           accessMutex; //!< this should be locked as soon as a frame is accessed. If processing accesses 2 frames, always lock the most current frame first! If all routines adhere to this guideline, no deadlocks will occur !!!
  bool                                   valid; //!< indicates if at least basic data (recTime, pos, speed, distance+intensity) is valid
  std::string                           sourcefile; //!< filename if available
  int64_t                               recTime; //!< timestamp of recording time (in nanoseconds)
  double                                timeDiffSec; //!< timedifference to last frame (in seconds)
//#define KOGMO_TIMESTAMP_TICKSPERSECOND 1000000000
//#define KOGMO_TIMESTAMP_NANOSECONDSPERTICK (1000000000/KOGMO_TIMESTAMP_TICKSPERSECOND)
  // egoMotionFilter
  mdefs::DMatrix                        positionHTM2v; //!< vehicle position as HTM. point in world cs pw can be transformed into vehicle cs by pv=HTM*pw
  mdefs::DMatrix                        positionHTM2w; //!< vehicle position as HTM. point in vehicle cs pv can be transformed into world cs by pw=HTM*pv
  mdefs::DMatrix                        egoEstimationHTM2w; //!< vehicle position as HTM. point in vehicle cs pv can be transformed into world cs by pw=HTM*pv
  double                                speed; //!< current vehicle speed in m/s
  LidarImage<double>                    distance; //!< distance in meter, DBL_MAX indicating invalid measurement
  LidarImage<double>                    distanceDiffLastFrame; //!< distance-difference to last frame, DBL_MAX indicating invalid measurement
  LidarImage<double>                    intensity; //!< reflectivity, values between 0.0 and 1.0
  LidarImage<mdefs::DVector>            point3D; //!< relative ego-based coordinates in meter, (DBL_MAX,DBL_MAX,DBL_MAX) indicating invalid measurement (in that case distance is also DBL_MAX)
  LidarImage<mdefs::DMatrix>            pointVariance3D; //!< the covariance of the point coordinate (invalidity indicated by point3D)
  LidarImage<double>                    distDerivativeH; //!< distance differences in meter, DBL_MAX if invalid
  LidarImage<double>                    distDerivativeV;
  LidarImage<double>                    connectivityH; //!< weighted connections, values between 0.0 and 1.0
  LidarImage<double>                    connectivityV;
  LidarImage<mdefs::DVector>            normal3D; //!< relative ego-based coordinates in meter, (0,0,0) indicating invalid measurement
  LidarImage<double>                    normalStdDevRAD; //!< the standard deviation of the angle of the normal vector in radian (zero by default)
  LidarImage<mdefs::DMatrix>            normalVariance3D; //!< the covariance of the normal vector's 3D coordinates (invalidity indicated by normal3D)
  LidarImage<double>                    normalConfidence; //!< a value between 0 (unreliable/not valid) and 1 (reliable) indicating how precise the normal3D estimate is. Can be interpreted as c*exp(-covar)
  LidarImage<double>                    segmentCritR; //!< segmentation criterion to the right, values between 0.0 and 1.0
  LidarImage<double>                    segmentCritD; // down
  LidarImage<double>                    segmentCritRU; // right and up
  LidarImage<double>                    segmentCritRD; // right and down
  LidarImage<double>                    segmentConnectH; //!< weighted connections derived from segmentCrit for segmentation (will be combined with connectivityH/V), values between 0.0 and 1.0
  LidarImage<double>                    segmentConnectV;
  LidarImage<double>                    segmentThreshH; //!< threshold used to decide about growing segment (track-ids might locally modify this otherwise constant threshold)
  LidarImage<double>                    segmentThreshV;
  LidarImage<bool>                      segmentGrowDecisionH; //!< decision of growing
  LidarImage<bool>                      segmentGrowDecisionV;
  LidarImage<LidarSegment>              segments; //!< pixel-wise labeling via pointer to segment object
  typedef std::list<PointCloudTrack::SPtr>                 TrackSPtrList;
  typedef std::vector<PointCloudTrack::SPtr>               TrackSPtrVec;
  typedef unsigned int                                     TrackIndex; //!< index in TrackSPtrVec
  typedef TwoLayerGraph<TrackIndex, TrackIndex>            TrackLinkageGraph;
  typedef TwoLayerGraph<TrackIndex, PointCloudTrack::IdT>  LTTrackLinkageGraph;
  typedef boost::shared_ptr< TrackLinkageGraph >           TrackLinkageGraphSPtr;
  typedef boost::shared_ptr< LTTrackLinkageGraph >         LTTrackLinkageGraphSPtr;
  struct                                         ExtTrackList {
    TrackSPtrVec              tracks;  //!< list of tracks
    TrackLinkageGraphSPtr     linksT1; //!< from "tracks" (above) to associated "tracks" in stt[+1] (derived from proj-tracks at creation time)
    TrackLinkageGraphSPtr     linksT2; //!< from "tracks" (above) to associated "tracks" in stt[+2] (derived from proj-tracks at creation time)
    LTTrackLinkageGraphSPtr   lttLinks; //!< from "tracks" (above) to associated ltt-tracks (derived from proj-tracks at creation time)
  };
  /*  with N=sttValidationSteps:
   *  and hence sttValidationSteps.size()==N+2
   *  shortTermTracks[0]       tracks created from segmentation (at time t), links refer to [1]+[2]
   *  shortTermTracks[1]       tracks from t-1 (registered), links refer to [2]&[3]
   *  shortTermTracks[2]       tracks from t-2 (registered), links refer to [3]&[4]
   *  ...
   *  shortTermTracks[N-1]     tracks from t-N+1 (registered), links refer to [N]&[N+1]
   *  -------------------------
   *  shortTermTracks[N]       tracks from t-N, already merged with LTTS. buffer only kept such that linkage graph from [N-1] has valid targets
   *  shortTermTracks[N+1]     tracks from t-N-2, already merged with LTTS. buffer only kept such that second linkage graph from [N-1] has valid targets
   *            both [N]&[N+1] are multi-subset of LTTs, i.e.
   *                           A) several entries might be equal (if original tracks were merged into the same track)
   *                           B) valid entries (!=NULL) are contained in LTTs
   *                           C) invalid entries (==NULL) derive from tracks that were deleted from LTTs
   *  this description is valid only at the END of processing!
   */
  const unsigned int                             sttValidationSteps; //!< number of validation steps before merging
  std::vector<ExtTrackList>                      shortTermTracks; //!< tracks created at [t_rel]. [0] from the current measurements, [1] from the last frame...
  TrackSPtrList                                  longTermTracks; //!< tracks that were verified over several frames, these tracks can be regarded as "valid". [0]=worldTrack
  PointCloudTrack::SPtr                          getLttWorldTrack(); //!< extracts world track from longTermTracks
  TrackSPtrList                                  removedTracks; //!< long-term-tracks that were removed in this frame (list exists for debugging purposes)
  LidarImage<TrackIndex>                         projLttTracks; //!< vector-index (stt[0]::lttTracksCopy) of projected tracks
  LidarImage<TrackIndex>                         projLCTracks; //!< vector-index (stt[+2]::tracks) of projected tracks
  LidarImage<TrackIndex>                         projTracks; //!< vector-index (stt[+1]::tracks) of projected tracks
  LidarImage<double>                             projTrackDist; //!< distance of projected track
  LidarImage<unsigned int>                       projTrackCnt; //!< number of tracks that project to same pixel (0,1,2), 2 meaning >=2
  LidarImage<TrackIndex>                         extProjTracks; //!< list-index of projected predicted track
  unsigned int                                   projSttIdx; //!< corresponding short-term-track buffer index of projected tracks
  struct                                         TrackLink {
    TrackIndex                index; // position in shortTermTracks[sstIdx+1] when called getMergedSTTLinks(sstIdx,..)
    PointCloudTrack::SPtr     track;
    TrackLinkageGraph::ValT   linkStrength;
    TrackLink(TrackIndex index_, PointCloudTrack::SPtr track_, TrackLinkageGraph::ValT linkStrength_)
      : index(index_), track(track_), linkStrength(linkStrength_) {};
    bool operator< (const TrackLink &other) const {
      return (track->getUID() < other.track->getUID()); // sort by UID, i.e. by object, automatically sorts per creation-time (world track 1st)
    }
  };
  typedef std::vector<TrackLink>                 TrackLinks;  //!< can be used for retrieving links for a given track
  TrackLinks                                     getMergedSTTLinks(unsigned int sstIdx, TrackIndex ti, bool includeWorldTrack = false); // same as "shortTermTracks[sstIdx].links->getConnectionsToL2(ti)" but with merging of results refering to the same track
  LidarImage<mdefs::DVector>            matchingFeatures; //!< some feature vector for each pixel (empty if not valid)
  LidarImage<int>                       matchedColNextFrame; //!< Matched pixel in next frame (-1/-1 indicating no matching found)
  LidarImage<int>                       matchedRowNextFrame;
  LidarImage<bool>                      regAvailPixels; //!< masks, which pixels are available for registration. once a track is registered, all corresponding pixels are masked, so they are not used further
  LidarImage<bool>                      visibleInNextFrame; //!< true if pixel is visible in next frame (static scene assumed)
  LidarImage<bool>                      visibleInPrevFrame; //!< true if pixel is visible in next frame (static scene assumed)
  LidarImage<mdefs::DVector>            tmpVec3D;
  LidarImage<double>                    tmpFloat1;
  LidarImage<double>                    tmpFloat2;
  LidarImage<double>                    tmpFloat3;
  LidarImage<int>                       tmpInt1;
  LidarImage<int>                       tmpInt2;
  
  void checkConsistency();

  std::string getTime();

private: // hide, so they cannot be used
  LidarFrame()
  : sttValidationSteps(0)
  {};
  LidarFrame& operator=(const LidarFrame &other) {(void)other; return *this;};
  void reset(); // set default values to all entries

};

#endif /*LIDARFRAME_H_*/
