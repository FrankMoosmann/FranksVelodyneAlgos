#include "LidarFrame.hpp"

#include <cfloat>
#include <climits>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include "kogmo_time.h"

using namespace std;
using namespace matrixTools;


LidarFrame::LidarFrame(const unsigned int horizSize, const unsigned int vertSize)
  : distance(horizSize, vertSize)
   ,distanceDiffLastFrame(horizSize, vertSize)
   ,intensity(horizSize, vertSize)
   ,point3D(horizSize, vertSize)
   ,pointVariance3D(horizSize, vertSize)
   ,distDerivativeH(horizSize-1, vertSize)
   ,distDerivativeV(horizSize, vertSize-1)
   ,connectivityH(horizSize-1, vertSize)
   ,connectivityV(horizSize, vertSize-1)
   ,normal3D(horizSize, vertSize)
   ,normalStdDevRAD(horizSize, vertSize)
   ,normalVariance3D(horizSize, vertSize)
   ,normalConfidence(horizSize, vertSize)
   ,segmentCritR(horizSize, vertSize)
   ,segmentCritD(horizSize, vertSize)
   ,segmentCritRU(horizSize, vertSize)
   ,segmentCritRD(horizSize, vertSize)
   ,segmentConnectH(horizSize-1, vertSize)
   ,segmentConnectV(horizSize, vertSize-1)
   ,segmentThreshH(horizSize-1, vertSize)
   ,segmentThreshV(horizSize, vertSize-1)
   ,segmentGrowDecisionH(horizSize-1, vertSize)
   ,segmentGrowDecisionV(horizSize, vertSize-1)
   ,segments(horizSize, vertSize)
   ,sttValidationSteps(3)
   ,projLttTracks(horizSize, vertSize)
   ,projLCTracks(horizSize, vertSize)
   ,projTracks(horizSize, vertSize)
   ,projTrackDist(horizSize, vertSize)
   ,projTrackCnt(horizSize, vertSize)
   ,extProjTracks(horizSize, vertSize)
   ,matchingFeatures(horizSize, vertSize)
   ,matchedColNextFrame(horizSize, vertSize)
   ,matchedRowNextFrame(horizSize, vertSize)
   ,regAvailPixels(horizSize, vertSize)
   ,visibleInNextFrame(horizSize, vertSize)
   ,visibleInPrevFrame(horizSize, vertSize)
   ,tmpVec3D(horizSize, vertSize)
   ,tmpFloat1(horizSize, vertSize)
   ,tmpFloat2(horizSize, vertSize)
   ,tmpFloat3(horizSize, vertSize)
   ,tmpInt1(horizSize, vertSize)
   ,tmpInt2(horizSize, vertSize)
{
  shortTermTracks.resize(sttValidationSteps+2);
  reset();
}

LidarFrame::~LidarFrame()
{
}

size_t LidarFrame::memorySize() const
{
  size_t sum = 0;
  sum += sizeof(boost::mutex);
  sum += sizeof(bool);
  sum += sizeof(int64_t);
  sum += sizeof(double);
  sum += sizeof(double)*positionHTM2v.size1()*positionHTM2v.size2();
  sum += sizeof(double)*positionHTM2w.size1()*positionHTM2w.size2();
  sum += sizeof(double)*egoEstimationHTM2w.size1()*egoEstimationHTM2w.size2();
  sum += sizeof(double);
  // TODO (9): all pointer-types (as DVector) will report a wrong size
  sum += distance.getMemorySize();
  sum += distanceDiffLastFrame.getMemorySize();
  sum += intensity.getMemorySize();
  sum += point3D.getMemorySize();
  sum += pointVariance3D.getMemorySize();
  sum += distDerivativeH.getMemorySize();
  sum += distDerivativeV.getMemorySize();
  sum += connectivityH.getMemorySize();
  sum += connectivityV.getMemorySize();
  sum += normal3D.getMemorySize();
  sum += normalStdDevRAD.getMemorySize();
  sum += normalVariance3D.getMemorySize();
  sum += normalConfidence.getMemorySize();
  sum += segmentCritR.getMemorySize();
  sum += segmentCritD.getMemorySize();
  sum += segmentCritRU.getMemorySize();
  sum += segmentCritRD.getMemorySize();
  sum += segmentConnectH.getMemorySize();
  sum += segmentConnectV.getMemorySize();
  sum += segmentThreshH.getMemorySize();
  sum += segmentThreshV.getMemorySize();
  sum += segmentGrowDecisionH.getMemorySize();
  sum += segmentGrowDecisionV.getMemorySize();
  sum += segments.getMemorySize();
//  sum += trackList.size()*sizeof(..);
  sum += projLttTracks.getMemorySize();
  sum += projLCTracks.getMemorySize();
  sum += projTracks.getMemorySize();
  sum += projTrackDist.getMemorySize();
  sum += projTrackCnt.getMemorySize();
  sum += extProjTracks.getMemorySize();
  // TODO (9): calculate memory consumption of tracks
  sum += matchingFeatures.getMemorySize();
  sum += matchedColNextFrame.getMemorySize();
  sum += matchedRowNextFrame.getMemorySize();
  sum += regAvailPixels.getMemorySize();
  sum += visibleInNextFrame.getMemorySize();
  sum += visibleInPrevFrame.getMemorySize();
  sum += tmpVec3D.getMemorySize();
  sum += tmpFloat1.getMemorySize();
  sum += tmpFloat2.getMemorySize();
  sum += tmpFloat3.getMemorySize();
  sum += tmpInt1.getMemorySize();
  sum += tmpInt2.getMemorySize();
  return sum;
}

void LidarFrame::reset()
{
  valid = false;
  recTime = 0;
  timeDiffSec = 0.0;
  positionHTM2v = ublas::identity_matrix<double>(4);
  positionHTM2w = ublas::identity_matrix<double>(4);
  speed = 0.0;
  distance.fill(DBL_MAX); // set externally
  distanceDiffLastFrame.fill(DBL_MAX); // by AlgoFrameFeatures
  intensity.fill(0.0f); // set externally
  point3D.fill(DVector(3,DBL_MAX)); // by AlgoFrameFeatures
  pointVariance3D.fill(DMatrix(3,3,0.0)); // by AlgoFrameFeatures
  distDerivativeH.fill(DBL_MAX); // by AlgoFrameFeatures
  distDerivativeV.fill(DBL_MAX); // by AlgoFrameFeatures
  connectivityH.fill(0.0f); // by AlgoFrameFeatures
  connectivityV.fill(0.0f); // by AlgoFrameFeatures
  normal3D.fill(DVector(3,0.0)); // by AlgoFrameFeatures
  normalStdDevRAD.fill(M_PI); // by AlgoFrameFeatures
  normalVariance3D.fill(DMatrix(3,3,0.0)); // by AlgoFrameFeatures
  normalConfidence.fill(0.0); // by AlgoFrameFeatures
  segmentCritR.fill(0.0); // by AlgoSegmentation
  segmentCritD.fill(0.0); // by AlgoSegmentation
  segmentCritRU.fill(0.0); // by AlgoSegmentation
  segmentCritRD.fill(0.0); // by AlgoSegmentation
  segmentConnectH.fill(0.0); // by AlgoSegmentation
  segmentConnectV.fill(0.0); // by AlgoSegmentation
  segmentThreshH.fill(0.5); // by AlgoSegmentation
  segmentThreshV.fill(0.5); // by AlgoSegmentation
  segmentGrowDecisionH.fill(false); // by AlgoSegmentation
  segmentGrowDecisionV.fill(false); // by AlgoSegmentation
//  segments.fill(LidarSegment()); // enabled and proxy disabled by AlgoSegmentation
  for (unsigned int bi=0; bi<shortTermTracks.size(); ++bi) {
    shortTermTracks[bi].tracks.clear(); // by AlgoRegistration
    shortTermTracks[bi].linksT1.reset(); // by AlgoRegistration
    shortTermTracks[bi].linksT2.reset(); // by AlgoRegistration
    shortTermTracks[bi].lttLinks.reset(); // by AlgoRegistration
  }
  longTermTracks.clear(); // by AlgoRegistration
  removedTracks.clear(); // by AlgoRegistration
  projLttTracks.fill(UINT_MAX); // by AlgoRegistration
  projLCTracks.fill(UINT_MAX); // by AlgoRegistration
  projTracks.fill(UINT_MAX); // by AlgoRegistration
  projTrackDist.fill(DBL_MAX); // by AlgoRegistration
  projTrackCnt.fill(0); // by AlgoRegistration
  extProjTracks.fill(UINT_MAX); // by AlgoRegistration
  projSttIdx = UINT_MAX; // by AlgoRegistration
//  matchingFeatures.fill(FeatureVector()); // empty FV by AlgoFrameFeatures
  matchedColNextFrame.fill(-1); // by AlgoRegistration
  matchedRowNextFrame.fill(-1); // by AlgoRegistration
  regAvailPixels.fill(true); // by AlgoRegistration
  visibleInNextFrame.fill(true);
  visibleInPrevFrame.fill(true);
  tmpVec3D.fill(DVector(3,0.0)); // borders must be initialized with 3D vectors
}

std::string LidarFrame::getTime()
{
  KogniMobil::kogmo_timestamp_string_t str;
  KogniMobil::kogmo_timestamp_to_string(recTime, str);
  return string(str);
}

// list must not contain NULL pointers!!!
PointCloudTrack::SPtr find(PointCloudTrack::IdT id, const LidarFrame::TrackSPtrList &list) {
  BOOST_FOREACH(PointCloudTrack::SPtr t, list) {
    if (t->getUID() == id)
      return t;
  }
  return PointCloudTrack::SPtr();
}

PointCloudTrack::SPtr LidarFrame::getLttWorldTrack()
{
  PointCloudTrack::SPtr rval;
  if (longTermTracks.size() > 0)
    rval = longTermTracks.front();
  return rval;
}

LidarFrame::TrackLinks LidarFrame::getMergedSTTLinks(unsigned int sstIdx, LidarFrame::TrackIndex ti, bool includeWorldTrack)
{
  typedef pair<TrackIndex, unsigned int>           ConnT;
  typedef std::vector<ConnT>                       ConnTList;
  typedef pair<PointCloudTrack::IdT, unsigned int> LTConnT;
  typedef std::vector<LTConnT>                     LTConnTList;
  LidarFrame::TrackLinks tl;

  PointCloudTrack::SPtr worldTrack = getLttWorldTrack();

  // usually: sstIdx == sttValidationSteps-1
  if ((sstIdx >= sttValidationSteps) || (shortTermTracks[sstIdx].linksT1.get() == NULL)) {
    if (includeWorldTrack && (worldTrack.get() != NULL))
      tl.push_back(TrackLink(-1, worldTrack, 0));
    return tl;
  }

  // 1a) get connections to t-1 tracks
  ConnTList conn = shortTermTracks[sstIdx].linksT1->getConnectionsToL2(ti);
  BOOST_FOREACH(ConnT c, conn) {
    PointCloudTrack::SPtr linkedTrack = shortTermTracks[sstIdx+1].tracks[c.first];
    if (linkedTrack.get() != NULL)
      tl.push_back(TrackLink(c.first, linkedTrack, c.second));
  }
  // 1b) get connections to t-2 tracks
  if (shortTermTracks[sstIdx].linksT2.get() == NULL) {
    cout << " WARN: no T2Links! " << flush;
  } else {
    ConnTList conn = shortTermTracks[sstIdx].linksT2->getConnectionsToL2(ti);
    BOOST_FOREACH(ConnT c, conn) {
      PointCloudTrack::SPtr linkedTrack = shortTermTracks[sstIdx+2].tracks[c.first];
      if (linkedTrack.get() != NULL)
        tl.push_back(TrackLink(UINT_MAX, linkedTrack, c.second));
    }
  }
  // 1c) get connections to long-term-tracks, need to check if it still exists
  if (shortTermTracks[sstIdx].lttLinks.get() == NULL) {
    cout << " WARN: no lttLinks! " << flush;
  } else {
    LTConnTList connltt = shortTermTracks[sstIdx].lttLinks->getConnectionsToL2(ti);
    bool oneValid = false;
    BOOST_FOREACH(ConnT c, connltt) {
      PointCloudTrack::SPtr linkedTrack = find(c.first, longTermTracks);
      if (linkedTrack.get() != NULL) { // track still exists
        tl.push_back(TrackLink(UINT_MAX, linkedTrack, c.second));
        oneValid = true;
      }
    }
    if ((!oneValid) && (connltt.size() > 0)) {
      cout << " WARN: none of " << connltt.size() << " lttLinks valid! " << flush;
    }
  }
  // 2) sort connections by linked track address
  std::sort(tl.begin(), tl.end());
  // 3) merge those entries that refer to the same track
  //    do this by creating a new list and pushing elements that summed up all equivalent links
  LidarFrame::TrackLinks newtl;
  bool foundWorld = false;
  for (unsigned int i=0; i<tl.size(); ++i) {
    if (tl[i].track.get() == NULL)
      continue; // if track was deleted during prediction (only happens at SSTSize-1)
    if (tl[i].track == worldTrack)
      foundWorld = true;
    if ((newtl.size()==0) || (tl[i].track.get() != newtl.back().track.get()))
      newtl.push_back(tl[i]);
    else
      newtl.back().linkStrength += tl[i].linkStrength;
  }
  if (!foundWorld && includeWorldTrack && (worldTrack.get() != NULL)) {
    newtl.push_back(TrackLink(-1, worldTrack, 0));
  }
  return newtl;
}


void LidarFrame::checkConsistency() 
{
  cout << endl << "---------------------------------------------------------------------------";
  cout << endl << "checking frame consistency at time " << getTime();

  int hsize, vsize;
  segments.getSize(hsize,vsize);
  // pass over image: check each segment pointer
  unsigned int imgProjPixCount = 0;
  unsigned int imgSegProxyCount = 0;
  unsigned int imgSegLeafCount = 0;
  unsigned int segRefPixCount = 0;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      if (projTracks.get(col,row) != UINT_MAX)
       ++imgProjPixCount;
      LidarSegment &s = segments.get(col,row);
      if (s.isValid()) {
        if (distance.get(col,row) == DBL_MAX)
          cout << endl << "ERROR: pixel " << col << "," << row << " is invalid but its segment-object valid!";
        if (s.isProxy())
          ++imgSegProxyCount;
        else {
          ++imgSegLeafCount;
          segRefPixCount += s.getPxCount();
        }
      }
    }
  }
  if (imgSegProxyCount+imgSegLeafCount != segRefPixCount)
    cout << endl << "ERROR: image has " << imgSegProxyCount+imgSegLeafCount << " valid segment-objects, but " << segRefPixCount << " pixels referenced by segments!";
  if ((int)segRefPixCount > vsize*hsize)
    cout << endl << "ERROR: image has " << vsize*hsize << " pixels, but segments reference " << segRefPixCount << " pixels!";

  unsigned int trkCount[shortTermTracks.size()];
  unsigned int trkNULLCount[shortTermTracks.size()];
  unsigned int trkPtCount[shortTermTracks.size()];
  double maxPosCovar = 0.0;
  double maxPosCovarPtCount = 0;
  unsigned int maxPosCovarBi = 0;
  unsigned int maxPosCovarTi = 0;
  for (unsigned int bi=0; bi<shortTermTracks.size(); ++bi) {
    trkCount[bi] = shortTermTracks[bi].tracks.size();
    trkNULLCount[bi] = 0;
    trkPtCount[bi] = 0;
    unsigned int ti = 0;
    BOOST_FOREACH(PointCloudTrack::SPtr track, shortTermTracks[bi].tracks) {
      ti++;
      if (track.get() == NULL) {
        trkNULLCount[bi]++;
        continue; // might happen for bi >= sttValidationSteps
      }
      trkPtCount[bi] += track->getPointCount();
      DMatrix tc;
      track->getCovar(tc);
      double trackmax = max(max(tc(3,3),tc(4,4)),tc(5,5));
      if (trackmax > maxPosCovar) {
        maxPosCovar = trackmax;
        maxPosCovarPtCount = track->getPointCount();
        maxPosCovarBi = bi;
        maxPosCovarTi = ti;
      }
      if ((track->mergeDecision != PointCloudTrack::MDKeep) && (track->trackMaxMergeScore.get() == NULL))
        cout << endl << "ERROR: track [" << track->getUID() << "] at sst[" << bi << "][" << ti-1 << "] has mergeDecision " << track->mergeDecision << " but NULL-link";
    }
  }
  PointCloudTrack::SPtr worldTrack = getLttWorldTrack();

  cout << endl << "- image has " << hsize*vsize << " pixels, " << imgProjPixCount << " from projected tracks [of buffer " << projSttIdx << "] and " << segRefPixCount << " from segments";
  cout << endl << "- image has " << imgSegLeafCount<< " valid segments (referenced by " << imgSegProxyCount << " proxies -> " << imgSegLeafCount+imgSegProxyCount << " valid segment objects in total)";
  cout << endl << "- frame has " ; for (unsigned int bi=0; bi<shortTermTracks.size(); ++bi) {if (bi>=sttValidationSteps) cout<<"{"; cout << trkCount[bi];     if (bi>=sttValidationSteps) cout<<"}"; cout << "+";} cout << "(ltt:" << longTermTracks.size() << ") tracks";
  cout << endl << "       with " ; for (unsigned int bi=0; bi<shortTermTracks.size(); ++bi) {if (bi>=sttValidationSteps) cout<<"{"; cout << trkPtCount[bi];   if (bi>=sttValidationSteps) cout<<"}"; cout << "+";} cout << "(ltt:?)"; if (worldTrack.get() != NULL) cout << "+(wt:" << worldTrack->getPointCount() << ") points";
  cout << endl << "    thereof " ; for (unsigned int bi=0; bi<shortTermTracks.size(); ++bi) {if (bi>=sttValidationSteps) cout<<"{"; cout << trkNULLCount[bi]; if (bi>=sttValidationSteps) cout<<"}"; cout << "+";} cout << "(ltt:?) tracks invalid";
  cout << endl << "- maximum uncertain short-term-track at stt[" << maxPosCovarBi << "][" << maxPosCovarTi << "]: " << maxPosCovarPtCount << " points with covar " << maxPosCovar;
  cout << endl << "---------------------------------------------------------------------------";
  cout << flush;
}
