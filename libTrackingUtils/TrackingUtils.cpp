#include "TrackingUtils.hpp"

#include <cmath>
#include <cstdlib>
#include <set>
#include <ANN/ANN.h> // (approximate) nearest neighbor
#include <png++/png.hpp>

#include "LidarFrame.hpp"
#include "MatrixDefs.hpp"

using namespace std;
using namespace matrixTools;

void saveToFile(const LidarImage<LidarSegment> &segments, std::string filename)
{
  int hsize, vsize;
  segments.getSize(hsize,vsize);
  png::image<png::gray_pixel_16> image(hsize,vsize);
  map<const LidarSegment*,png::gray_pixel_16> segToPixValMap;
  segToPixValMap[NULL] = 0;
  png::gray_pixel_16 currPixVal = 1;
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      if (segments(col,row).isValid()) {
        BOOST_AUTO(itp,segToPixValMap.insert(make_pair(segments(col,row).getPtr(), currPixVal)));
        if (itp.second) ++currPixVal; // element was newly inserted
        image[row][col] = itp.first->second;
      } else {
        image[row][col] = 0;
      }
    }
  }
  image.write(filename);
  cout << endl << "saved " << segToPixValMap.size()-1 << " segments" << flush;
}

void loadFromFile(LidarImage<LidarSegment> &segments, std::string filename)
{
  int hsize, vsize;
  segments.getSize(hsize,vsize);
  png::image<png::gray_pixel_16> image;
  image.read(filename);
  if (((int)image.get_height() != vsize) || ((int)image.get_width() != hsize))
    throw std::range_error("image sizes don't match");
  map<png::gray_pixel_16,LidarSegment*> pixValToSegMap;
  pixValToSegMap[0] = NULL;
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      segments(col,row).reset();
    }
  }
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      if (image[row][col] == 0) {
        segments(col,row).setValid(false);
      } else {
        segments(col,row).setValid(true);
        BOOST_AUTO(itp,pixValToSegMap.insert(make_pair(image[row][col], segments(col,row).getPtr())));
        itp.first->second->merge(segments(col,row).getPtr());
        segments(col,row).addPix(ColRowPair(col,row));
        // missing: neighborhood relation
      }
    }
  }
  cout << endl << "loaded " << pixValToSegMap.size()-1 << " segments" << flush;
}

void distanceFromPng(LidarImage<double> &distance, const PngDistanceImage &dstImg)
{
  int hsize, vsize;
  distance.getSize(hsize,vsize); // should be representative to take resolution information out of distance image
  if (hsize != (int)dstImg.width())
    throw range_error("TrackingUtils::loadFromFile: objects have different horizontal resolution: "+to_string((long long int)hsize)+"/"+to_string((long long int)dstImg.width()));
  if (vsize != (int)dstImg.height())
    throw range_error("TrackingUtils::loadFromFile: objects have different vertical resolution"+to_string((long long int)vsize)+"/"+to_string((long long int)dstImg.height()));
  // Distance Image Generation
  distance.fill(DBL_MAX);
  double *distMemBlock = distance.getAdr(0, 0);
  dstImg.getDistances(distMemBlock, 2); // two bytes pixel boundary (one left, one right)
  double *cDst = distMemBlock;
  for (int i=0; i<hsize*vsize; ++i, ++cDst)
    if (*cDst == 0.0) *cDst = DBL_MAX;
  // Now, invalid distance-values are MAX and invalid intensity values are 0.0
}

///////////////////////////////////////////////////////
////////////////    Helper Class      /////////////////
///////////////////////////////////////////////////////
/*!
 * \class ClosestSegmentSearch
 * \brief allows to search for the segment which is closest to some point
 */
class ClosestSegmentSearch {
public:
  ClosestSegmentSearch(LidarFrame *frame, set<LidarSegment*>::iterator first, set<LidarSegment*>::iterator end);
  virtual ~ClosestSegmentSearch();
  LidarSegment* findClosestSegment(const DVector &searchPoint);
private:
  ANNpointArray dataPts;
  LidarSegment  **idxMap;
  ANNkd_tree    *kdTree;
  ANNpoint      queryPt; // query point
  ANNidxArray   nnIdx; // query-result: index
  ANNdistArray  dists; // query-result: distance of closest neighbor
};
ClosestSegmentSearch::ClosestSegmentSearch(LidarFrame *frame, set<LidarSegment*>::iterator first, set<LidarSegment*>::iterator end)
{
  // copy 3D data into kd-Tree specific variables
  int nDim = 3;                                     // dimensionality of Search space, x,y,z
  unsigned int nPts = 0;
  BOOST_AUTO(curr,first);
  while (curr!=end) { // get number of points
    LidarSegment* seg = *curr;
    nPts += seg->getPxCount();
    ++curr;
  };
  dataPts = annAllocPts(nPts, nDim);// data points
  idxMap = new LidarSegment*[nPts];
  unsigned int i = 0;
  LidarSegment::ColRowIterator crCurr; LidarSegment::ColRowIterator crEnd; unsigned int crCount;
  for (curr=first; curr!=end; ++curr) {
    LidarSegment* seg = *curr;
    seg->getColRows(crCurr, crEnd, crCount);
    while (crCurr != crEnd) {
      DVector &p = frame->point3D.get(crCurr->first,crCurr->second);
      dataPts[i][0] = p(0);
      dataPts[i][1] = p(1);
      dataPts[i][2] = p(2);
      idxMap[i] = seg;
      ++i;
      ++crCurr;
    }
  }
  // Build kd-Tree
  kdTree = new ANNkd_tree(dataPts, nPts, nDim);
  // Allocate Variables for Queries
  queryPt = annAllocPt(nDim);         // query point
  nnIdx = new ANNidx[1];              // nearest neighbor index for 1 result
  dists = new ANNdist[1];             // nearest neighbor dists for 1 result
}
ClosestSegmentSearch::~ClosestSegmentSearch()
{
  delete [] dists;
  delete [] nnIdx;
  delete kdTree;
  delete [] idxMap;
  annClose(); // done with ANN
}

LidarSegment* ClosestSegmentSearch::findClosestSegment(const DVector &searchPoint)
{
  // use kd-Tree to search closest point
  queryPt[0] = searchPoint(0);
  queryPt[1] = searchPoint(1);
  queryPt[2] = searchPoint(2);
  kdTree->annkSearch(queryPt, 1  // number of neighbors searched
                    ,nnIdx       // nearest neighbors (returned)
                    ,dists);     // distance (returned)
  if (nnIdx[0] < 0 ) {
    // this should not happen
    cerr << endl << "ClosesSegmentSearch returned NULL-pointer" << flush;
    return NULL;
  } else {
    return idxMap[nnIdx[0]];
  }
}

///////////////////////////////////////////////////////
////////////////  Association Functions  //////////////
///////////////////////////////////////////////////////

//void associateSegmentsToTracks(Frame &frame)
//cerr << "." << flush;
//unsigned int segCount = 0;
//for (int row = 0; row < vsize; ++row) {
//  for (int col = 0; col < hsize; ++col) {
//    LidarSegment &seg = frame->segments.get(col,row);
//    if (seg.getPxCount() >= minSegmentSize)
//      unsigned int tidx = frame->projTrackIdx.get(col,row);
//      if (tidx != UINT_MAX) { // valid track was assigned to this pixel
//        LidarSegment *sp = seg.getPtr(); // get leaf pointer
//        PointCloudTrack *t = frame->trackList[tidx].get();
//        sp->linkWithTrack(t);
//        t->linkWithSegment(sp);
//      }
//    }
//  }
//}


//void LidarSegment::linkWithTrack(PointCloudTrack *t)
//{
//  LidarSegment* s = getPtr();
//  if (s->tracks.find(t) == s->tracks.end())
//    s->tracks.insert( pair<PointCloudTrack*, unsigned int>(t,1));
//  else
//    s->tracks[t]++;
//}
//
//void LidarSegment::removeLowCountLinks(double unlinkThresh)
//{
////  // remove track-links with count/maxCount < thresh (low count relative to strongest link)
//  LidarSegment* s = getPtr();
//  if (s->tracks.size() <= 1) return; // max. 1 links, so no link to remove
//  // 1) find maximum
//  unsigned int maxCount = 0;
//  typedef std::pair<PointCloudTrack*, unsigned int> pair_t;
//  BOOST_FOREACH(pair_t tp, s->tracks) {
//    maxCount = max(maxCount, tp.second);
//  }
//  // 2) delete low links
//  BOOST_AUTO(ti1, s->tracks.begin()); // points to 1st element
//  BOOST_AUTO(ti2, ti1); ++ti2; // points to 2nd element (which indeed exists)
//  while (ti1 != s->tracks.end()) {
//    if (((double)ti1->second/(double)maxCount) < unlinkThresh) {
//      ti1->first->unlinkSegment(s);
//      s->tracks.erase(ti1);
//      cout << "r" << flush;
//    }
//    ti1 = ti2;
//    ++ti2;
//  }
//}
//
//void LidarSegment::replaceTrack(PointCloudTrack *orig, PointCloudTrack *repl)
//{
//  LidarSegment* s = getPtr();
//  BOOST_AUTO(it, s->tracks.find(orig));
//  if (it != s->tracks.end()) {
//    s->tracks.insert(pair<PointCloudTrack*, unsigned int>(repl,it->second));
//    s->tracks.erase(it);
//  }
//}
//
//void LidarSegment::mergeAndCreateTracks(LidarFrame *frame, const double currEgoXVelocity, const DSMatrix &initVariance, const DSMatrix &predictionVariance)
//{
//  if (isProxy()) {
////    cerr << endl << "called LidarSegment::mergeTracks on a proxy! redirecting" << flush;
//    getPtr()->mergeAndCreateTracks(frame, currEgoXVelocity, predictionVariance, initVariance);
//    return;
//  }
//  // if no track associated -> create
//  if (tracks.empty()) {
//    DVector initState(12);
//    initState(0) = 0.0; // yaw(RAD)
//    initState(1) = 0.0; // pitch(RAD)
//    initState(2) = 0.0; // roll(RAD)
//    // Position: take coordinates of 1 point as Coordinate-System(CS) base, set orientation parallel to ego-CS
//    ColRowPair cr = *(colRows.begin()); // segment has at least one pixel
//    DVectorRange(initState,ublas::range(3,6)) = frame->point3D.get(cr); // x,y,z
//    // Speed: Assume zero ground-level-speed. As Orientation of CS is parallel to car and tracking is performed relative to ego-car, set only x-component
//    initState(6) = 0.0; // roll-rate
//    initState(7) = 0.0; // pitch-rate
//    initState(8) = 0.0; // yaw-rate
//    initState(9) = -currEgoXVelocity; // x-velocity = -ego-velocity
//    initState(10) = 0.0; // y-velocity
//    initState(11) = 0.0; // z-velocity
//    PointCloudTrack *newTrack = new PointCloudTrack(initState, initVariance, predictionVariance);
//    tracks.insert(pair<PointCloudTrack*, unsigned int>(newTrack,1)); // link segment with track
//    newTrack->linkWithSegment(this);
//    newTrack->setIndex(frame->trackList.size());
//    frame->trackList.push_back(PointCloudTrack::SPtr(newTrack));
//  } else {
//    // if several tracks associated -> merge
//    if (tracks.size() > 1) {
//      // search biggest track and select this as target (alternatively search oldest??)
//      BOOST_AUTO(tiBiggest,tracks.begin());
//      for (BOOST_AUTO(ti,tracks.begin()); ti!=tracks.end(); ++ti) {
//        if (ti->first->getPointCount()*ti->first->getAge() > tiBiggest->first->getPointCount()*tiBiggest->first->getAge())
//          tiBiggest = ti;
//      }
//      unsigned int mergedTrkCount = 0;
//      unsigned int mergedPtCount = 0;
//      unsigned int OrigPtCount = tiBiggest->first->getPointCount();
//      for (BOOST_AUTO(ti,tracks.begin()); ti!=tracks.end(); ++ti) {
//        if (ti == tiBiggest) continue; // don't merge track into itself
//        mergedPtCount += ti->first->getPointCount();
//        unsigned int tid = ti->first->getIndex();
//        tiBiggest->first->merge(ti->first); // merge track objects
//        ++mergedTrkCount;
//        // remove track from list and decrement index of all following tracks
//        frame->trackList.erase(frame->trackList.begin()+tid);
//        for (BOOST_AUTO(tit,frame->trackList.begin()+tid); tit != frame->trackList.end(); ++tit) {
//          (*tit)->setIndex((*tit)->getIndex()-1);
//        }
//      }
//      if (mergedTrkCount > 0)
//        cerr << endl << "segment[" << pxCount << "]: merged " << mergedTrkCount << " tracks with " << mergedPtCount << " points into " << OrigPtCount << " points" << flush;
//    }
//  }
//  assert(tracks.size()==1 && "LidarSegment::mergeAndCreateTracks did not result in 1 associated segment");
//}
//
//PointCloudTrack* LidarSegment::getFirstTrack()
//{
//  LidarSegment* s = getPtr();
//  if (s->tracks.empty())
//    return NULL;
//  return s->tracks.begin()->first;
//}
//
//void LidarSegment::getTracks(LidarSegment::TrackAssocIterator &begin, LidarSegment::TrackAssocIterator &end, unsigned int &count)
//{
//  LidarSegment* s = getPtr(); // this segment might by a proxy itself, so first get target pointer
//  begin = s->tracks.begin();
//  end = s->tracks.end();
//  count = s->tracks.size();
//}

//void PointCloudTrack::splitTrack(PointCloudTrack::SPtr track, LidarFrame *frame, std::vector<PointCloudTrack::SPtr> &tmpList, unsigned int &newId)
//{
//  if (segments.size() < 2) return; // max 1 segment referenced, so splitting impossible
//  // TOOODO: frequently a track is split, but one of the new tracks will not get any point assigned. check what to do!
//
//  // for each segment 2..n create new track
//  cerr << endl << "track[" << getPointCount() << "," << projPixCount <<"]: splitted into " << segments.size() << " parts" << flush;
//  map<LidarSegment*, PointCloudTrack*> newMapping; // stores connection of segments and new tracks
//  ClosestSegmentSearch segSearch(frame, segments.begin(), segments.end()); // kd-tree for splitting point cloud
//  // first segment is connected with current track, so no new track has to be created:
//  BOOST_AUTO(sit,segments.begin());
//  LidarSegment* frstSeg = (*sit)->getPtr();
//  ++sit;
//  // for all other segments, create new track object:
//  for (; sit != segments.end(); ++sit) {
//    LidarSegment* seg = (*sit)->getPtr();
//    PointCloudTrack* newTrack = new PointCloudTrack(x, P, Q);
//    newMapping.insert(pair<LidarSegment*, PointCloudTrack*>(seg,newTrack));
//    tmpList.push_back(SPtr(newTrack));
//    newTrack->setIndex(newId++);
//    newTrack->segments.insert(seg);
//    newTrack->segLastFrame.insert(segLastFrame.begin(),segLastFrame.end());
//    seg->replaceTrack(this, newTrack);
//    newTrack->age = age;
//    newTrack->lastUpdateCounter = lastUpdateCounter;
//    newTrack->expectingPrediction = expectingPrediction;
//    newTrack->projPixCount = 0; // should be calculated after splitting,
//  }
//  segments.clear(); // unlink this track with all other segments...
//  segments.insert(frstSeg); // ...except the first
//
//  // distribute track-points into new tracks
//  DMatrix R; DVector t; getRt2EgoCS(R,t);
//  DVector pt;
//  BOOST_AUTO(pit, pointCloudNGrid.begin());
//  BOOST_AUTO(pitend, pointCloudNGrid.end());
//  while (pit != pitend){
//    Surface::SPtr s = pit->s;
//    pt = prod(R,s->position) + t; // transform to Ego-CS
//    LidarSegment* seg = segSearch.findClosestSegment(s->position);
//    BOOST_AUTO(newtrkIt, newMapping.find(seg));
//    if (newtrkIt != newMapping.end()) { // associated new track found
//      BOOST_AUTO(newtrk,newtrkIt->second);
//      //newtrk->pointCloudXRelCS.splice(newtrk->pointCloudXRelCS.begin(),pointCloudXRelCS,vit);
//      newtrk->pointCloudPGrid.insert(SurfacePProxy(s));
//      newtrk->pointCloudNGrid.insert(SurfaceNProxy(s));
//      pointCloudPGrid.remove(SurfacePProxy(s));
////      pit = pointCloudNGrid.erase(pit);
//    } else {
//      // as "this" was not added to the map, this will happen, if point hits the current segment
//      // --> do nothing but increment the iterator
//      ++pit;
//    }
//  }
//  cerr << "[" << getPointCount();
//  for (BOOST_AUTO(p, newMapping.begin()); p!=newMapping.end(); ++p) {
//    cerr << "," << p->second->getPointCount();
//  }
//  cerr << "]...done" << flush;
//}




//void PointCloudTrack::linkWithSegment(LidarSegment* s)
//{
//  segments.insert(s);
//}
//
//void PointCloudTrack::unlinkSegment(LidarSegment* s)
//{
//  segments.erase(s);
//}
//
//void PointCloudTrack::mergeNonConnectedSegments()
//{
//  if (segments.size() < 2) return; // max 1 segment referenced, so merging impossible
//  // TOOODO: merge segments only if sufficient overlap with track!
//  // for each segment...
//  for (BOOST_AUTO(sit1,segments.begin()); sit1 != segments.end(); ++sit1) {
//    LidarSegment* seg1 = ((*sit1)->getPtr());
//    BOOST_AUTO(sit2,sit1); sit2++;
//    // ...check all other segments, if they can be merged
//    unsigned int mergedSegCount = 0;
//    unsigned int mergedPixCount = 0;
//    while (sit2 != segments.end()) {
//      LidarSegment* seg2 = ((*sit2)->getPtr());
//      if (seg1 == seg2) { // might be the same, if already merged by previous track that references both, too
//        BOOST_AUTO(sittmp,sit2); ++sit2;
//        segments.erase(sittmp);
//      } else {
//        if ((seg1->hasNeighbor(seg2)) || (seg2->hasNeighbor(seg1))) { // they are neighbors, i.e. have been separated by segmentation -> no merge
//          ++sit2;
//        } else {
//          mergedSegCount++;
//          mergedPixCount += seg2->getPxCount();
//          seg1->merge(seg2); // turns seg2 into proxy pointing to seg1
//          BOOST_AUTO(sittmp,sit2); ++sit2;
//          segments.erase(sittmp);
//        }
//      }
//    }
//    if (mergedSegCount > 0)
//      cerr << endl << "track[" << getPointCount() << "," << projPixCount << "]: merged " << mergedSegCount << " segments with " << mergedPixCount << " pixels into segment of " << seg1->getPxCount() << " pixels" << flush;
//  }
//}


