/*
 * Map.cpp
 *
 *  Created on: Nov 8, 2010
 *      Author: moosmann
 */

#include "Map.hpp"

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <HomogeneousTransformationMatrix.hpp>
#include <ICPChen.hpp>

using namespace std;
using namespace matrixTools;
namespace HTM = HomogeneousTransformationMatrix;
namespace fs = boost::filesystem;

Map::Map(double mapResol_, double maxDist_) :
  mapResol(mapResol_)
, cellDiag(sqrt(3*mapResol*mapResol))
, maxDist(maxDist_)
, weightZ(2.0)
, treeDim(6)
, icpNbCorresp(3)
, icpRegCoeff(0.00001)
, nbBlocks2Unwarp(30) // 10° each
, filterNCellRad(5)
, filterMinWSum(1.5) // dist=0.5->exp(-dist)=0.6
, filterMinNormalConfidence(0.9)
, filterNbNSearched(15)
, adeptNbNSearched(6)
, nbNMinFound(3)
, nMaxDist(filterNCellRad * mapResol)
, minNormalCoincidence(0.985) // cos(10°)=0,985
, pointCloudGrid(new GridT(mapResol)) // choose grid resolution equal to map resolution -> 1 point per cell
{
  nbAddedScans = 0;
  lastNbAddedScans = 0;
  dataPts = NULL;
  idxMap = NULL;
  kdTree = NULL;
  queryPt = annAllocPt(treeDim);
  resultArraySize = 0;
  nbValidResults = 0;
  nnIdx = NULL;
  dists = NULL;
  neighbors = NULL;
  weights = NULL;
  subPtsHTM = DIdMatrix(4);
  icp = NULL;

  nbAddingChangedIdx = 0;
  nbAddingDidntChangeIdx = 0;

  csi = end();
}

Map::~Map() {
  for (GridT::iterator gi=pointCloudGrid->begin(); gi!=pointCloudGrid->end(); ++gi) {
    SurfaceProxy &sp = *gi; // get single element associated with cell-container
    if ((cellDeleteNotifier))
      cellDeleteNotifier(sp.s->position[0],sp.s->position[1],sp.s->position[2],sp.s->intensity, sp.hitCount);
  }
  deinitKdTree();
}

void Map::deinitKdTree() {
  delete[] weights;
  weights = NULL;
  delete[] neighbors;
  neighbors = NULL;
  delete[] dists;
  dists = NULL;
  delete[] nnIdx;
  nnIdx = NULL;
  nbValidResults = 0;
  resultArraySize = 0;
  delete kdTree;
  kdTree = NULL;
  delete[] idxMap;
  idxMap = NULL;
  if (dataPts)
    annDeallocPts(dataPts);
  dataPts = NULL;
}

void Map::reinitializeSearch(const DVector &center, double maxDistance) {
  // TODO (5): optimize (may search directly in hash-map?) as this needs 24% of the total computation time)
  if (pointCloudGrid->size() == 0)
    return;
  // don't reinitialize if 1 frame was added and search was already reinitialized
  if ((nbAddedScans == 1) && (lastNbAddedScans == 1))
    return;
  // only reinitialize search every 2nd scan
  if (((nbAddedScans - 1 <= lastNbAddedScans)) && !(nbAddedScans == 1))
    return;
  cout << "reinitializing search..." << flush;
  deinitKdTree();

  // copy data into kd-Tree specific variables
  unsigned int maxPts = pointCloudGrid->size();
  dataPts = annAllocPts(maxPts, treeDim);// data points
  idxMap = new Surface::SPtr[maxPts];
  unsigned int nPts = 0; // actual number of data points
  BOOST_FOREACH(const SurfaceProxy &s, *pointCloudGrid) {
    if (ublas::norm_2(center - s->position) < maxDistance) {
      dataPts[nPts][0] = s->position[0];
      dataPts[nPts][1] = s->position[1];
      dataPts[nPts][2] = s->position[2]*weightZ;
      dataPts[nPts][3] = s->normal[0];//*s.normalConfidence;
      dataPts[nPts][4] = s->normal[1];//*s.normalConfidence;
      dataPts[nPts][5] = s->normal[2];//*s.normalConfidence;
      idxMap[nPts] = s.s;//&(*s);
      ++nPts;
    }
  }
  // Build kd-Tree
  kdTree = new ANNkd_tree(dataPts, nPts, treeDim);
  lastNbAddedScans = nbAddedScans;
}

Map::SearchResults Map::findClosestNeighbors(const DVector &searchPoint, const DVector &searchNormal, /*double normalConfidence,*/
    unsigned int nnCount) const {
  if (nnCount > resultArraySize) { // increase result buffer (reallocate)
    delete[] weights;
    weights = new double[nnCount];
    delete[] neighbors;
    neighbors = new Surface::SPtr[nnCount];
    delete[] dists;
    dists = new ANNdist[nnCount];
    delete[] nnIdx;
    nnIdx = new ANNidx[nnCount];
    resultArraySize = nnCount;
  }
  // use kd-Tree to search closest points
  queryPt[0] = searchPoint(0);
  queryPt[1] = searchPoint(1);
  queryPt[2] = weightZ * searchPoint(2);
  queryPt[3] = searchNormal(0);//*normalConfidence;
  queryPt[4] = searchNormal(1);//*normalConfidence;
  queryPt[5] = searchNormal(2);//*normalConfidence;
  kdTree->annkSearch(queryPt, nnCount // number of neighbors searched
      , nnIdx // nearest neighbors (returned)
      , dists); // distance (returned)
  // post-process results
  nbValidResults = 0;
  for (unsigned int idx = 0; idx < nnCount; ++idx) {
    if (nnIdx[idx] < 0)
      break;
    Surface::SPtr s = idxMap[nnIdx[idx]];
    neighbors[idx] = s;
    dists[idx] = norm_2(searchPoint - s->position); // store real euclidean distance
    assert(dists[idx] >= 0.0 && "Map::findClosestNeighbors: dists[idx] <= 0.0");
    weights[idx] = dists[0] / max(1e-8,dists[idx]);// * s->normalConfidence; // 1.0;
    weights[idx] = min(1.0, weights[idx]);
    nbValidResults++;
  }

  return SearchResults(*this);
}

void Map::registerScan(LFrameSPtr frame, LFrameSPtr lastFrame, unsigned int sampleCnt, bool unwarp) {
  cout << "." << flush;
  reg1(frame, lastFrame, sampleCnt, unwarp);
  if ((!icp) || (pointCloudGrid->empty()))
    return; // nothing to register
  icp->iterateUntilConvergence(2,10,0.01);
  reg3(frame);
}

void Map::reg1(LFrameSPtr frame, LFrameSPtr lastFrame, unsigned int sampleCnt, bool unwarp) {
  cout << endl << "registering scan..." << flush;
  if (pointCloudGrid->empty()) {
    cout << "not necessary" << flush;
    return; // nothing to register
  }

  DMatrix scanInitTransform = ublas::prod(lastFrame->positionHTM2w, frame->diffHTMToLast);
  DVector scanCenter(3); for (int i=0; i<3;++i) scanCenter[i] = scanInitTransform(i,3);
  reinitializeSearch(scanCenter, maxDist+10.0);

  // debug:
  cout << "currTD=" << frame->timeDiffToLast/1000000 << "ms..." << flush;
  cout << "icp-init: " << HTM::HTM_2_YawPitchRollXYZ(frame->positionHTM2w)
       << " + " << HTM::HTM_2_YawPitchRollXYZ(frame->diffHTMToLast) << " " << flush;

  // sub-sample:
  pts.clear();
  nrs.clear();
  nCs.clear(); //nSs.clear();
  // determine valid pixel number and sampling ratio:
  int hSize, vSize;
  frame->point3D.getSize(hSize, vSize);
  unsigned int nbPix = hSize*vSize;
  unsigned int nbPix2 = nbPix/2;
  unsigned int nbValidPixUp = 0;
  unsigned int nbValidPixLo = 0;
  unsigned int pixIdx = 0;
  BOOST_FOREACH(const double &p, frame->distance) { // iterates over rows
    if (p < maxDist) {
      if (pixIdx < nbPix2)
        ++nbValidPixUp;
      else
        ++nbValidPixLo;
    }
    ++pixIdx;
  }
  unsigned int sampleCntUp = min(3*sampleCnt/4, nbValidPixUp);
  unsigned int sampleCntLo = min(1*sampleCnt/4, nbValidPixLo);
  double sampFacUp = (double) sampleCntUp / (double)nbValidPixUp;
  double sampFacLo = (double) sampleCntLo / (double)nbValidPixLo;
  double accSamplingUp = 1.0f;
  double accSamplingLo = 1.0f;
  // determine unwarping parameters:
  unsigned int hbSize = hSize / nbBlocks2Unwarp; // might be rounded off
  double fDec = (double)hbSize / (double)hSize; // = 1/realNbBlocks
  DVector poseDiff = HTM::HTM_2_YawPitchRollXYZ(HTM::invert_HTM(frame->moveHTM));
  DMatrix Rci; DVector tci;
  if (unwarp)
    cout << "unwarping by " << poseDiff << "..." << flush;
  subPtsHTM = scanInitTransform; // points should be specified relative to end-position
  double fac = 1.0;
  for (int col = 0; col < hSize; ++col) {
    if ((unwarp) && (col % hbSize == 0)) {
//      cout << "  f=" << fac << "->" << fac*poseDiff << flush;
      HTM::YawPitchRollXYZ_2_Rt(fac*poseDiff, Rci, tci);
      fac = max(0.0,fac-fDec); // decrease influence of last registration while approaching end position
    }
    double *accSampling = &accSamplingUp;
    double *sampFac = &sampFacUp;
    for (int row = 0; row < vSize; ++row) {
      if (row == 32) {
        accSampling = &accSamplingLo;
        sampFac = &sampFacLo;
      }
      if (frame->distance.get(col, row) < maxDist) {
        *accSampling += *sampFac; // only use valid pixels
        if (*accSampling >= 1.0f) {
          *accSampling -= 1.0f;
          if (unwarp) { // unwarp using lastRegCorrection
            pts.push_back(ublas::prod(Rci,frame->point3D.get(col, row))+tci);
            nrs.push_back(ublas::prod(Rci,frame->normal3D.get(col, row)));
          } else {
            pts.push_back(frame->point3D.get(col, row));
            nrs.push_back(frame->normal3D.get(col, row));
          }
          nCs.push_back(frame->normalConfidence.get(col, row));
        //nSs.push_back(frame->normalStdDevRAD.get(col, row));
        }
      }
    }
  }

  // register:
  delete icp;
  cout << "with " << pts.size() << " samples..." << flush;
  typedef list<DVector>::iterator Pt_it_t;
  typedef list<double>::iterator D_it_t;
  BOOST_AUTO(pB, pts.begin());
  BOOST_AUTO(pE, pts.end());
  BOOST_AUTO(nB, nrs.begin());
  BOOST_AUTO(nE, nrs.end());
  BOOST_AUTO(nCB, nCs.begin());
  BOOST_AUTO(nCE, nCs.end());
  ICP::NeighborSearch<SearchResults>::ncpSearch ncpFunc = boost::bind(&Map::findClosestNeighbors, this, _1, _2, _3);
  icp = new ICP::ICPChen<Pt_it_t, SearchResults>(pB, pE, nB, nE, ncpFunc, icpNbCorresp, icpRegCoeff);
  DMatrix R; DVector t;
  HTM::HTM_2_Rt(scanInitTransform,R,t);
  if (fabs(1.0-determinant(R)) > 1e-8) {
    cout << "fixing R.." << flush;
    R = orthonormalize(R);
  }
  icp->reset(R,t);
}

void Map::reg2() {
  if ((!icp) || (pointCloudGrid->empty()))
    return; // nothing to register
  if (icp->iterate())
    icp->getEstimate(subPtsHTM); // for visualization only
  cout << "." << flush;
}

void Map::reg3(LFrameSPtr frame) {
  if ((!icp) || (pointCloudGrid->empty()))
    return; // nothing to register
  // apply:
  if (icp->iterate()) {
    DMatrix scanFinalTransform;
    icp->getEstimate(scanFinalTransform);
    //    cout << "htm: " << frame->positionHTM2w << flush;
    //    lastRegCorrection = ublas::prod(frame->positionHTM2v, scanFinalTransform);
    frame->positionHTM2w = scanFinalTransform;
    frame->positionHTM2v = HTM::invert_HTM(frame->positionHTM2w);
  }
  delete icp;
  icp = NULL;
  cout << "finished" << flush;
}

void Map::addToMap(LFrameSPtr frame, bool adept, bool unwarp, bool remOldMap, bool remFarPts) {

  // check if current position has enough advanced from last position
  if (trajectory.size() > 0) {
    const ScanPose &oldSP = trajectory.back();
    const ScanPose newSP(frame->recTime, frame->positionHTM2w);
    double moveDiff = ublas::norm_2(oldSP.getWorldPos()-newSP.getWorldPos());
    if (moveDiff < 0.05) {
      cout << endl << "skipping scan for map-integration (diff to old position is " << moveDiff << "m)" << flush;
      return;
    }
  }

  cout << endl << "adding scan..." << flush;

  // check how many pixels are valid
  unsigned int nbPixelsValid = 0;
  BOOST_FOREACH(const double &val, frame->distance) {
    if (val != DBL_MAX)
      ++nbPixelsValid;
  }
  cout << nbPixelsValid*100/frame->distance.size()/2 << "% pix ok..." << flush;
  // only clear map if frame has enough valid pixels
  if ((remOldMap) && (nbPixelsValid > 0.7*frame->distance.size()))
    pointCloudGrid->clear();
  if (pointCloudGrid->empty()) {
    adept = false;
    remFarPts = false;
  }
  // prepare
  trajectory.push_back(ScanPose(frame->recTime, frame->positionHTM2w));
  insTrajectory.push_back(ScanPose(frame->recTime, frame->insPositionHTM2w));
  DVector scanCenter(3); for (int i=0; i<3;++i) scanCenter[i] = frame->positionHTM2w(i,3);
  reinitializeSearch(scanCenter, maxDist+10.0);
  Surface::PIT poseIdx = trajectory.size()-1;
  DVector poseDiff(6);
  DMatrix currT;
  DMatrix R;
  DVector t;
  if (unwarp) {
    poseDiff = HTM::HTM_2_YawPitchRollXYZ(HTM::invert_HTM(frame->moveHTM));
    cout << "unwarping by=" << poseDiff << "..." << flush;
  } else {
    HTM::HTM_2_Rt(frame->positionHTM2w, R, t);
  }
  //  cout << " with R/t: " << R << " / " << t << "..." << flush;
  int hSize, vSize;
  frame->point3D.getSize(hSize, vSize);
  unsigned int nbPtsAdded = 0;
  unsigned int hbSize = hSize / nbBlocks2Unwarp; // might be rounded off
  double fDec = (double)hbSize / (double)hSize; // = 1/realNbBlocks
  double fac = 1.0;
  nbAddingChangedIdx = 0;
  nbAddingDidntChangeIdx = 0;
  // loop over image and add surfaces
  for (int col = 0; col < hSize; ++col) {
    if ((unwarp) && (col % hbSize == 0)) {
      // change current R/t
//      cout << "  f=" << fac << flush;
      currT = HTM::YawPitchRollXYZ_2_HTM(fac*poseDiff);
      currT = ublas::prod(frame->positionHTM2w,currT);
      HTM::HTM_2_Rt(currT, R, t);
      fac = max(0.0,fac-fDec); // decrease influence of last position while approaching end position
    }
    for (int row = 0; row < vSize; ++row) {
      DVector &pt = frame->point3D.get(col, row);
      double sDist = norm_2(pt); // distance to scanner
      if (sDist < maxDist) {
        // change to world coordinates
        DVector po = ublas::prod(R, frame->tmpVec3D.get(col, row)) + t;
        DVector p = ublas::prod(R, pt) + t;
        DVector n = ublas::prod(R, frame->normal3D.get(col, row));
        // TODO (5): shift "new" to "addToMap(s,adept)" so memory is only allocated if the surface is accepted
        Surface::SPtr s(new Surface(poseIdx, po, p, n, frame->intensity.get(col, row), frame->normalConfidence.get(col, row), sDist));
        frame->surface.set(col, row, s);
        if (addToMap(s, adept))
          ++nbPtsAdded;
      } else {
        frame->surface.set(col, row, Surface::SPtr()); // assign empty pointer
      }
    }
  }
  ++nbAddedScans;
  cout << "added " << nbPtsAdded << " points -> map-size is now " << pointCloudGrid->size() << flush;
  cout << ", adaptation changed index " << nbAddingChangedIdx*100/(nbAddingChangedIdx+nbAddingDidntChangeIdx+1) << "% of the time" << flush;

  // remove surfaces that are out-of-sight
  if (remFarPts) {
    GridT::grid_iterator  gi = pointCloudGrid->gridBegin();
    GridT::grid_iterator  gE = pointCloudGrid->gridEnd();
    while (gi != gE) {
      BOOST_AUTO(cellcenter, pointCloudGrid->getCenter(gi->first));
      double centerdist = sqrt(pow(scanCenter[0]-cellcenter[0],2)+pow(scanCenter[1]-cellcenter[1],2)+pow(scanCenter[2]-cellcenter[2],2));
      if (centerdist+cellDiag > maxDist+10.0) {
        SurfaceProxy &sp = gi->second.get(); // get single element associated with cell-container
        if ((cellDeleteNotifier))
          cellDeleteNotifier(sp.s->position[0],sp.s->position[1],sp.s->position[2],sp.s->intensity, sp.hitCount);
        gi = pointCloudGrid->erase(gi);
      } else {
        ++gi;
      }
    }
  }
}

bool Map::addToMap(Surface::SPtr s, bool adept) {
  double nDist = 0.0;
  iterator snOrig = pointCloudGrid->find_lazy_neighbor(SurfaceProxy(s), &nDist); // neighbor before adept
  if (adept) {
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
    DVector &p = s->position; // reference -> will be adepted
    DVector &n = s->normal;
    DVector q[adeptNbNSearched];
    if ((s->distFromScanner < maxDist) && (s->normalConfidence > 0.8)) { // point valid (not DBL_MAX) and normal confidence high enough
      // TODO (5): optimize: currently this needs 15% of the total calculation time
      SearchResults neighborIt = findClosestNeighbors(p, n, adeptNbNSearched);
      unsigned int nbNValid = 0;
      double aSum = 0.0; // LS-solution for adaptation
      double bSum = 0.0; // LS-solution for adaptation
      double avgNormalCoincidence = 0.0;
      for (unsigned int i = 0; i < adeptNbNSearched; ++i) { // looping over neighbors, sorted by distance
        q[i] = neighborIt.getCorrespPoint();
        DVector di = p - q[i];
        double dist = norm_2(di);
        if (dist > nMaxDist)
          break;
        const DVector &mi = neighborIt.getCorrespNormal();
        double a = ublas::inner_prod(n, mi);
        double b = ublas::inner_prod(di, mi);
        double w = a * exp(-dist);
        if ((avgNormalCoincidence + a) / (double) (nbNValid + 1) < minNormalCoincidence) // drops below threshold using this correspondence
          break;
        avgNormalCoincidence += a;
        aSum += a * w * a;
        bSum += a * w * b;
        nbNValid++;
        neighborIt.next();
      }
      if (nbNValid >= nbNMinFound) {
        avgNormalCoincidence /= (double) nbNValid;
        double c = -bSum / aSum; // LS-adaption-factor
        if (c < 0.5) {
          p += n * c;
  //        avgAccept += nbNValid;
  //        ++adepted;
        } else {
  //        ++rejectNC;
        }
      } else {
  //      ++rejectND;
      }
    }
  }
  iterator sn = pointCloudGrid->find_lazy_neighbor(SurfaceProxy(s), &nDist); // neighbor after adept

  // check whether to add or not depending on neighbor (maybe delete neighbor)
  // grid-resolution was chosen equal to desired map resolution -> 1 point per cell
  if (sn == snOrig)
    ++nbAddingDidntChangeIdx;
  else
    ++nbAddingChangedIdx;
  if ((sn != pointCloudGrid->end())) {
    // cell already contains a point, thus choose which one is better suited. criteria:
    // 1) better to be closer to the cell center
    const Surface *n = &(**sn);
    vector<double> center = pointCloudGrid->getCellCenter(sn);
    double sdcc = 0.0; // distance to center of s
    double ndcc = 0.0; // distance to center of n
    for (unsigned int i = 0; i < 3; ++i) {
      sdcc += fabs(center[i] - s->position[i]);
      ndcc += fabs(center[i] - n->position[i]);
    } // distanced are now in (0..1,5*mapResol]
    // cell-center distance difference [-1..1]: -1 = farer away, 0 = equal distance, 1 = CLOSE to center
    double ccdd = (ndcc - sdcc) / 1.5;
    // 2) better to have a smaller distance to the scanner
    // relative scanner distance (-inf..1]: <0 = farer away, 0 = equal distance, 1 = CLOSE to scanner
    double rsd = (n->distFromScanner - s->distFromScanner) / n->distFromScanner;
    // 3) better to have a higher normal confidence
    // normal confidence increase [-1..1]: -1 = lower confidence, 1 = HIGHER confidence
    double nci = s->normalConfidence - n->normalConfidence;
    // accept when scanner distance or normal confidence or cell-center-distance is reduced by 30%
    bool useNewMeasurement = (ccdd + rsd + nci > 0.3);

    // increase hit count
    SurfaceProxy &snsp = *sn;
    ++snsp.hitCount;
    // add new measurement
    if (useNewMeasurement) {
      sn->s = s;
      return true;
    }
  } else {
    // cell does not exist -> add point
    SurfaceProxy sp(s);
    pointCloudGrid->insert(sp);
    return true;
  }
  return false;
}

Map::iterator Map::modifyPos(Map::iterator si, DVector pos) {
  si->s->position = pos;
  double nDist = 0.0;
  iterator snn = pointCloudGrid->find_lazy_neighbor(SurfaceProxy(si->s), &nDist);
  if (snn != si) { // cell changed
    // just call method to add surface, will insert new cell or decide about replace existing surface
    addToMap(si->s, false); // TODO (9): could be optimized as neighbor-search is carried out twice
    // erase existing cell, return iterator to next element
    return pointCloudGrid->erase(si);
  } else {
    // if cell did not change -> nothing to do
    return ++si;
  }
}

void Map::filter() {
  filter_init();
  filter_finish();
}

void Map::filter_init() {
  // for each point:
  // if its normal is different to neighbors but they agree on the direction, remove point or adjust normal
  // check the "original" position with neighbors "original" position and adjust!

  // force "reinitializeSearch" to correctly filter map
  lastNbAddedScans = 0;
  reinitializeSearch(DZeroVector(3), DBL_MAX); // full search radius
  cout << endl << "filtering map..." << flush;
  //  list<iterator> changed;
  //  unsigned int nbChanged = 0;
  //  list<iterator> toDelete;
  //  iterator csE = end();
  //  double avgC = 0.0;

  csi = begin();
  nbDeleted = 0;
  nbChanged = 0;
  avgC = 0.0;
}

//  for (iterator csi = begin(); csi != csE; ++csi) {

void Map::filter_next() {
  //  cout << "." << flush;
  // find neighbors
  SurfaceProxy &s = *csi;
  DVector avgM = DZeroVector(3);
  double avgNc = 0.0;
  double sumW = 0.0;
  double sumWd = 0.0;
  const DVector &p = s->origPosition;
  const DVector &n = s->normal;
  SearchResults neighborIt = findClosestNeighbors(s->position, s->normal, filterNbNSearched);
  unsigned int nbNValid = 0;
  //    double avgNormalCoincidence = 0.0;
  double aSum = 0.0;
  double bSum = 0.0;
  for (unsigned int i = 0; i < filterNbNSearched; ++i) { // looping over neighbors, sorted by distance
    const DVector &qi = neighborIt.getCorrespOrigPoint(); // getCorrespPoint() ? in that case cannot overwrite "position" below
    const DVector &mi = neighborIt.getCorrespNormal();
    double nci = neighborIt.getNormConf();
    DVector di = p - qi;
    double dist = norm_2(di);
    if (dist > nMaxDist)
      break;
    if (dist > 0.0) { // skip search point as (probable) first correspondence
      double ai = ublas::inner_prod(n, mi); // =1 if normals are in same direction
      double bi = ublas::inner_prod(di, mi); // =dist if other point lies in normal direction
      double wid = exp(-dist);
      double wi = (0.5 + 0.5 * fabs(bi / dist)) * max(0.0, ai) * exp(-dist);
      sumW += wi;
      sumWd += wid;
      avgNc += wid * nci;
      avgM += wid * mi;
      //      avgNormalCoincidence += wi*ai;
      aSum += ai * wi * ai;
      bSum += ai * wi * bi;
      nbNValid++;
    }
    if (!neighborIt.next())
      break;
  }
  if ((sumW >= 0.1) && (sumWd >= filterMinWSum) && (nbNValid >= nbNMinFound)) {
    // if enough neighbors found, adept point
    avgNc /= sumWd;
    double avgMLength = norm_2(avgM);
    if (avgMLength > 0.0) // avoid div by zero
      avgM /= norm_2(avgM); // norm to length 1
    if ((avgNc > filterMinNormalConfidence) && (ublas::inner_prod(avgM, n) < minNormalCoincidence)) {
      // remove point if its normal is different to neighbors normals and normal confidence is high
      csi = pointCloudGrid->erase(csi);
      ++nbDeleted;
    } else {
      // adapt point based on "original point" coordinates
      double c = -bSum / aSum; // LS-adaption-factor
      avgC += c;
      csi = modifyPos(csi, p + c * n);
      ++nbChanged;
      //        cout << endl << "adapted with " << nbNValid << " neighbors by " << c << flush;
    }
  } else {
    // not enough neighbors -> continue with next point
    ++csi;
  }
}
//  }

void Map::filter_finish() {
  while (csi != end())
    filter_next();
  if (nbDeleted > 0)
    cout << "deleted " << nbDeleted << " points..." << flush;
  if (nbChanged > 0)
    cout << "changed " << nbChanged << " points in avg by " << avgC / (double) nbChanged << "..." << flush;

  // force "reinitializeSearch", as some elements don't exist anymore / were moved
  lastNbAddedScans = 0;
  cout << "done" << flush;
}

void Map::saveToFile(std::string filename) const {
  cout << "trying to save map to file \"" << filename << "\"..." << flush;
  std::ofstream ofs(filename.c_str());
  assert(ofs.good() && "Map::saveToFile: problem opening the file for writing");
  string ext = extension(fs::path(filename));
  if (ext == "txt") {
    boost::archive::text_oarchive  archive(ofs);
    archive << boost::serialization::make_nvp("Map", *this);
  } else if (ext == "xml") {
    boost::archive::xml_oarchive  archive(ofs);
    archive << boost::serialization::make_nvp("Map", *this);
  } else {
    boost::archive::binary_oarchive  archive(ofs);
    archive << boost::serialization::make_nvp("Map", *this);
  }
  cout << "done" << flush;
}

void Map::loadFromFile(std::string filename) {
  cout << "trying to load map from file \"" << filename << "\"..." << flush;
  std::ifstream ifs(filename.c_str());
  assert(ifs.good() && "Map::loadFromFile: problem opening the file for reading");
  
  pointCloudGrid->clear(); // clear all scans stored until now
  string ext = extension(fs::path(filename));
  if (ext == "txt") {
    boost::archive::text_iarchive  archive(ifs);
    archive >> boost::serialization::make_nvp("Map", *this);
  } else if (ext == "xml") {
    boost::archive::xml_iarchive  archive(ifs);
    archive >> boost::serialization::make_nvp("Map", *this);
  } else {
    boost::archive::binary_iarchive  archive(ifs);
    archive >> boost::serialization::make_nvp("Map", *this);
  }
  cout << "done" << flush;
  csi = end();
}

void exportTrajectory(const Map::TrajT &traj, string filename)
{
  ofstream out(filename.c_str());
  BOOST_FOREACH(const ScanPose &pose, traj) {
    out << pose.getTimestamp() << " ";
    const matrixTools::DMatrix& htm = pose.getHTM2w();
    for (unsigned int row=0; row<4;++row)
      for (unsigned int col=0; col<4;++col)
        out << htm(row,col) << " ";
    out << endl;
  }
}

void importTrajectory(Map::TrajT &traj, string filename)
{
  ifstream in(filename.c_str());
  if (!in.good())
    throw runtime_error("Map::importTrajectory: problem opening the file for reading");
  traj.clear();
  string line;
  int64_t timestamp;
  matrixTools::DMatrix htm2w(4,4);
  while (getline(in, line)) { // read one line into string
    istringstream istr;   // streaming the string containing a line of the calib_file
    istr.str(line);
    istr >> timestamp;
    for (unsigned int row=0; row<4;++row)
      for (unsigned int col=0; col<4;++col)
        istr >> htm2w(row,col);
    traj.push_back(ScanPose(timestamp, htm2w));
  }
}

void Map::exportTrajectories(string basefilename)
{
  exportTrajectory(trajectory, basefilename + ".traj");
  exportTrajectory(insTrajectory, basefilename + ".ins.traj");
}

void Map::importTrajectories(string basefilename)
{
  importTrajectory(trajectory, basefilename + ".traj");
  importTrajectory(insTrajectory, basefilename + ".ins.traj");
}

void Map::exportMap(string filename, unsigned int resolutionFactor, ExportType type)
{
  //save trajectory
  if (type == Surfaces) {
    exportTrajectory(trajectory, filename + ".traj");
    exportTrajectory(insTrajectory, filename + ".ins.traj");
  }

  //save map
  ofstream out(filename.c_str());
  Map::const_iterator mapIt = pointCloudGrid->begin();
  Map::const_iterator mapEnd = pointCloudGrid->end();
  unsigned int nbTotal = 0;
  unsigned int nbExported = 0;
  while (mapIt != mapEnd) {
    const Surface s = **mapIt;
    bool skip = false;
    if (resolutionFactor > 1) { // check whether to store surface
      for (unsigned int d=0; (d<3) && (!skip); ++d) {
        int ii = (int)(s.position[d]/mapResol);
        if (s.position[d] < 0) ii-=1;  // correct neg. indices due to rounding
        if (ii%resolutionFactor != 0)
          skip = true;
      }
    }
    ++nbTotal;
    if (!skip) {
      ++nbExported;
      // write information into file:
      if (type == Surfaces)
        out << s.poseIdx << " ";
      out << s.position[0] << " ";
      out << s.position[1] << " ";
      out << s.position[2] << " ";
      if (type == PointIntens)
        out << s.intensity << " ";
      if (type == Surfaces) {
        out << s.origPosition[0] << " ";
        out << s.origPosition[1] << " ";
        out << s.origPosition[2] << " ";
        out << s.normal[0] << " ";
        out << s.normal[1] << " ";
        out << s.normal[2] << " ";
        out << s.normalConfidence << " ";
        out << s.distFromScanner << " ";
        out << mapIt->hitCount << " ";
        out << s.intensity << " ";
        out << endl;
      }
    }
    ++mapIt;
  }
  cout << "exported " << nbExported << " of " << nbTotal << " (" << nbExported*100/nbTotal << "%)" << flush;
}

void Map::importMap(string filename)
{
  //load trajectory
  importTrajectory(trajectory, filename + ".traj");
  importTrajectory(insTrajectory, filename + ".ins.traj");
  
  //load map
  ifstream in(filename.c_str());
  if (!in.good())
    throw runtime_error("Map::importMap: problem opening the file for reading");
  pointCloudGrid->clear(); // clear all scans stored until now
  string line;
  Surface tmp;
  unsigned int hitCnt = 1;
  while (getline(in, line)) { // read one line into string
    istringstream istr;   // streaming the string containing a line of the calib_file
    istr.str(line);
    istr >> tmp.poseIdx;
    istr >> tmp.position[0];
    istr >> tmp.position[1];
    istr >> tmp.position[2];
    istr >> tmp.origPosition[0];
    istr >> tmp.origPosition[1];
    istr >> tmp.origPosition[2];
    istr >> tmp.normal[0];
    istr >> tmp.normal[1];
    istr >> tmp.normal[2];
    istr >> tmp.normalConfidence;
    istr >> tmp.distFromScanner;
    if (!istr.eof())
      istr >> hitCnt;
    if (!istr.eof())
      istr >> tmp.intensity;
    pointCloudGrid->insert(SurfaceProxy(Surface::SPtr(new Surface(tmp)),hitCnt));
  }
  csi = end();
}

void Map::importPointCloud(string filename, bool withIntensity)
{
  //load trajectory
  trajectory.clear();
  insTrajectory.clear();

  //load map
  ifstream in(filename.c_str());
  if (!in.good())
    throw runtime_error("Map::importPointCloud: problem opening the file for reading");
  pointCloudGrid->clear(); // clear all scans stored until now
  Surface tmp;
  tmp.poseIdx = 0;
  tmp.normal = DZeroVector(3);
  tmp.normalConfidence = 0.0;
  tmp.distFromScanner = 0.0;
  unsigned int hitCnt = 1;
  unsigned int i=0;
  cout << "reading..." << flush;
  while (!in.eof()) {
    if ((++i)%3000 == 0) cout << "." << flush;
    in >> tmp.position[0];
    in >> tmp.position[1];
    in >> tmp.position[2];
    if (withIntensity)
      in >> tmp.intensity;
    tmp.origPosition = tmp.position;
    pointCloudGrid->insert(SurfaceProxy(Surface::SPtr(new Surface(tmp)),hitCnt));
  }
  cout << "done" << flush;
  csi = end();
}
