#include "AlgoRegistration.hpp"

#include <cfloat>
#include <cmath>
#include <set>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/filesystem.hpp>
#include "HomogeneousTransformationMatrix.hpp"
#include "ICPBesl.hpp" // ICP implementation
#include "ICPChen.hpp" // ICP implementation
#include "ICPLinearized.hpp" // ICP implementation
#include "ICPSwitch.hpp" // ICP implementation
#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;
using namespace boost::numeric::ublas;
namespace HTM = HomogeneousTransformationMatrix;

SearchResultsBase::SearchResultsBase(SearchResultsBase* sr_) : sr(sr_) {};
SearchResultsBase::~SearchResultsBase() {delete sr;};
bool SearchResultsBase::next() {return sr->next();};
mdefs::DVector SearchResultsBase::getCorrespPoint() {return sr->getCorrespPoint();};
mdefs::DVector SearchResultsBase::getCorrespNormal() {return sr->getCorrespNormal();};
double SearchResultsBase::getNormConf() {return sr->getNormConf();};
double SearchResultsBase::getNormStdDevRAD() {return sr->getNormStdDevRAD();};
mdefs::DSMatrix SearchResultsBase::getPCovar() {return sr->getPCovar();};
mdefs::DSMatrix SearchResultsBase::getNCovar() {return sr->getNCovar();};
double SearchResultsBase::getDist() {return sr->getDist();};
double SearchResultsBase::getWeight() {return sr->getWeight();};

unsigned int debugBreak = 70;

ProjectedNeighborSearch::ProjectedNeighborSearch(FrameSPtr frame_, const LidarImageProjector *projector_)
  : frame(frame_)
   ,projector(projector_)
{
  ParameterHeap* params = ParameterHeap::get();
  maxDist = -params->regCorrespDistWithHalfWeight / log(0.5); // -> weight=0.5 at regCorrespDistWithHalfWeight
  weightZ = params->matchWeightZ;
  projDistTolerance = params->projectionDistTolerance;
  useNormalWeight = params->regIcpUseNormalWeight;
  useDistDiffWeight = params->regIcpUseDistDiffWeight;
  useOcclWeight = params->regIcpUseOcclWeight;
  usePixMask = true;
  occlWeight = params->regIcpOcclWeight;
  if (useNormalWeight) cout << " using NormalWeight " << flush;
  if (useDistDiffWeight) cout << " using useDistDiffWeight " << flush;
  if (useOcclWeight) cout << " using useOcclWeight (" << occlWeight << ") " << flush;

  // copy 3D data into kd-Tree specific variables
  int hsize, vsize;
  frame->point3D.getSize(hsize,vsize);
  const int nDim = 6; // dimensionality of Search space, x,y,z
  //int maxPts = hsize*vsize; // maximum number of points, used for memory allocation
  int nPts = 0; // actual number of data points
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      if (frame->point3D.get(col,row)(0) != DBL_MAX)
        ++nPts;
    }
  }
  dataPts.resize(boost::extents[nPts][nDim]);
  idxMap = new int[nPts*2];
  int idx = 0;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      DVector &p = frame->point3D.get(col,row);
      if ((p(0) != DBL_MAX)) {
        dataPts[idx][0] = p(0);
        dataPts[idx][1] = p(1);
        dataPts[idx][2] = p(2);
        DVector &n = frame->normal3D.get(col,row);
        double nMult = 2.0 * max(0.1,frame->normalConfidence.get(col,row)); /* change here -> change below */
        assert((n.size()==3) && "ProjectedNeighborSearch::findClosestNeighbor: normal has size != 3");
        assert(((n(0) != DBL_MAX)&&(n(1) != DBL_MAX)&&(n(2) != DBL_MAX)) && "ProjectedNeighborSearch: normal is DBL_MAX");
        dataPts[idx][3] = n(0)*nMult;
        dataPts[idx][4] = n(1)*nMult;
        dataPts[idx][5] = n(2)*nMult;
        idxMap[2*idx + 0] = col;
        idxMap[2*idx + 1] = row;
        ++idx;
      }
    }
  }
  if (idx != nPts) throw std::logic_error("ProjectedNeighborSearch: idx != nPts");
  dataPtsStretch.resize(boost::extents[nPts][nDim]);
  for (int pt = 0; pt < nPts; ++pt) {
    for (int d = 0; d < nDim; ++d)
      dataPtsStretch[pt][d] = dataPts[pt][d];
    dataPtsStretch[pt][2] *= weightZ;
  }
  // Build kd-Tree
  kdTree6D = new kdtree2::KDTree(dataPts,true,6); // KDTreeArray& data_in, bool rearrange_in = true, int only_use_ndim_in=-1
  kdTree6D->sort_results = true;
  kdTree3D = new kdtree2::KDTree(dataPts,true,3); // KDTreeArray& data_in, bool rearrange_in = true, int only_use_ndim_in=-1
  kdTree3D->sort_results = true;
  kdTreeStretch = new kdtree2::KDTree(dataPtsStretch,true,6); // KDTreeArray& data_in, bool rearrange_in = true, int only_use_ndim_in=-1
  kdTreeStretch->sort_results = true;
  kdTree = kdTree6D;
  // Allocate Variables for Queries
  queryPt.resize(nDim);         // query point
}

ProjectedNeighborSearch::~ProjectedNeighborSearch()
{
  delete kdTreeStretch;
  delete kdTree3D;
  delete kdTree6D;
  delete [] idxMap;
}


bool ProjectedNeighborSearch::SearchResults::next() {
  ++index;
  correspValid = false;
  if (index >= ns.result.size()) return false;
  return true;
}

void ProjectedNeighborSearch::SearchResults::assertCorrespValid() {
  if (index >= ns.result.size()) {
    throw range_error("ProjectedNeighborSearch::SearchResults: access on search result > resultSize");
  }
  if (!correspValid) { // only do calculation if not calculated before
    // 1) find closest point-neighbor
    int colNN, rowNN;
    ns.getColRow(index, colNN, rowNN); // might throw logic_error on pixMask==false
    correspPoint = ns.frame->point3D.get(colNN,rowNN);
    correspNormal = ns.frame->normal3D.get(colNN,rowNN);
    correspNormConf = ns.frame->normalConfidence.get(colNN,rowNN);
    correspNormStdDev = ns.frame->normalStdDevRAD.get(colNN,rowNN);
    correspPCovar = ns.frame->pointVariance3D.get(colNN,rowNN);
    correspNCovar = ns.frame->normalVariance3D.get(colNN,rowNN);
    ns.result[index].dis = norm_2(ns.searchPoint-correspPoint); // store real distance
    //ns.dists[index] = sqrt(ns.dists[index]); // unsquare distance, not possible since normal is contained in feature space
    if ((pmode == PROJECT) && (norm_1(correspNormal) > 0.1)) { // normal is valid -> project search point to plane of correspondence, otherwise just return correspondence
      double projRelTarget = inner_prod(correspNormal, correspPoint-ns.searchPoint); // inner_prod(n,v)=|n|*|v|*cos(angle(n,v))=|n|*|v|*(|vp|/|v|)=|vp|  ,since |n|=1
      DVector relCorresp = projRelTarget * correspNormal; // relative to searchPoint
      correspPoint = relCorresp + ns.searchPoint; // turn into absolute coordinates
    }
    assert(!std::isnan(correspPoint(0)) && "ProjectedNeighborSearch::assertCorrespValid: correspPoint is NAN");
    assert(!std::isnan(correspPoint(1)) && "ProjectedNeighborSearch::assertCorrespValid: correspPoint is NAN");
    assert(!std::isnan(correspPoint(2)) && "ProjectedNeighborSearch::assertCorrespValid: correspPoint is NAN");
    assert(!std::isnan(correspNormal(0)) && "ProjectedNeighborSearch::assertCorrespValid: correspNormal is NAN");
    assert(!std::isnan(correspNormal(1)) && "ProjectedNeighborSearch::assertCorrespValid: correspNormal is NAN");
    assert(!std::isnan(correspNormal(2)) && "ProjectedNeighborSearch::assertCorrespValid: correspNormal is NAN");
    ns.correspondenceBuffer.push_back(ns.searchPoint);
    ns.correspondenceBuffer.push_back(correspPoint);
    correspValid = true;
  }
}

mdefs::DVector ProjectedNeighborSearch::SearchResults::getCorrespPoint() {
  assertCorrespValid();
  return correspPoint;
}

mdefs::DVector ProjectedNeighborSearch::SearchResults::getCorrespNormal() {
  assertCorrespValid();
  return correspNormal;
}

double ProjectedNeighborSearch::SearchResults::getNormConf() {
  assertCorrespValid();
  return correspNormConf;
}

double ProjectedNeighborSearch::SearchResults::getNormStdDevRAD() {
  assertCorrespValid();
  return correspNormStdDev;
}

mdefs::DSMatrix ProjectedNeighborSearch::SearchResults::getPCovar() {
  assertCorrespValid();
  return correspPCovar;
}

mdefs::DSMatrix ProjectedNeighborSearch::SearchResults::getNCovar() {
  assertCorrespValid();
  return correspNCovar;
}

double ProjectedNeighborSearch::SearchResults::getDist() {
  assertCorrespValid();
  assert(!std::isnan(ns.result[index].dis) && "ProjectedNeighborSearch::findClosestNeighbor: returned distance is NAN");
  return ns.result[index].dis;
}

double ProjectedNeighborSearch::SearchResults::getWeight() {
  assertCorrespValid();
  if (ns.result[index].dis > ns.maxDist) // 3D euclidean distance
    return 0.0;
  double weight = 1.0;
  weight *= exp(-(ns.result[index].dis/ns.maxDist)); // always weight by absolute distance
  //^^for maxDist=10m (dist->weight): 10m->=0.36, 2m->=0.82, 1m->=0.90, 0.5m->=0.95, 0.2m->=0.98, 0.1m->=0.99
  if (index > 0) { // weight by relative distance (closest correspondence gets w=1.0, other neighbors < 1.0)
    if (ns.result[index].dis > 0.0)
      weight *= ns.result[0].dis/ns.result[index].dis;
  }
  if (ns.useNormalWeight) { // penalizes correspondences with different normal direction and correspondences along a surface (good for ICP convergence)
    weight *= 0.05 + 0.95*max(0.0,inner_prod(correspNormal, ns.searchNormal)); // parallelness of normal vectors, if opposite direction (inner_prod < 0) then set to 0.0
    //weight *= abs(inner_prod(searchNormal, (searchPoint-correspPoint))); // perpendicularity of normal-correspondence
  }
  if (ns.useOcclWeight) { // penalizes if searchPoint is not visible (if activated)
    if (ns.pixPlaneDist < -ns.projDistTolerance)
      weight *= ns.occlWeight;
  }
  if (ns.useDistDiffWeight) { // penalizes if inter-frame-pixel-distance remains unchanged
    int correspCol, correspRow;
    ns.getColRow(index, correspCol, correspRow); //!< access the kd-tree query result at the specified index and convert into col/row
    double distDiff = ns.frame->distanceDiffLastFrame.get(correspCol,correspRow);
    if (distDiff != DBL_MAX)
      weight *= (1-exp(-fabs(distDiff)));
  }
  assert(!std::isnan(weight) && "ProjectedNeighborSearch::setWeight: weight is NAN");
  return weight;
}

SearchResultsBase ProjectedNeighborSearch::findClosestNeighbor(const DVector &searchPoint, const DVector &searchNormal, double nConf, unsigned int nnCount) const
{
  findNeighbors(searchPoint, searchNormal, nConf, nnCount); //!< calls the kd-tree search function, storing results in the kd-tree query arrays
  return SearchResultsBase(new SearchResults(*this, SearchResults::DONTPROJECT));
}

SearchResultsBase ProjectedNeighborSearch::findClosestProjectedNeighbor(const DVector &searchPoint, const DVector &searchNormal, double nConf, unsigned int nnCount) const
{
  findNeighbors(searchPoint, searchNormal, nConf, nnCount); //!< calls the kd-tree search function, storing results in the kd-tree query arrays
  return SearchResultsBase(new SearchResults(*this, SearchResults::PROJECT));
}

void ProjectedNeighborSearch::findNeighbors(const DVector &searchPt, const DVector &searchNrm, double nConf, unsigned int nnCount) const
{
  assert((searchPt.size()==3) && "ProjectedNeighborSearch::findClosestNeighbor: point has size != 3");
  assert((searchNrm.size()==3) && "ProjectedNeighborSearch::findClosestNeighbor: normal has size != 3");
  assert((searchPt(0)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
  assert((searchPt(1)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
  assert((searchPt(2)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
  assert((searchNrm(0)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
  assert((searchNrm(1)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
  assert((searchNrm(2)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
  assert(!std::isnan(searchPt(0)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is NAN");
  assert(!std::isnan(searchPt(1)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is NAN");
  assert(!std::isnan(searchPt(2)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPt is NAN");
  assert(!std::isnan(searchNrm(0)) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is NAN");
  assert(!std::isnan(searchNrm(1)) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is NAN");
  assert(!std::isnan(searchNrm(2)) && "ProjectedNeighborSearch::findClosestNeighbor: searchNrm is NAN");
  searchPoint = searchPt;
  searchNormal = searchNrm;
  // use kd-Tree to search closest point
  queryPt[0] = searchPoint(0);
  queryPt[1] = searchPoint(1);
  queryPt[2] = searchPoint(2);
  if (kdTree == kdTreeStretch)
    queryPt[2] *= weightZ;
  double nMult = 2.0 * max(0.1,nConf); /* change here -> change above */
  queryPt[3] = searchNormal(0)*nMult;
  queryPt[4] = searchNormal(1)*nMult;
  queryPt[5] = searchNormal(2)*nMult;
  kdTree->n_nearest(queryPt,nnCount,result); // search for 10 of them.
  //result[k].dis // its square Euclidean distance
  //result[k].idx // index to the data
  pixPlaneDist = -DBL_MAX;
  if (projector) {
    int colPC, rowPC; double distToCam; // project search point to image pixel
    if (projector->getImageIndexRel(searchPoint(0), searchPoint(1), searchPoint(2), colPC, rowPC, distToCam)) {
      double pixDist = frame->distance.get(colPC,rowPC);
      if (pixDist != DBL_MAX) {
        if (nConf > 0.7) {
          DVector& pixPt = frame->point3D.get(colPC,rowPC);
          DVector& pixNr = frame->normal3D.get(colPC,rowPC);
          pixPlaneDist = inner_prod(searchPoint-pixPt, pixNr);
        } else {
          pixPlaneDist = pixDist - distToCam;
        }
      }
    }
  } else {
    pixPlaneDist = 0.0;
  }
}

void ProjectedNeighborSearch::getColRow(unsigned int index, int &col, int &row) const {
  assert((index < result.size()) && "ProjectedNeighborSearch::getColRow: index out of range");
  if (result[index].idx < 0 )
    throw logic_error("ProjectedNeighborSearch::findClosestNeighbor: invalid correspondence of kd-Tree");
  col = idxMap[2*result[index].idx+0];
  row = idxMap[2*result[index].idx+1];
  if ((usePixMask) && (frame->regAvailPixels.get(col,row) == false)) // check if neighbor is valid and still available for use
    throw logic_error("ProjectedNeighborSearch::findClosestNeighbor: invalid correspondence: pixel was unmasked");
}




TrackNeighborSearch::TrackNeighborSearch(PointCloudTrack::SPtr track)
{
  ParameterHeap* params = ParameterHeap::get();
  maxDist = -params->regCorrespDistWithHalfWeight / log(0.5); // -> weight=0.5 at regCorrespDistWithHalfWeight

  // copy 3D data into kd-Tree specific variables
  size_t nPts = track->getPointCount();
  if (nPts > 0) {
    const int nDim = 6; // dimensionality of Search space, x,y,z
    PointCloudTrack::PointConstIterator ptIt, ptEnd;
    track->getPoints(ptIt, ptEnd);
    surfaces.resize(nPts);
    dataPts.resize(boost::extents[nPts][nDim]);
    size_t idx = 0;
    while (ptIt != ptEnd) {
      surfaces[idx] = ptIt->s.get();
      dataPts[idx][0] = surfaces[idx]->position(0);
      dataPts[idx][1] = surfaces[idx]->position(1);
      dataPts[idx][2] = surfaces[idx]->position(2);
      double nMult = 2.0 * max(0.1,surfaces[idx]->normalConfidence); /* change here -> change below */
      dataPts[idx][3] = surfaces[idx]->normal(0)*nMult;
      dataPts[idx][4] = surfaces[idx]->normal(1)*nMult;
      dataPts[idx][5] = surfaces[idx]->normal(2)*nMult;
      ++idx;
      ++ptIt;
    }
    if (idx != nPts) throw std::logic_error("TrackNeighborSearch: idx != nPts");
    // Build kd-Tree
    kdTree = new kdtree2::KDTree(dataPts,true,nDim); // KDTreeArray& data_in, bool rearrange_in = true, int only_use_ndim_in=-1
    kdTree->sort_results = true;
    // Allocate Variables for Queries
    queryPt.resize(nDim);         // query point
  } else {
    kdTree = NULL;
  }
}

TrackNeighborSearch::~TrackNeighborSearch()
{
  delete kdTree;
}

SearchResultsBase TrackNeighborSearch::findNeighbors(const DVector &searchPt, const DVector &searchNrm, double nConf, unsigned int nnCount) const
{
  if (kdTree) {
    assert((searchPt.size()==3) && "TrackNeighborSearch::findClosestNeighbor: point has size != 3");
    assert((searchNrm.size()==3) && "TrackNeighborSearch::findClosestNeighbor: normal has size != 3");
    assert((searchPt(0)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
    assert((searchPt(1)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
    assert((searchPt(2)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchPt is DBL_MAX");
    assert((searchNrm(0)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
    assert((searchNrm(1)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
    assert((searchNrm(2)!=DBL_MAX) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is DBL_MAX");
    assert(!std::isnan(searchPt(0)) && "TrackNeighborSearch::findClosestNeighbor: searchPt is NAN");
    assert(!std::isnan(searchPt(1)) && "TrackNeighborSearch::findClosestNeighbor: searchPt is NAN");
    assert(!std::isnan(searchPt(2)) && "TrackNeighborSearch::findClosestNeighbor: searchPt is NAN");
    assert(!std::isnan(searchNrm(0)) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is NAN");
    assert(!std::isnan(searchNrm(1)) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is NAN");
    assert(!std::isnan(searchNrm(2)) && "TrackNeighborSearch::findClosestNeighbor: searchNrm is NAN");
    searchPoint = searchPt;
    searchNormal = searchNrm;
    // use kd-Tree to search closest point
    queryPt[0] = searchPoint(0);
    queryPt[1] = searchPoint(1);
    queryPt[2] = searchPoint(2);
    double nMult = 2.0 * max(0.1,nConf); /* change here -> change above */
    queryPt[3] = searchNormal(0)*nMult;
    queryPt[4] = searchNormal(1)*nMult;
    queryPt[5] = searchNormal(2)*nMult;
    kdTree->n_nearest(queryPt,nnCount,result); // search for 10 of them.
    //result[k].dis // its square Euclidean distance
    //result[k].idx // currResIndex to the data
    for (unsigned int i=0; i<result.size(); ++i)
      result[i].dis = norm_2(searchPoint - surfaces[result[i].idx]->position); // store real distance
  } else { // might not be valid if established on 0 points
    result.resize(0);
  }
  return SearchResultsBase(new SearchResults(*this));
}

bool TrackNeighborSearch::SearchResults::next() {
  ++currResIndex;
  if (currResIndex >= ns.result.size()) return false;
  return true;
}

mdefs::DVector TrackNeighborSearch::SearchResults::getCorrespPoint() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->position;
}

mdefs::DVector TrackNeighborSearch::SearchResults::getCorrespNormal() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->normal;
}

double TrackNeighborSearch::SearchResults::getNormConf() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->normalConfidence;
}

double TrackNeighborSearch::SearchResults::getNormStdDevRAD() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->normalStdDevRAD;
}

mdefs::DSMatrix TrackNeighborSearch::SearchResults::getPCovar() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->ptCovar;
}

mdefs::DSMatrix TrackNeighborSearch::SearchResults::getNCovar() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  return ns.surfaces[ns.result[currResIndex].idx]->nrCovar;
}

double TrackNeighborSearch::SearchResults::getDist() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  assert(!std::isnan(ns.result[currResIndex].dis) && "TrackNeighborSearch::findClosestNeighbor: returned distance is NAN");
  return ns.result[currResIndex].dis;
}

double TrackNeighborSearch::SearchResults::getWeight() {
  if (currResIndex >= ns.result.size())
    throw range_error("TrackNeighborSearch::SearchResults: access on search result > resultSize");
  if (ns.result[currResIndex].dis > ns.maxDist) // 3D euclidean distance
    return 0.0;
  double weight = 1.0;
  weight *= exp(-(ns.result[currResIndex].dis/ns.maxDist)); // always weight by absolute distance
  //^^for maxDist=10m (dist->weight): 10m->=0.36, 2m->=0.82, 1m->=0.90, 0.5m->=0.95, 0.2m->=0.98, 0.1m->=0.99
  if (currResIndex > 0) { // weight by relative distance (closest correspondence gets w=1.0, other neighbors < 1.0)
    if (ns.result[currResIndex].dis > 0.0)
      weight *= ns.result[0].dis/ns.result[currResIndex].dis;
  }
  assert(!std::isnan(weight) && "TrackNeighborSearch::setWeight: weight is NAN");
  return weight;
}




FeatureMatching::FeatureMatching(FrameSPtr currFrame, FrameSPtr lastFrame)
  : cf(currFrame)
   ,lf(lastFrame)
{
  ParameterHeap* params = ParameterHeap::get();
  featFac = params->matchWeightFeat;
  vertFac = params->matchWeightZ;
  maxDistanceSqr = pow(params->matchMaxFSpaceDist,2);
  // copy data into kd-Tree specific variables
  int hsize, vsize;
  cf->point3D.getSize(hsize,vsize);
  int nDim = 3+4;                                   // dimensionality of Search space, x,y,z+Features
  int maxPts = hsize*vsize;                         // maximum nuber of points, used for memory allocation
  dataPts = annAllocPts(maxPts, nDim);// data points
  idxMap = new int[maxPts*2];
  int nPts = 0; // actual number of data points
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      DVector &vf = cf->matchingFeatures.get(col,row);
      DVector &pt = cf->point3D.get(col,row);
      if ((pt(0) != DBL_MAX) && (vf.size() == 4)) {
        dataPts[nPts][0] = pt(0);
        dataPts[nPts][1] = pt(1);
        dataPts[nPts][2] = pt(2)*vertFac;
        dataPts[nPts][3] = featFac*vf(0);
        dataPts[nPts][4] = featFac*vf(1);
        dataPts[nPts][5] = featFac*vf(2);
        dataPts[nPts][6] = featFac*vf(3);
        idxMap[2*nPts+0] = col;
        idxMap[2*nPts+1] = row;
        ++nPts;
      }
    }
  }
  // Build kd-Tree
  kdTree = new ANNkd_tree(dataPts, nPts, nDim);
  // Allocate Variables for Queries
  queryPt = annAllocPt(nDim);         // query point
  nnIdx = new ANNidx[1];              // nearest neighbor index for 1 result
  dists = new ANNdist[1];             // nearest neighbor dists for 1 result
}
FeatureMatching::~FeatureMatching()
{
  delete [] dists;
  delete [] nnIdx;
  delete kdTree;
  delete [] idxMap;
//  annClose(); // done with ANN
}
void FeatureMatching::findClosestNeighbor(int colL, int rowL, int &colC, int &rowC, double &dist) const
{
  DVector &vf = lf->matchingFeatures.get(colL,rowL);
  DVector &pt = lf->point3D.get(colL,rowL);
  if ((pt(0)==DBL_MAX) || (pt(1)==DBL_MAX) || (pt(2)==DBL_MAX)) // this will happen if point is invalid
    throw logic_error("ProjectedNeighborSearch::findClosestNeighbor: searchPoint is DBL_MAX");
  assert(!std::isnan(pt(0)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPoint is NAN"); // this should never happen
  assert(!std::isnan(pt(1)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPoint is NAN");
  assert(!std::isnan(pt(2)) && "ProjectedNeighborSearch::findClosestNeighbor: searchPoint is NAN");
  assert((vf.size()==4) && "ProjectedNeighborSearch::findClosestNeighbor: size(featureVector) != 4"); // this will happen if features were not calculated
  assert((vf(0)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is DBL_MAX"); // this should never happen
  assert((vf(1)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is DBL_MAX");
  assert((vf(2)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is DBL_MAX");
  assert((vf(3)!=DBL_MAX) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is DBL_MAX");
  assert(!std::isnan(vf(0)) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is NAN"); // this should never happen
  assert(!std::isnan(vf(1)) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is NAN");
  assert(!std::isnan(vf(2)) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is NAN");
  assert(!std::isnan(vf(3)) && "ProjectedNeighborSearch::findClosestNeighbor: featureVector is NAN");
  queryPt[0] = pt(0);
  queryPt[1] = pt(1);
  queryPt[2] = pt(2)*vertFac;
  queryPt[3] = featFac*vf(0);
  queryPt[4] = featFac*vf(1);
  queryPt[5] = featFac*vf(2);
  queryPt[6] = featFac*vf(3);
  kdTree->annkSearch(queryPt, 1  // number of neighbors searched
                    ,nnIdx       // nearest neighbors (returned)
                    ,dists);     // distance (returned)
  if (nnIdx[0] < 0)
    throw logic_error("ProjectedNeighborSearch::findClosestNeighbor: invalid correspondence of kd-Tree");
  if (dists[0] > maxDistanceSqr)
    throw logic_error("ProjectedNeighborSearch::findClosestNeighbor: correspondence exceeds maxDistance");
  dist = sqrt(dists[0]);       // unsquare distance
  colC = idxMap[2*nnIdx[0]+0];
  rowC = idxMap[2*nnIdx[0]+1];
}


TrackRegistration::TrackRegistration(FrameSPtr cf_, FrameSPtr lf_, const LidarImageProjector &projector_)
  : cf(cf_)
  , lf(lf_)
  , projector(projector_)
  , pnSearch(cf,&projector)
  , fMatcher(NULL)
  , enPtPl()
  , enPtPt()
  , enComb(enPtPl, enPtPt)
  , matchFeatFac(ParameterHeap::get()->matchWeightFeat)
  , matchVertFac(ParameterHeap::get()->matchWeightZ)
  , matchMaxDistance(ParameterHeap::get()->matchMaxFSpaceDist)
  , useFMatch(ParameterHeap::get()->regUseFeatMatch)
  , fMatchAvgDistThresh(ParameterHeap::get()->regFMatchAvgDistThresh)
  , icp4DPtCountThresh(ParameterHeap::get()->reg4DPtCountThresh)
  , icpAvgDistThresh(ParameterHeap::get()->regFMatchAvgDistThresh)
  , icpRegCoeff(ParameterHeap::get()->regIcpRegularizationCoeff)
  , icpNbCorresp(ParameterHeap::get()->regIcpNbCorresp)
  , icpErrOk(ParameterHeap::get()->regIcpMinErrOK)
  , icpMahalDistState(ParameterHeap::get()->regMaxRelStateDiff)
  , projDistTolerance(ParameterHeap::get()->projectionDistTolerance)
  , lidarMeasStdDevPerMeter(ParameterHeap::get()->lidarMeasStdDevPerMeter)
  , maxSampleDist(ParameterHeap::get()->framegenMaxDistance*0.8)
{
  const double minSamplingPtCount(ParameterHeap::get()->regMinSamplingCnt);
  const bool trackingEnabled(!ParameterHeap::get()->disableTracking);

  cout << endl << "creating registrator with ";
  if (useFMatch) cout << "activated"; else cout << "no";
  cout << " feature matching" << endl;

  // subsample each track, creating list "subsampledTracks"
  // this list is sorted after trackBufferIdx and pointCount
  SubTrackSPtrComp subtrackcompobj;
  int ti = 0;
  cout << "subsampling " << cf->longTermTracks.size() << " ltt-tracks..." << flush;
  BOOST_FOREACH(PointCloudTrack::SPtr track, cf->longTermTracks) {
    if (ti != 0) { // priority to non-world-tracks (otherwise world greedily masks all pixels out)
      cout << "[" << track->getUID() << "]" << flush;
      subsampledTracks.push_back(generate(track, cf->sttValidationSteps, ti, minSamplingPtCount));
    }
    ++ti;
  }
  subsampledTracks.sort(subtrackcompobj);
  PointCloudTrack::SPtr worldTrack = cf->getLttWorldTrack();
  if (worldTrack.get() != NULL) {
    cout << "[" << worldTrack->getUID() << "]" << flush;
    //subsampledTracks.push_back(generate(worldTrack, cf->sttValidationSteps, 0, minSamplingPtCount));
    cout << endl << endl << " ########### special world preparations!!!!! ##########" << endl;
    subsampledTracks.push_back(generateFromWorld(cf, worldTrack, cf->sttValidationSteps, 0, maxSampleDist));
  }
  // track in bi=cf->sttValidationSteps are contained in longTermTracks --> don't double-register
  for (int bi=cf->sttValidationSteps-1; bi>=0; --bi) { // newest tracks last -> projection image contains newest tracks (used for track linkage)
    cout << "subsampling " << cf->shortTermTracks[bi].tracks.size() << " stt-tracks..." << flush;
    ti = 0;
    SubTrackLst tmpSubtrackBuffer; // collect in set
    BOOST_FOREACH(PointCloudTrack::SPtr track, cf->shortTermTracks[bi].tracks) {
      if ((trackingEnabled) || (worldTrack.get() == NULL)) {
        tmpSubtrackBuffer.push_back(generate(track, bi, ti, minSamplingPtCount));
        ++ti;
      } else {
        track->update();
      }
    }
    tmpSubtrackBuffer.sort(subtrackcompobj);
    subsampledTracks.splice(subsampledTracks.end(), tmpSubtrackBuffer);
  }
  cout << "...done." << endl;

  // set iterator to first track
  currTrackRegistered = subsampledTracks.end();
  nextTrackToRegister = subsampledTracks.begin();
  lastTrackBufferIdx = INT_MAX;
  nbRegistered = 0;
  cf_->projSttIdx = UINT_MAX;
}

TrackRegistration::~TrackRegistration()
{
  flushTrackProjListAndExtend();
  if (fMatcher) delete fMatcher;
}

TrackRegistration::SubTrackSPtr TrackRegistration::generate(PointCloudTrack::SPtr track, int bi, int ti, double minSamplingPtCount)
{
  const double VIPWeight = ParameterHeap::get()->regIcpVIPWeight;
  unsigned int pcount = track->getPointCount();
  unsigned int subSampledPointCount = (pcount <= minSamplingPtCount) ? pcount : minSamplingPtCount + ceil(2*sqrt(pcount - minSamplingPtCount));
  subSampledPointCount = max(subSampledPointCount,1u); // minimum size 1
//    cout << "[" << pcount << "->" << subSampledPointCount << "]" << flush;

  SubTrackSPtr subtrack(new TrackSubsampled());
  subtrack->track = track;
  subtrack->trackBufferIdx = bi;
  subtrack->trackIdx = ti;
  subtrack->trackPtCnt = pcount;
  track->getRt2EgoCS(subtrack->RInit, subtrack->tInit);
  if (fabs(1.0-determinant(subtrack->RInit)) > 1e-8) {
    cout << "fixing R.." << flush;
    subtrack->RInit = orthonormalize(subtrack->RInit);
  }
  //subtrack->visibilityRatio = 0.0;
  subtrack->pointsTrackCS.resize(3,subSampledPointCount);
  subtrack->pointNormalsTrackCS.resize(3,subSampledPointCount);
  subtrack->initialWeights.resize(subSampledPointCount);
  subtrack->normalConfidence.resize(subSampledPointCount);
  subtrack->normalStdDevRAD.resize(subSampledPointCount);
  subtrack->pointCovarCS.reserve(subSampledPointCount);
  subtrack->normalCovarCS.reserve(subSampledPointCount);

  list<PointCloudTrack::Surface::SPtr> plist;
  if (subSampledPointCount == track->getPointCount()) {
    // don't subsample, take all points (avoids the same points to be chosen several times and others not)
    PointCloudTrack::PointIterator b,e;
    track->getPoints(b,e);
    while (b != e) {
      plist.push_back(b->s);
      ++b;
    }
  } else {
    // sample half of the points according to normal direction and half according to location (preserved details
    back_insert_iterator< list<PointCloudTrack::Surface::SPtr> > backInsIter(plist);
    unsigned int nbSampled = track->samplePointsUniformly(0.5*subSampledPointCount, backInsIter);
    nbSampled += track->sampleNormalsUniformly(subSampledPointCount - nbSampled, backInsIter);
    if (nbSampled != subSampledPointCount)
      cerr << endl << "WARNING: sampling (" << nbSampled << ") did not yield desired point count (" << subSampledPointCount << ")" << flush;
  }
  unsigned int indx1 = 0;
  BOOST_FOREACH(PointCloudTrack::Surface::SPtr s, plist) {
    DCMatrixCol(subtrack->pointsTrackCS, indx1) = s->position;
    DCMatrixCol(subtrack->pointNormalsTrackCS, indx1) = s->normal;
    subtrack->initialWeights(indx1) = s->isVIP ? VIPWeight : 1.0;
    subtrack->normalConfidence(indx1) = s->normalConfidence;
    subtrack->normalStdDevRAD(indx1) = s->normalStdDevRAD;
    subtrack->pointCovarCS.push_back(s->ptCovar);
    subtrack->normalCovarCS.push_back(s->nrCovar);
    ++indx1;
    if (indx1 == subSampledPointCount) break; // enough points sampled
  }
  //subtrack->visibilityRatio = (double)subtrack->visiblePts / (double)subSampledPointCount;
  return subtrack;
}

bool pixValid(int col, int row, FrameSPtr cf, double maxDist)
{
  LidarFrame::TrackIndex ti = cf->projLttTracks.get(col, row);
  assert(!std::isnan(cf->distance.get(col, row)) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->point3D.get(col, row)[0]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->point3D.get(col, row)[1]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->point3D.get(col, row)[2]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->normal3D.get(col, row)[0]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->normal3D.get(col, row)[1]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  assert(!std::isnan(cf->normal3D.get(col, row)[2]) && "TrackRegistration::generateFromWorld: unexpected nan in cf->distance");
  return    (cf->distance.get(col, row) < maxDist)
         && ((ti == UINT_MAX) || (ti == 0));
}

TrackRegistration::SubTrackSPtr TrackRegistration::generateFromWorld(FrameSPtr cf, PointCloudTrack::SPtr track, int bi, int ti, double maxDist)
{
  const unsigned int sampleCntMax = 4000;
  // determine valid pixels of frame
  int hSize, vSize;
  cf->point3D.getSize(hSize, vSize);
  unsigned int nbPix = hSize*vSize;
  unsigned int nbPix2 = nbPix/2;
  unsigned int nbValidPixUp = 0;
  unsigned int nbValidPixLo = 0;
  unsigned int pixIdx = 0;
  for (int col = 0; col < hSize; ++col) {
    for (int row = 0; row < vSize; ++row) {
      if (pixValid(col, row, cf, maxDist)) {
        if (pixIdx < nbPix2)
          ++nbValidPixUp;
        else
          ++nbValidPixLo;
      }
      ++pixIdx;
    }
  }
  unsigned int sampleCntUp = min(3*sampleCntMax/4, nbValidPixUp);
  unsigned int sampleCntLo = min(1*sampleCntMax/4, nbValidPixLo);
  unsigned int sampleCnt = sampleCntUp+sampleCntLo;
  double sampFacUp = (double) sampleCntUp / (double)nbValidPixUp;
  double sampFacLo = (double) sampleCntLo / (double)nbValidPixLo;
  double accSamplingUp = 1.0f;
  double accSamplingLo = 1.0f;
  // generate subtrack-object
  SubTrackSPtr subtrack(new TrackSubsampled());
  subtrack->track = track;
  subtrack->trackBufferIdx = bi;
  subtrack->trackIdx = ti;
  subtrack->trackPtCnt = track->getPointCount();
  track->getRt2TrkCS(subtrack->RInit, subtrack->tInit); // it's the INVERSE as compared to other sampled tracks!
  if (fabs(1.0-determinant(subtrack->RInit)) > 1e-8) {
    cout << "fixing R.." << flush;
    subtrack->RInit = orthonormalize(subtrack->RInit);
  }
  subtrack->pointsTrackCS.resize(3,sampleCnt);
  subtrack->pointNormalsTrackCS.resize(3,sampleCnt);
  subtrack->initialWeights.resize(sampleCnt);
  subtrack->normalConfidence.resize(sampleCnt);
  subtrack->normalStdDevRAD.resize(sampleCnt);
  subtrack->pointCovarCS.reserve(sampleCnt);
  subtrack->normalCovarCS.reserve(sampleCnt);

  // determine unwarping parameters and subsample:
//  unsigned int hbSize = hSize / nbBlocks2Unwarp; // might be rounded off
//  double fDec = (double)hbSize / (double)hSize; // = 1/realNbBlocks
//  DVector poseDiff = HTM::HTM_2_YawPitchRollXYZ(HTM::invert_HTM(frame->moveHTM));
//  DMatrix Rci; DVector tci;
//  if (unwarp)
//    cout << "unwarping by " << poseDiff << "..." << flush;
//  subPtsHTM = scanInitTransform; // points should be specified relative to end-position
//  double fac = 1.0;
  unsigned int indx1 = 0;
  for (int col = 0; col < hSize; ++col) {
//    if ((unwarp) && (col % hbSize == 0)) {
//      HTM::YawPitchRollXYZ_2_Rt(fac*poseDiff, Rci, tci);
//      fac = max(0.0,fac-fDec); // decrease influence of last registration while approaching end position
//    }
    double *accSampling = &accSamplingUp;
    double *sampFac = &sampFacUp;
    for (int row = 0; row < vSize; ++row) {
      if (indx1 == sampleCnt) break; // enough points sampled
      if (row == 32) {
        accSampling = &accSamplingLo;
        sampFac = &sampFacLo;
      }
      if (pixValid(col, row, cf, maxDist)) {
        *accSampling += *sampFac; // only use valid pixels
        if (*accSampling >= 1.0f) {
          *accSampling -= 1.0f;
//          if (unwarp) { // unwarp using lastRegCorrection
//            pts.push_back(ublas::prod(Rci,frame->point3D.get(col, row))+tci);
//            nrs.push_back(ublas::prod(Rci,frame->normal3D.get(col, row)));
//          } else {
            DCMatrixCol(subtrack->pointsTrackCS, indx1) = cf->point3D.get(col, row);
            DCMatrixCol(subtrack->pointNormalsTrackCS, indx1) = cf->normal3D.get(col, row);
//          }
          subtrack->initialWeights(indx1) = 1.0; // no VIPs
          subtrack->normalConfidence(indx1) = cf->normalConfidence.get(col, row);
          subtrack->normalStdDevRAD(indx1) = cf->normalStdDevRAD.get(col, row);
          subtrack->pointCovarCS.push_back(cf->pointVariance3D.get(col, row));
          subtrack->normalCovarCS.push_back(cf->normalVariance3D.get(col, row));
          ++indx1;
        }
      }
    }
  }

  return subtrack;
}


ICP::EnergyFunction::Surface TrackRegistration::convert(const PointCloudTrack::Surface &s)
{
  ICP::EnergyFunction::Surface ret;
  ret.p = s.position;
  ret.n = s.normal;
  ret.nConfidence = s.normalConfidence;
  ret.nStdDevRAD = s.normalStdDevRAD;
  ret.pCovar3D = s.ptCovar;
  ret.nCovar3D = s.nrCovar;
  return ret;
}

ICP::EnergyFunction::Surface TrackRegistration::convert(SearchResultsBase &sr)
{
  ICP::EnergyFunction::Surface ret;
  ret.p = sr.getCorrespPoint();
  ret.n = sr.getCorrespNormal();
  ret.nConfidence = sr.getNormConf();
  ret.nStdDevRAD = sr.getNormStdDevRAD();
  ret.pCovar3D = sr.getPCovar();
  ret.nCovar3D = sr.getNCovar();
  return ret;
}

std::pair<double,double> TrackRegistration::calculateError(PointCloudTrack::SPtr track, const mdefs::DVector &state)
{
  DMatrix R; DVector t;
  PointCloudTrack::state2Rt_T2E(state, R, t);
  PointCloudTrack::PointIterator pcurr, pend;
  track->getPoints(pcurr, pend);
  ICP::EnergyFunction *en1 = new ICP::Point2PlaneEnergy();
  ICP::EnergyFunction *en2 = new ICP::Point2PointEnergy();
//  ICP::EnergyFunction *en3 = new ICP::NormalConfLinCombEnergy(*en1, *en2);
  double err1 = 0.0;
  double err2 = 0.0;
  double cumw = 0.0;
  for (; pcurr != pend; ++pcurr) { // for each point in track
    try {
      ICP::EnergyFunction::Surface s = convert(*(pcurr->s));
      s.p = prod(R,s.p) + t; // transform to Ego-CS
      s.n = prod(R,s.n); // transform to Ego-CS
      // s.pCovar3D = RT * covar * R
      // s.nCovar3D = R * covar * RT
      SearchResultsBase sr = pnSearch.findClosestNeighbor(s.p, s.n, s.nConfidence, 1);
      ICP::EnergyFunction::Surface sNN = convert(sr);
      err1 += en1->calculate(s, sNN) * sr.getWeight();
      err2 += en2->calculate(s, sNN) * sr.getWeight();
      cumw += sr.getWeight();
    } catch (std::exception &e) { // do nothing, just ignore this
    }
  }
  if (cumw < 0.001)
    return make_pair(DBL_MAX, DBL_MAX);
  else
    return make_pair(err1/cumw, err2/cumw);
}

// given a center in the state-CS, this returns the translated state with the CS-origin at the center
DVector centerState(DVector state6D, const DVector &center)
{
  DMatrix Rot; DVector trans;
  HTM::YawPitchRollXYZ_2_Rt(state6D(2), state6D(1), state6D(0), state6D(3), state6D(4), state6D(5), Rot, trans);
  trans += ublas::prod(Rot, center);
  //HTM::Rt_2_YawPitchRollXYZ(Rot, trans, state6D(2), state6D(1), state6D(0), state6D(3), state6D(4), state6D(5));
  //since only translation changes, apply this manually (=faster)
  state6D(3) = trans(0);
  state6D(4) = trans(1);
  state6D(5) = trans(2);
  return state6D;
}

bool TrackRegistration::nextTrack()
{
  if (nextTrackToRegister == subsampledTracks.end()) return false;

  SubTrackSPtr subtrack = *nextTrackToRegister;
  PointCloudTrack::SPtr track = subtrack->track;
  double tNormConf = track->getAvgNormalConfidence();
  ParameterHeap* params = ParameterHeap::get();

  // subsampledTracks contains in order: [ltt-moving-tracks, world-track], [old-stts], ..., [newest-stts]
  // lastTrackBufferIdx is initially INT_MAX, the following if should hence execute when the first track is registered
  if (lastTrackBufferIdx != subtrack->trackBufferIdx) {
    cout << endl << "------------- clearing projection buffer (trackBufferIdx changes from " << lastTrackBufferIdx << " to " << subtrack->trackBufferIdx << ")" << flush;
    flushTrackProjListAndExtend(); // project unsuccessful-ICP-registration tracks, will fill cf->extProjTracks
    if (lastTrackBufferIdx == (int)cf->sttValidationSteps) { // switching from long-term-tracks to old-ssts
      cf->projLttTracks.copyShallowFrom(cf->extProjTracks); // store projection buffer
    }
    if (lastTrackBufferIdx == 1) { // switching from [1] to [0]
      cf->projLCTracks.copyShallowFrom(cf->extProjTracks); // store projection buffer
    }
    // flushTrackProjListAndExtend() is also called just before segmentation is invoked
    // hence, when segmentation is invoked, we have:
    // - cf->projLttTracks  containing track-indices to cf->longTermTracks
    // - cf->projLCTracks   containing track-indices to cf->shortTermTracks[1]
    // - cf->extProjTracks  containing track-indices to cf->shortTermTracks[0]
    lastTrackBufferIdx = subtrack->trackBufferIdx;
    cf->projSttIdx = subtrack->trackBufferIdx;
    cf->regAvailPixels.fill(true); // each point in next frame can only be used once for correspondence
    cf->projTracks.fill(UINT_MAX);
    cf->projTrackDist.fill(DBL_MAX);
    cf->projTrackCnt.fill(0);
    cf->matchedColNextFrame.fill(-1);
    cf->matchedRowNextFrame.fill(-1);
  }


  cout << endl << "          track[" << track->getUID() << "](" << subtrack->trackPtCnt << "->" << subtrack->pointsTrackCS.size2() << ")" << flush;
  /////////////////////////////////////////////////////////////////////
  ///////// Step 1: Setup ICP                        //////////////////
  /////////////////////////////////////////////////////////////////////
  unsigned int minIter(4); // changing this value might require changing buffer access below
  unsigned int maxIter(10);
  double nConfThresh(0.7);

  unsigned int segSubSampledPointCount = subtrack->pointsTrackCS.size2();
  DCMatrixColConstIterator pB(subtrack->pointsTrackCS, 0);
  DCMatrixColConstIterator pE(subtrack->pointsTrackCS, segSubSampledPointCount);
  DCMatrixColConstIterator nB(subtrack->pointNormalsTrackCS, 0);
  DCMatrixColConstIterator nE(subtrack->pointNormalsTrackCS, segSubSampledPointCount);
  BOOST_AUTO(wB, subtrack->initialWeights.begin());
  BOOST_AUTO(wE, subtrack->initialWeights.end());
  BOOST_AUTO(nCB, subtrack->normalConfidence.begin());
  BOOST_AUTO(nCE, subtrack->normalConfidence.end());
  BOOST_AUTO(nSB, subtrack->normalStdDevRAD.begin());
  BOOST_AUTO(nSE, subtrack->normalStdDevRAD.end());
  BOOST_AUTO(pcvB, subtrack->pointCovarCS.begin());
  BOOST_AUTO(pcvE, subtrack->pointCovarCS.end());
  BOOST_AUTO(ncvB, subtrack->normalCovarCS.begin());
  BOOST_AUTO(ncvE, subtrack->normalCovarCS.end());
  ICP::NeighborSearch<SearchResultsBase>::ncpSearch ncpFunc;

  //// Besl and Chen don't have good covariance-estimators integrated!
  ////  ICP::ICPBesl<DCMatrixColConstIterator> icp(cpwFunc, pB, pE, nB, nE);
  ////  ICP::ICPChen<DCMatrixColConstIterator> icp(cpnFunc, pB, pE, nB, nE);
  ////  ICP::ICPChen<DCMatrixColConstIterator,ProjectedNeighborSearch::SearchResults> icp(pB, pE, nB, nE, ncpFunc, icpNbCorresp, icpRegCoeff);
  typedef ICP::ICPLinearized< std::vector<mdefs::DMatrix>::iterator, DCMatrixColConstIterator, DVector::const_iterator, SearchResultsBase> LinICPT;

  TrackNeighborSearch *wtNNSearch = NULL;
  ICP::EnergyFunction *enLin = NULL;
  LinICPT::EstimationMode esMode = LinICPT::Mode3R3T;
  bool specialWorldTreatment = ((track->getUID() == 1ull) && (subtrack->trackBufferIdx == (int)cf->sttValidationSteps));
  if (specialWorldTreatment) { // world track
    cout << endl << endl << " ########### special world treatment!!!!! ##########" << endl;
//    if (! ParameterHeap::get()->disableLocalization)
    wtNNSearch = new TrackNeighborSearch(track);
    ncpFunc = boost::bind(&TrackNeighborSearch::findNeighbors, wtNNSearch, _1, _2, _3, _4);
    enLin = &enPtPl;
  } else {
    //  ICP::cpSearch  corFunc = boost::bind(&ProjectedNeighborSearch::findClosestProjectedNeighbor, &projNeighbSearch, _1, _2, _3);
    //  ICP::cpwSearch corWeightFunc = boost::bind(&ICPBesl::threshWeightCorresp, _1, _2, _3, _4, corFunc, 5.0);
    //  ICP::cpnwSearch cpnwFunc = boost::bind(&ProjectedNeighborSearch::findClosestProjectedNeighbor, &pnSearch, _1, _2, _3, _4, _5, _6);
    //  ICP::cpwSearch cpwFunc = boost::bind(&ICP::cpnwSearch2cpwSearch, _1, _2, _3, _4, _5, cpnwFunc);
    //  ICP::cpnSearch cpnFunc = boost::bind(&ICP::cpnwSearch2cpnSearch, _1, _2, _3, _4, _5, cpnwFunc);
    ncpFunc = boost::bind(&ProjectedNeighborSearch::findClosestNeighbor, &pnSearch, _1, _2, _3, _4);
    if (tNormConf > nConfThresh) {
      cout << " nConf:" << tNormConf << " --> using PtPl" << flush;
      enLin = &enPtPl;
      pnSearch.useNormalWeighting(params->regIcpUseNormalWeight);
    } else {
      cout << " nConf:" << tNormConf << " --> using PtPt" << flush;
      enLin = &enPtPt;
      pnSearch.useNormalWeighting(false);
      minIter *= 2;
      maxIter *= 2;
    }
    esMode = (subtrack->trackPtCnt < icp4DPtCountThresh) ? LinICPT::Mode1R3T : LinICPT::Mode3R3T;
  }

  LinICPT  icp(pB, pE, nB, nE, wB, wE, nCB, nCE, nSB, nSE, pcvB, pcvE, ncvB, ncvE, ncpFunc, *enLin, icpNbCorresp, icpRegCoeff, esMode);
//  const unsigned int nbIterSwitchEnergy(1);
//  BOOST_AUTO(sw, boost::bind(&ProjectedNeighborSearch::useNormalSearch, &pnSearch, _1));
//  ICP::ICPSwitch< std::vector<mdefs::DMatrix>::iterator, DCMatrixColConstIterator, DVector::const_iterator, ProjectedNeighborSearch::SearchResults>
//    icp(pB, pE, nB, nE, nCB, nCE, nSB, nSE, pcvB, pcvE, ncvB, ncvE, ncpFunc, *en1, *en2, nbIterSwitchEnergy, sw, true, icpNbCorresp, icpRegCoeff);
  icp.reset(subtrack->RInit, subtrack->tInit);
  track->histRegICPError.clear();
  subtrack->icpCorrespondenceBuffer.clear();

  /////////////////////////////////////////////////////////////////////
  ///////// Step 2: Feature Matching if prediction failed /////////////
  /////////////////////////////////////////////////////////////////////
  if (track->getAge() == 0) {
    cout << " A0 " << flush;
    pnSearch.useHorizSearch(true); // implicitly sets useNormalSearch(false)
  } else {
    pnSearch.useHorizSearch(false); // implicitly sets useNormalSearch(false)
    pnSearch.useNormalSearch(tNormConf > nConfThresh);
  }
  pnSearch.correspondenceBuffer.clear(); // DEBUG
  icp.establishCorrespondences();
  subtrack->icpCorrespondenceBuffer.push_back(std::list<mdefs::DVector>()); // add empty list
  subtrack->icpCorrespondenceBuffer.back().swap(pnSearch.correspondenceBuffer); // move all contents of correspSearch into this empty list
  track->histRegICPError.push_back(icp.getLastAvgWeightDist());
  track->histRegUsedFeatureMatching = false;
  if ((useFMatch) && (icp.getLastAvgWeightDist() > fMatchAvgDistThresh)) {
    throw range_error("TrackRegistration::nextTrack: Feature Matching was requested but is deactivated in code");
  }
  /////////////////////////////////////////////////////////////////////
  ///////// Step 3: Refine estimation of R+t by running ICP ///////////
  /////////////////////////////////////////////////////////////////////
  const double maxRelAvgDistChange(0.01);
  //const double minExpErr(0.0005); // expected error (always accept if result-dist smaller)
  const double maxRelErrorInc(1.2); // allow maximum 20% error increase (only if dist > minExpErr)
  const unsigned int buffSize(3);
  double avgDistChange[buffSize] = {1+maxRelAvgDistChange,1+maxRelAvgDistChange,1+maxRelAvgDistChange};
  unsigned int avgDistIdx = 0;
//  double firstAvgDist = icp.getLastAvgDist(); // this error is independent of the energy used
  double firstAvgDist = icp.getLastAvgWeightDist(); // this error is independent of the energy used
//  double firstAvgDist = icp.getLastAvgWeightError(); // as the error is minimized, this should be used as stop criterion
  double oldAvgDist = firstAvgDist;
  double currAvgDist = 0.0;
  bool icpSucceeded = false;
  cout << " fd=" << firstAvgDist << " " << flush;
  for (unsigned int i=1; i<maxIter; ++i) {
    icpSucceeded = icp.minimizeError(); // correspondences were already established
    icpSucceeded = icp.establishCorrespondences() && icpSucceeded;
    if (icpSucceeded)
      cout << "." << flush;
    else
      cout << "x" << flush;
//    currAvgDist = icp.getLastAvgDist();
    currAvgDist = icp.getLastAvgWeightDist();
//    currAvgDist = icp.getLastAvgWeightError();
    if (i == 1) { firstAvgDist = currAvgDist; cout << " fd=" << firstAvgDist << " " << flush; }
    subtrack->icpCorrespondenceBuffer.push_back(std::list<mdefs::DVector>()); // add empty list
    subtrack->icpCorrespondenceBuffer.back().swap(pnSearch.correspondenceBuffer); // move all contents of correspSearch into this empty list
    if (std::isnan(currAvgDist)) cout << " AD=NAN" << flush;
    //cout << icp.t << " " << icp.R << endl;
//    cout << ":d=" << newAvgDst << flush;
    track->histRegICPError.push_back(icp.getLastAvgWeightDist());
    avgDistChange[avgDistIdx] = (currAvgDist == 0.0) ? 0.0 : abs(oldAvgDist-currAvgDist)/(oldAvgDist+currAvgDist);
    cout << " d=" << currAvgDist << "(delta:" << avgDistChange[avgDistIdx] << ") " << flush;
    if (i+1 >= minIter) {
      if (   (avgDistChange[avgDistIdx] < maxRelAvgDistChange) // if avgDist doesn't change much
          && (avgDistChange[(avgDistIdx-1+buffSize)%buffSize] < maxRelAvgDistChange)
          && (avgDistChange[(avgDistIdx-2+buffSize)%buffSize] < maxRelAvgDistChange)
          && (icpSucceeded)) {
        cout << " (i=" << i << ")" << flush;
        break; // stop iterating, exit
      }
      if ((currAvgDist > icpErrOk) && (currAvgDist > maxRelErrorInc*firstAvgDist)) { // at most 20% error increase
        cout << " (i=" << i << ")" << flush;
        icpSucceeded = false;
        break;
      }
    }
    oldAvgDist = currAvgDist;
    avgDistIdx = (avgDistIdx+1)%buffSize;
  }

  DVector z(6);
  icp.getEstimate(z);

  delete wtNNSearch;
  if (specialWorldTreatment) { // world track
    // invert ICP estimate
    DMatrix Rot; DVector trans;
    icp.getEstimate(Rot, trans);
    HTM::invert_Rt(Rot, trans);
    HTM::Rt_2_YawPitchRollXYZ(Rot, trans, z(2), z(1), z(0), z(3), z(4), z(5));
  }

  DVector trackState = track->getState();
  DMatrix trackCovar = track->getCovar();
  if (icpSucceeded) {
//    if ((specialWorldTreatment)
//        && (   (boost::filesystem::path(cf->sourcefile).filename() == "scan001470.png")
//            || (boost::filesystem::path(cf->sourcefile).filename() == "scan001360.png"))) {
//      cout << "Model center" << icp.getModelCenter() << flush;
//    }
    // if the measurement is a very unprobable state, mark it as invalid
    DVector modelCenter = specialWorldTreatment ? icp.getTransModelCenter() : icp.getModelCenter();
    DVector trackZState = DVectorRange(trackState, ublas::range(0,6));
    DMatrix trackZCovar = DMatrixRange(trackCovar, ublas::range(0,6), ublas::range(0,6));
    DSMatrix trackZCovarI = invSym(trackZCovar, 1e-12);
    //DVector stateDiff = z - trackZState; // when center is far away, small rotations already lead to big differences in translation
    DVector stateDiff = centerState(z, modelCenter) - centerState(trackZState, modelCenter); // center-based calculation does not have this problem
    // normal distribution: exp(-0.5 * (z-mu)^T * Covar^-1 * (z-mu))/(pow(2*M_PI,ndim/2)*sqrt(det(Covar)));
    // 68% are within 1sigma, 90% are within 1.64sigma, 95% within 2sigma, 99.7% within 3sigma
    // calculate probability with mahalanobis distance: sqrt((z-mu)^T * Covar^-1 * (z-mu))
    // MD is eqivalent to (z-mu)/stdDev and hence 1 at stdDev
    DVector tmpVec = ublas::prod(stateDiff, trackZCovarI);
    double distMahaSq = ublas::inner_prod(tmpVec, stateDiff);
    if (sqrt(distMahaSq) > icpMahalDistState)
      icpSucceeded = false;
//    if (specialWorldTreatment) // TODO: remove HACK!
//      icpSucceeded = true;
    //DMatrix tmpMat = ublas::prod(trackZCovar, trackZCovarI);
//    double detCovar = determinant(trackZCovar);
//    double probab = exp(-0.5*distMahaSq)/(pow(2*M_PI,6/2)*sqrt(detCovar));
    // probab should be in [0,1] but due to numerical instabilities it might be > 1
//    if (probab < 0.3)
//      icpSucceeded = false;
    cout << endl << "z       " << z;
    cout << endl << "x       " << trackZState;
    cout << endl << "covar   " << trackZCovar;
    cout << endl << "covarI  " << trackZCovarI;
    //cout << endl << "I6      " << tmpMat;
    cout << endl << "x-z     " << stateDiff;
    cout << endl << "[[[[ dist^2  " << distMahaSq << " ]]]] " << endl;
//    cout << endl << "dist^2  " << distMahaSq;
//    cout << endl << "det     " << detCovar;
//    cout << endl << "pow     " << pow(2*M_PI,6/2);
//    cout << endl << "[[[[ p  " << probab << " ]]]] " << endl;
  }

  if (!icpSucceeded) {
    cout << " failed" << flush;
    track->update();
  } else {
    cout << " success" << flush;
    // set variance of measurement = var_registration + var_association
    // var_registration is obtained directly from ICP plus the measurement variance of the Lidar scanner (at 10m)
    // var_association is calculated from mahalanobis-dist(=std.dev.) from resulting translation vector and track-variance
    //                 and the number of points contained in the track to weight this association-variance
    DMatrix icpCovar; // delta[rx-roll,ry-pitch,rz-yaw,tx,ty,tz]
//    icpCovar = icp.getAnalyticalVariance();
//    icpCovar = icp.getEstimatedVariance();
    icpCovar = /*0.5**/(icp.getEstimatedVariance() + icp.getAnalyticalVariance());
//    icpCovar = icp.getSampledVariance(true); // correspondences re-searched
//    icpCovar /= 4*track->pointNormalDistributionRatio();
//    double addCovar = pow(norm_2(icpTrans),2)*exp(-(double)(segSubSampledPointCount-5)/50.0);
    if ((ParameterHeap::get()->disableLocalization) && (track->getUID() == 1)) { //TODO(9): increase efficiency by completely skipping ICP and kdTree-construction
      z = DVectorRange(trackState, ublas::range(0,6));
      icpCovar = DIdMatrix(6) * ParameterHeap::get()->lidarDistStdDev;
    }
    track->update(z,icpCovar);
    // Debug:
//    DVector finalState;
//    track->getState(finalState);
//    DMatrix finalCovar;
//    track->getCovar(finalCovar);
//    cout << endl <<"...update:"
//             << " z=" << z
//             << ", covarDiag=" << matrix_vector_range<DMatrix>(icpCovar, range (0, 6), range (0, 6))
//             << ", covar=" << icpCovar
//             << ", wAvgD=" << icp.getLastAvgWeightDist()
//             << ", avgW=" << icp.getAvgWeight()
//             << ", validW=" << icp.getSignificantWeightRatio()*100 << "%"
//             << ", outliers=" << icp.getOutlierRatio()*100 << "%";
//    cout << endl << "--->";
//    cout << " state=" << finalState;
//    cout << " covar=" << finalCovar << flush;
//    cout << "...done" << flush;
  }

  /////////////////////////////////////////////////////////////////////
  ///////// Step 4 (Postprocessing): Unmask correspondences ///////////
  /////////////////////////////////////////////////////////////////////
  if (icpSucceeded) {
    projectToImage(subtrack, true); //removeInvalidPoints
  } else {
    trackProjList.push_back(subtrack);
  }
  
  currTrackRegistered = nextTrackToRegister;
  nextTrackToRegister++;
  nbRegistered++;
  return (nextTrackToRegister != subsampledTracks.end());
  //return ((nextTrackToRegister != subsampledTracks.end()) && (icpSucceeded));
}

void TrackRegistration::projectToImage(SubTrackSPtr subtrack, bool icpSuccessful)
{
  // Mark all corresponding points in new frame, so they are not used for further tracks
  // In the same step, project TrackIDs
  PointCloudTrack::SPtr track = subtrack->track;
  PointCloudTrack::PointIterator pcurr, pend;
  track->getPoints(pcurr, pend);
  DMatrix R; DVector t;
  track->getRt2EgoCS(R, t); // points must be transformed in Ego-CS
  //unsigned int nbErased = 0;
  while (pcurr != pend) { // for each point in track
    DVector pt = pcurr->s->position;
    DVector nr = pcurr->s->normal;
    double& nc = pcurr->s->normalConfidence;
    assert(!std::isnan(pt(0)) && "pt(0) is NAN");
    assert(!std::isnan(pt(1)) && "pt(1) is NAN");
    assert(!std::isnan(pt(2)) && "pt(2) is NAN");
    assert(!std::isnan(nr(0)) && "nr(0) is NAN");
    assert(!std::isnan(nr(1)) && "nr(1) is NAN");
    assert(!std::isnan(nr(2)) && "nr(2) is NAN");
    pt = prod(R,pt) + t; // transform to Ego-CS
    nr = prod(R,nr); // transform to Ego-CS
    assert(!std::isnan(pt(0)) && "pt'(0) is NAN");
    assert(!std::isnan(pt(1)) && "pt'(1) is NAN");
    assert(!std::isnan(pt(2)) && "pt'(2) is NAN");

    double ptDist; int ptCol, ptRow;
    bool incrementedPCurr = false;
    if (projector.getImageIndexRel(pt(0), pt(1), pt(2), ptCol, ptRow, ptDist)) { // projection yields valid pixel
      double& nConf = nc;
//      double& nConf = cf->normalConfidence.get(ptCol,ptRow);
//      double  nConf = min(nc, cf->normalConfidence.get(ptCol,ptRow));
      double& pixDst = cf->distance.get(ptCol,ptRow);
      DVector& pixPt = cf->point3D.get(ptCol,ptRow); // only valid if dist!=DBL_MAX
//      DVector& pixNr = cf->normal3D.get(ptCol,ptRow); // only valid if dist!=DBL_MAX
      double& pixPDist = cf->projTrackDist.get(ptCol,ptRow);
      unsigned int& pixPCnt = cf->projTrackCnt(ptCol,ptRow);
      bool pixValid = (pixDst != DBL_MAX);
      const double nConfThresh = 0.3; // used twice
      double pixPlaneDist = (nConf < nConfThresh) ? // positive if in front of plane
                            pixDst - ptDist :
                            inner_prod(pt-pixPt, nr); // assume that normal vector of map is relevant
//                            inner_prod(pt-pixPt, pixNr); // assume that normal vector at pixel is relevant
      // if projected point so far closest point on that pixel and lies not far behind current scene point -> overwrite projection values
      if (   (pixValid)
          && (icpSuccessful || (pixPCnt == 0))
          && (ptDist < pixPDist) // point is in front of previously projected points
          && (pixPlaneDist <= projDistTolerance*ptDist/10) // 0.2m tolerance @ 10m
          && (pixPlaneDist > -projDistTolerance*ptDist/10)
          ) {
        cf->projTracks(ptCol,ptRow) = subtrack->trackIdx;
        pixPDist = ptDist; // overwrite with smaller value
        if (icpSuccessful) {
          cf->regAvailPixels(ptCol, ptRow) = false; // mark that matched point was already used
          pixPCnt += 1;
        }
      }
      // if point is way in front of plane, remove it (point was not but should have been seen from scanner)
      double remThreshAngle = (track->getUID() == 1) ? 30. : 10.;
      double remThresh = cos((90.-remThreshAngle)/180.*M_PI);
      double cosPtPlAngle = fabs(inner_prod(-pt/norm_2(pt), nr));
//      double cosPtPlAngle = fabs(inner_prod(-pixPt/norm_2(pixPt), pixNr));
      if (icpSuccessful
          && (pixValid)
          && (pixPlaneDist > projDistTolerance*ptDist/10)
          && ((nConf < nConfThresh) || (cosPtPlAngle > remThresh)) // skip if flat area with small angle towards scanner
          && (subtrack->trackBufferIdx == (int)cf->sttValidationSteps)) { // only for long-term-tracks
        //if ((track->getUID() == 1) && (track->getAge() > debugBreak) && (nConf >= nConfThresh))
        //  nbErased++;
        pcurr = track->erasePoint(pcurr);
        incrementedPCurr = true;
      }
    }  // end: if projection lies within image
    if (!incrementedPCurr)
      ++pcurr;
  } // end: for all points in segment
  //if (nbErased > 0)
  //  cout << "...erased " << nbErased;
  cout << "...maskOK" << flush;
}

void TrackRegistration::flushTrackProjListAndExtend()
{
  cout << endl << "projecting remaining " << trackProjList.size() << " tracks" << flush;
  BOOST_FOREACH(SubTrackSPtr subtrack, trackProjList)
    projectToImage(subtrack, false);
  trackProjList.clear();
  extendProjTrackIDs();
}

// interpolate missing pixels and increase borders
void TrackRegistration::extendProjTrackIDs()
{
  int hsize, vsize;
  cf->projTracks.getSize(hsize,vsize);
  LidarFrame::TrackIndex trk[4];

  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      LidarFrame::TrackIndex toAssign = cf->projTracks(col,row);
      if ((toAssign == UINT_MAX) && (cf->distance(col,row) != DBL_MAX)) {
        // no track assigned but pixel valid --> search 4-neighborhood
        trk[0] = cf->projTracks(col-1,row);
        trk[1] = cf->projTracks(col+1,row);
        trk[2] = cf->projTracks(col,row-1);
        trk[3] = cf->projTracks(col,row+1);
        sort(trk, trk+4);
        if (trk[1] == trk[2]) {
          toAssign = trk[1]; // due to sort 0 and 3 can only be same if equal to 1 and 2
        } else {
          if ((trk[0] == trk[1]) && (trk[2] != trk[3])) {
            toAssign = trk[0];
          } else {
            if ((trk[2] == trk[3]) && (trk[0] != trk[1]))
              toAssign = trk[2];
          }
        }
      }
      cf->extProjTracks.set(col,row, toAssign);
    } // loop over cols
  } // loop over rows
}

void calculateMatches(FrameSPtr lf, FrameSPtr cf)
{
  FeatureMatching fm(cf, lf);

  int hsize, vsize;
  lf->distance.getSize(hsize,vsize);
  int colN, rowN; double dist;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      try {
        fm.findClosestNeighbor(col, row, colN, rowN, dist); // returns col,row,dist of found neighbor. throws on error
        if ((colN >= 0) && ((lf->visibleInNextFrame.get(col,row)) || (cf->visibleInPrevFrame.get(colN,rowN)))) {
          lf->matchedColNextFrame.set(col,row,colN);
          lf->matchedRowNextFrame.set(col,row,rowN);
        } else {
          lf->matchedColNextFrame.set(col,row,-1);
          lf->matchedRowNextFrame.set(col,row,-1);
        }
      } catch (exception &e) {
      }
    }
  }
}


