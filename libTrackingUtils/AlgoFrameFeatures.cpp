#include "AlgoFrameFeatures.hpp"

#include <stdexcept>
#include <boost/foreach.hpp>

#include "ParameterHeap.hpp"
#include "LidarImageFeatures.hpp"

using namespace std;
using namespace matrixTools;


void postprocessDistanceData(FrameSPtr frame, const LidarImageProjector &projector)
{
  if (!frame->valid)
    throw logic_error("postprocessDistanceData: frame data is not valid");
  ParameterHeap* params = ParameterHeap::get();

  LidarImageFeatures::postprocessDistanceImage(
      frame->distance,
      params->framegenMaxDistance,
      params->framegenInterpolHPixThresh,
      params->framegenInterpolVPixThresh,
      params->framegenInterpolMetPixThresh,
      params->framegenSmoothDstImg,
      params->framegenSmoothDstThresh,
      &frame->tmpFloat1);
  LidarImageFeatures::transformTo3D(
      frame->distance,
      frame->point3D,
      projector,
      &frame->pointVariance3D,
      params->lidarDistStdDev,
      params->lidarHAngStdDevRAD,
      params->lidarVAngStdDevRAD);

  ////////////////////////////////////////////////////////////////
  //////////////////   filter below-ground  //////////////////////
  ////////////////////////////////////////////////////////////////
//  for (int row = 0; row < vsize; ++row) {
//    for (int col = 0; col < hsize; ++col) {
//      if (frame->point3D.get(col,row)(2) < -0.5) {
//      }
//    }
//  }
}

void calculateDerivative(FrameSPtr frame)
{
  LidarImageFeatures::derivative(frame->distance, frame->distDerivativeH, frame->distDerivativeV);
}

void calculateDistDiff(FrameSPtr lf, FrameSPtr cf)
{
  LidarImageFeatures::difference(cf->distance, lf->distance, cf->distanceDiffLastFrame);
}

void calculateConnections(FrameSPtr frame)
{
  ParameterHeap* params = ParameterHeap::get();

  LidarImageFeatures::connectionWeigths(
      frame->distance,
      frame->distDerivativeH,
      frame->distDerivativeV,
      frame->connectivityH,
      frame->connectivityV,
      params->connMaxDst,
      params->connMaxDstRatio,
      params->connFacA,
      params->connFacB,
      params->connFacC,
      params->connFacD);
  if (params->segVariantTilt == 409) { // IV09 method --> binarize connections
    int hsize, vsize;
    // horizontally
    frame->connectivityH.getSize(hsize, vsize);
    for (int row = 0; row < vsize; ++row) {
      for (int col = 0; col < hsize; ++col) {
        frame->connectivityH(col, row) = (frame->connectivityH(col, row) >= 0.5) ? 1 : 0;
      }
    }
    // vertically
    frame->connectivityV.getSize(hsize, vsize);
    for (int row = 0; row < vsize; ++row) {
      for (int col = 0; col < hsize; ++col) {
        frame->connectivityV(col, row) = (frame->connectivityV(col, row) >= 0.5) ? 1 : 0;
      }
    }
  }
}

void calculateNormals(FrameSPtr frame, const LidarImageProjector &projector, bool full, bool verbose)
{
  ParameterHeap* params = ParameterHeap::get();


  LidarImage<double> *sdi = full ? &frame->normalStdDevRAD : NULL;
  LidarImage<matrixTools::DMatrix> *nvi = full ? &frame->normalVariance3D : NULL;
  LidarImageFeatures::normals(
      frame->distance,
      frame->point3D,
      frame->normal3D,
      &frame->connectivityH,
      &frame->connectivityV,
      sdi,
      nvi,
      &frame->normalConfidence,
      &frame->tmpVec3D,
      &frame->tmpFloat1,
      &projector,
      2*params->lidarMeasStdDevConst,
      params->featureConnectWeightCreate,
//      params->featureConnectWeightSmooth,
      params->lidarHAngResolRAD,
      params->lidarVAngResolRAD,
      params->lidarDistStdDev,
      params->lidarHAngStdDevRAD,
      params->lidarVAngStdDevRAD,
      params->featurenNbConfMedianPasses,
      params->featureConfNeighb,
      params->featureConfMin,
      verbose);
}

void calculateMatchingFeatures(FrameSPtr frame)
{
  int hsize, vsize;
  frame->distance.getSize(hsize,vsize);
  const double limit = 5.0; // some constant for non-valid distances

  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      DVector fv(4);
      double d = frame->distance.get(col,row);
      if (d != DBL_MAX) {
        double dl = frame->distance.get(col-1,row);
        double dr = frame->distance.get(col+2,row);
        double du = frame->distance.get(col,row-1);
        double dd = frame->distance.get(col,row+2);
        fv(0) = max(min(dl-d, limit), -limit);
        fv(1) = max(min(dr-d, limit), -limit);
        fv(2) = max(min(du-d, limit), -limit);
        fv(3) = max(min(dd-d, limit), -limit);
//        fv.addElement( dl!=DBL_MAX ? dl - d : limit); 
//        fv.addElement( dr!=DBL_MAX ? d - dr : limit);
//        fv.addElement( du!=DBL_MAX ? du - d : limit);
//        fv.addElement( dd!=DBL_MAX ? d - dd : limit);
      }
      frame->matchingFeatures.set(col,row,fv);
    }
  }
}

void visibilityNextFrame(FrameSPtr cf, FrameSPtr nf, LidarImageProjector &projector)
{
  LidarImage<bool> &visibilityResult = cf->visibleInNextFrame;
  visibilityOtherFrame(cf, nf, projector, visibilityResult);
//  BOOST_FOREACH(LidarSegment::SPtr currSeg, cf->segmentList) {
//    unsigned int nbVisible = 0;
////    typedef std::pair<int,int> ColRowPair;
//    BOOST_FOREACH(ColRowPair p, currSeg->points) { // for each point in segment
//      int col = p.first;
//      int row = p.second;
//      if (visibilityResult.get(col,row))
//        ++nbVisible;
//    }
//    currSeg->ratioVisibleNextFrame = (double)nbVisible / (double)currSeg->ptCount;
//    currSeg->ptCountVisibleNextFrame = nbVisible;
//  }
}

void visibilityPreviousFrame(FrameSPtr cf, FrameSPtr pf, LidarImageProjector &projector)
{
  LidarImage<bool> &visibilityResult = cf->visibleInPrevFrame;
  visibilityOtherFrame(cf, pf, projector, visibilityResult);
}

void visibilityOtherFrame(FrameSPtr cf, FrameSPtr of, LidarImageProjector &projector, LidarImage<bool> &visibilityResult)
{
  // idea: make 2 passes over the image:
  // 1) project points of current frame to other frame, calculating pixel-wise min-distance
  // 2) deciding per pixel about its visibility by comparing projected distance to min-distance

  int hsize, vsize;
  cf->distance.getSize(hsize,vsize); // should be representative to take resolution information out of distance image
  ParameterHeap* params = ParameterHeap::get();
  double distThresh = params->projectionDistTolerance;
  LidarImage<int>   &projColNextFrame = cf->tmpInt1; // stores the column of associated pixel in next frame
  LidarImage<int>   &projRowNextFrame = cf->tmpInt2; // stores the row of associated pixel in next frame
  LidarImage<double> &projDistNextFrame = cf->tmpFloat1; // stores the distance relative to next frame
  LidarImage<double> &minDistNextFrame = cf->tmpFloat2; // pixel-wise minimum value for next frame, to decide about visibility
  
  // 1st pass over image: project 3D coordinates to new frame and calculate min-dist
  minDistNextFrame = of->distance; // assign distance values of next frame
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      int hi, vi; 
      double dist;
      DVector p = cf->point3D.get(col, row);
      p.resize(4);
      p(3)=1.0;
      p = prod(cf->positionHTM2w,p); // transform to world coordinates
      p = prod(of->positionHTM2v,p); // transform to vehicle coordinates of other frame
      if (projector.getImageIndexRel(p(0), p(1), p(2), hi, vi, dist)) {
        projColNextFrame.set(col,row,hi);
        projRowNextFrame.set(col,row,vi);
        projDistNextFrame.set(col,row,dist);
        if (dist < minDistNextFrame.get(hi,vi))
          minDistNextFrame.set(hi,vi,dist);
      } else {
        projDistNextFrame.set(col,row,DBL_MAX);
      }
    } // end: for each column
  } // end: for each row

  // 2nd pass over image: decide about visibility
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double dstPNF = projDistNextFrame.get(col,row);
      if (dstPNF == DBL_MAX)
        visibilityResult.set(col,row,false);
      else {
        int colPF = projColNextFrame.get(col,row);
        int rowPF = projRowNextFrame.get(col,row);
        double dstNF = of->distance.get(colPF,rowPF);
        double normalZ = fabs(of->normal3D.get(colPF,rowPF)(2));
        double thresh = distThresh * (1+0.2*normalZ*normalZ*dstPNF);
        // TODO (6): use correct plane-projection
//          DVector& pixPt = ns.frame->point3D.get(colPC,rowPC);
//          DVector& pixNr = ns.frame->normal3D.get(colPC,rowPC);
//          double pixPlaneDist = inner_prod(ns.searchPoint-pixPt, pixNr);
//          if (pixPlaneDist < -ns.projDistTolerance)
        if ((dstNF == DBL_MAX) || ((dstPNF - minDistNextFrame.get(colPF,rowPF)) > thresh)) {
          visibilityResult.set(col,row,false);
        } else {
          visibilityResult.set(col,row,true);
        }
      }
    } // end: for each column
  } // end: for each row
}


