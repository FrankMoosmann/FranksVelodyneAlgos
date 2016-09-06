#include "AlgoSegmentation.hpp"

#include <cfloat>
#include <cmath>
#include <climits>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "LidarImageFeatures.hpp"
#include "AlgoFrameFeatures.hpp"
#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;

// returns a "connection score" (based on local convexity). Together with "connectivity" this is compared against a threshold to segment the image
double connectWith(FrameSPtr frame, int col, int row, int col2, int row2, ParameterHeap* params,
                   double convexThresh, double twistThresh, double normalSimThresh)
{
  double       &convexFac     = params->segConvFac;
  double       &twistFac      = params->segTwistFac;
  bool         &useNormConf   = params->segUseNormConf;
  bool         &segNcdSquare  = params->segNcdSquare;
  bool         &crit1Active   = params->segCrit1Active;
  bool         &crit2Active   = params->segCrit2Active;
  bool         &crit3Active   = params->segCrit3Active;
  unsigned int &variantTilt = params->segVariantTilt;
  unsigned int &variantTwist = params->segVariantTwist;

  // prepare variables:
  DVector &n1 = frame->normal3D.get(col,row);
  DVector &n2 = frame->normal3D.get(col2,row2);
  double n1x = n1(0);
  double n1y = n1(1);
  double n1z = n1(2);
  double n2x = n2(0);
  double n2y = n2(1);
  double n2z = n2(2);

  if (n1x+n1y+n1z == 0.0) return 0.0; // length 0 meaning normal vector does not exist -> always reject
  if (n2x+n2y+n2z == 0.0) return 0.0; // length 0 meaning normal vector does not exist -> always reject

  DVector &p1 = frame->point3D.get(col,row);
  DVector &p2 = frame->point3D.get(col2,row2);
  double d12x = p2(0) - p1(0);
  double d12y = p2(1) - p1(1);
  double d12z = p2(2) - p1(2);
  double d12l = length(d12x, d12y, d12z);

  // calculate connection cost:
  double retVal = 0.0; 
  
  // both surfaces approx. vertical
  if (crit1Active) {
    retVal = max(retVal,  min(1.0f-fabs(n1z),1.0f-fabs(n2z)) ); // = 0.5 at 30°
  }

  // similar normal direction
  if (crit2Active) {
    double dp2n = dotProduct(n1x, n1y, n1z, n2x, n2y, n2z); // 1 if the same, 0 if perpendicular, -1 if opposite
    retVal = max(retVal, LidarImageFeatures::sigmoidLikeSoftThresh(1-dp2n, normalSimThresh, 10.0) ); // currently 0.5 @ 10deg
  }

  // convexity
  if (crit3Active) {
    double c1x, c1y, c1z, c2x, c2y, c2z;
    crossProduct(n1x, n1y, n1z, d12x, d12y, d12z, c1x, c1y, c1z);
    crossProduct(n2x, n2y, n2z, -d12x, -d12y, -d12z, c2x, c2y, c2z);
    normalize(c1x, c1y, c1z); // cross products must even be normalized if both input vectors are normalized: |a x b|=|a||b|sin(alpha)
    normalize(c2x, c2y, c2z);
    double cosN1D = dotProduct(n1x, n1y, n1z,  d12x,  d12y,  d12z)/d12l; // n is normalized to length 1, but d12 is not
    double cosN2D = dotProduct(n2x, n2y, n2z, -d12x, -d12y, -d12z)/d12l;

    double convVal = DBL_MAX;

    // check tilting angle
    switch (variantTilt) {
    case 1:
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN1D, convexThresh, convexFac));
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN2D, convexThresh, convexFac));
// the below variant should be equivalent, since sigmoidLikeSoftThresh is monotonically falling
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(max(cosN1D,cosN2D), convexThresh, convexFac));
      break;

    case 2:
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN1D/d12l, convexThresh, convexFac));
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN2D/d12l, convexThresh, convexFac));
// the below variant should be equivalent, since sigmoidLikeSoftThresh is monotonically falling
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(max(cosN1D,cosN2D)/d12l, convexThresh, convexFac));
      break;

    case 3:
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN1D*d12l, convexThresh*d12l, convexFac));
//      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(cosN2D*d12l, convexThresh*d12l, convexFac));
      // the below variant should be equivalent, since sigmoidLikeSoftThresh is monotonically falling
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(max(cosN1D,cosN2D)*d12l, convexThresh*d12l, convexFac));
      break;
      
    case 409: // [IV09] version of segmentation (hard thresholds). needs crit3Active=true, variantTilt=409, variantTwist=0 and useNormConf=false
              // IV settings: params->segConvThresh=6, params->segNormThresh=2, params->segTwistThres=60
//    double convexThresh = cos((90.0f - params->segConvThresh)*M_PI/180); // convert to RAD
//    double normalSimThresh = 1-cos(params->segNormThresh*M_PI/180); // convert to RAD
      double nSimOrigRAD = acos(1-normalSimThresh);
      double nSimIV = 1 - d12l*cos(M_PI/2 - nSimOrigRAD);
//    double twistThresh = cos(params->segTwistThresh*M_PI/180); // convert to RAD
      if ( (   (cosN1D <= convexThresh)
            && (cosN2D <= convexThresh))
        || (dotProduct(n1x, n1y, n1z, n2x, n2y, n2z) >= nSimIV)
        || (   (fabs(n1z) < twistThresh) // misuse of this variable
            && (fabs(n2z) < twistThresh))
         )
        convVal = 1.0;
      else
        convVal = 0.0;
      retVal = convVal;
      break;
    }
  
    // check twisting angle
    switch (variantTwist) {
    case 1:
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(1.0-abs(dotProduct(c1x, c1y, c1z, c2x, c2y, c2z)), twistThresh, twistFac));
      break;

    case 2:
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(abs(dotProduct(c2x, c2y, c2z, n1x, n1y, n1z)), twistThresh, twistFac));
      convVal = min(convVal, LidarImageFeatures::sigmoidLikeSoftThresh(abs(dotProduct(c1x, c1y, c1z, n2x, n2y, n2z)), twistThresh, twistFac));
      break;
    }

    retVal = max(retVal,  convVal);
  
    // up to now, the return value is based upon the normal directions
    // this is only meaningful, if the normal confidence is high
    // --> cope for the case of bad normals, and blend the cases linearly
    if (useNormConf) {
      double nc1 = frame->normalConfidence.get(col,row);
      double nc2 = frame->normalConfidence.get(col2,row2);
    //  retVal =  ((  nc1)*(  nc2)) * retVal // both normals good -> use above calculated value
    //           +((1-nc1)*(  nc2)) * 0.0  // first normal bad, second good -> don't connect
    //           +((  nc1)*(1-nc2)) * 0.0  // first normal good, second bad -> don't connect
    //           +((1-nc1)*(1-nc2)) * 1.0; // both normals bad -> connect
      double ncm = (nc1+nc2)/2.0; // mean, [0..1]
      double ncd = fabs(nc1-nc2); // difference, [0..1]
      if (segNcdSquare)
        ncd = ncd*ncd;
      retVal = (1-ncd)*( ncm*retVal + (1-ncm)*1.0);
    }
  }

  return retVal;
};

// returns a value between 0 (for x<=0) and 1 (for x>=xMax). Inbetween linear relation
inline double linSat(double x, double xMax)
{
  if (x<=0.0f)
    return 0.0f;
  if (x>=xMax)
    return 1.0f;
  return x/xMax;
}

// modifies threshold [0.0..1.0]->[0.0..1.0] depending on projected tracks
void adjustSegThresh(FrameSPtr frame, int col, int row, int col2, int row2, double &thresh, const double threshNorm)
{
  LidarFrame::TrackIndex ti1 = frame->extProjTracks.get(col,row);
  LidarFrame::TrackIndex ti2 = frame->extProjTracks.get(col2,row2);
  if ((ti1 != UINT_MAX) && (ti2 != UINT_MAX) && (frame->projSttIdx == 0)) {
    PointCloudTrack *track1 = frame->shortTermTracks[0].tracks[ti1].get();
    PointCloudTrack *track2 = frame->shortTermTracks[0].tracks[ti2].get();
    if (track1 == track2) { // same track
      thresh = (thresh + 0.0f)*0.8f; // move towards 0 // TODO (8): decrease threshold in dependence of position uncertainty/age of track
    } else { // different tracks
      double pct1 = (double)(track1->getPointCount());
      double pct2 = (double)(track2->getPointCount());
      DVector t1,t2;
      track1->getVelocity(t1);
      track2->getVelocity(t2);
      t1 *= 0.1; t2 *= 0.1; // turn velocity into translation vector
      double nt1 = norm_2(t1);
      double nt2 = norm_2(t2);
      double mnt1t2 = max(nt1,nt2);
      // merge segments with similar translation vector t >> 0, keep separate if t1,t2 >> 0 are not similar
      double thresh2 = fabs(mnt1t2)< 0.00001 ? 0.5f : 0.5*norm_2(t1-t2)/mnt1t2;
      double weight2 = linSat(mnt1t2, threshNorm);
      // TODO (8): use uncertainty of translation estimation
      // keep large segments separate, merge only smaller ones
      double segSizeWeight = 0.1;
      thresh2 = (linSat(min(pct1,pct2), 50.0)*segSizeWeight + thresh2*weight2) / (segSizeWeight+weight2);
      weight2 = max(weight2, segSizeWeight);
      // now adjust segmentation threshold:
      thresh = thresh*(1-weight2) + thresh2*weight2; // linear combination
//      thresh = 1.0;
//      thresh = sqrt(thresh); // move towards 1
    }
  }
};


double getSegAvgDist(const LidarSegment &seg, FrameSPtr frame)
{
  LidarSegment::ColRowConstIterator cr,crEnd;
  unsigned int pixCnt;
  seg.getColRows(cr, crEnd, pixCnt);
  double d = 0.0;
  while (cr != crEnd) {
    d += frame->distance.get(*cr);
    ++cr;
  }
  d /= (double)pixCnt;
  //cout << " d=" << d << flush;
  return d;
}

void segment(FrameSPtr frame, bool verbose)
{
  ParameterHeap* params = ParameterHeap::get();
  double convexThresh = cos((90.0f - params->segConvThresh)*M_PI/180); // convert to RAD
  double twistThresh = cos(params->segTwistThresh*M_PI/180); // convert to RAD
  double normalSimThresh = 1-cos(params->segNormThresh*M_PI/180); // convert to RAD
  double minSegScore = params->segMinScore;
  double linSatVelWeightThresh = params->segLinSatVelWeightThresh;
  unsigned int minSegmentPtCnt = params->segMinSizePt;
  // horizPtCnt = sideLength/(d*tan(samplingAngle))
  // ptCnt = (sideLength^2 / tan(samplingHAngle) / tan(samplingVAngle)) / d^2
  double minSegSizeFac = params->segMinSizeM*params->segMinSizeM/tan(params->lidarHAngResolRAD)/tan(params->lidarVAngResolRAD);
  bool useTrackIDs = params->segUseTrackIDs;
  bool useMultiscale = params->segUseMultiscale;
  bool useTriangles = params->segUseTriangles;
  int hsize, vsize;
  frame->distance.getSize(hsize,vsize);
  
  if (verbose) cout << endl << "calculating segments..." << flush;

  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      frame->segments(col,row).reset(); // reset segment
      frame->segments(col,row).setValid((frame->distance.get(col,row) != DBL_MAX));
    }
  }

  if (verbose) cout << "." << flush;
  // 1) first pass over image: calculate segmentation criterion
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      frame->segmentCritR(col,row)  = (col >= hsize-1)                       ? 0.0 : connectWith(frame, col+1, row  , col, row, params, convexThresh, twistThresh, normalSimThresh);
      frame->segmentCritD(col,row)  = (row >= vsize-1)                       ? 0.0 : connectWith(frame, col  , row+1, col, row, params, convexThresh, twistThresh, normalSimThresh);
      frame->segmentCritRU(col,row) = ((col >= hsize-1) || (row <= 0))       ? 0.0 : connectWith(frame, col+1, row-1, col, row, params, convexThresh, twistThresh, normalSimThresh);
      frame->segmentCritRD(col,row) = ((col >= hsize-1) || (row >= vsize-1)) ? 0.0 : connectWith(frame, col+1, row+1, col, row, params, convexThresh, twistThresh, normalSimThresh);
      if (useMultiscale) {
        if ((col > 0) && (col < hsize-2))
          frame->segmentCritR(col,row) = min(frame->segmentCritR(col,row),
                                             max(connectWith(frame, col-1, row, col+1, row, params, convexThresh, twistThresh, normalSimThresh),
                                                 connectWith(frame, col,   row, col+2, row, params, convexThresh, twistThresh, normalSimThresh)));
        if ((row > 0) && (row < vsize-2))
          frame->segmentCritD(col,row) = min(frame->segmentCritD(col,row),
                                             max(connectWith(frame, col, row-1, col, row+1, params, convexThresh, twistThresh, normalSimThresh),
                                                 connectWith(frame, col, row,   col, row+2, params, convexThresh, twistThresh, normalSimThresh)));
      }
    }
  }
  // 1a) connecting nodes to their upper and left neighbors depending on the segmentation criterion
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double ch = frame->segmentCritR.get(col, row);
      double cv = frame->segmentCritD.get(col, row);
      double h,v;
      if (!useTriangles) {
        // simple version:
        h = 1.0;
        v = 1.0;
      } else {
        // extended version taking into account 4 triangles:
        // each triangle-score is min(side1, side2, side3),
        // but all 4 triangles share one side (e.g. ch), min of this is moved to the end to save time
        // connection right:
        h = 0.0;
        h = max(h, min(frame->segmentCritD.get(col, row-1),frame->segmentCritRD.get(col, row-1))); // upper left triangle
        h = max(h, min(frame->segmentCritD.get(col+1, row-1),frame->segmentCritRU.get(col, row))); // upper right triangle
        h = max(h, min(frame->segmentCritD.get(col+1, row),frame->segmentCritRD.get(col, row))); // lower right triangle
        h = max(h, min(cv,frame->segmentCritRU.get(col, row+1))); // lower left triangle
        // connection down:
        v = 0.0;
        v = max(v, min(ch,frame->segmentCritRU.get(col, row+1))); // upper right triangle
        v = max(v, min(frame->segmentCritR.get(col, row+1),frame->segmentCritRD.get(col, row))); // lower right triangle
        v = max(v, min(frame->segmentCritR.get(col-1, row),frame->segmentCritRD.get(col-1, row))); // upper left triangle
        v = max(v, min(frame->segmentCritR.get(col-1, row+1),frame->segmentCritRU.get(col-1, row+1))); // lower left triangle
      }
      frame->segmentConnectH(col,row) = min(h, ch);
      frame->segmentConnectV(col,row) = min(v, cv);
    }
  }

  // this HACK replaces the following 2) method and creates a single segment
//  LidarSegment *onlySeg = NULL;
//  for (int row = 0; row < vsize; ++row) {
//    for (int col = 0; col < hsize; ++col) {
//      if (frame->distance.get(col,row) == DBL_MAX) continue;
//      LidarSegment *s = frame->segments.getAdr(col,row);
//      if (onlySeg == NULL) {
//        onlySeg = s;
//      } else {
//        onlySeg->merge(s);
//      }
//      onlySeg->addPix(ColRowPair(col,row));
//    }
//  }


  // 2) second pass over image: assign pixel-wise segments
  //    This effectively builds an "inverse" segmentation-forest with each pixel being a root node. Segments may link to others, if they were merged.
  //    The leafs effectively represent the segments, all root nodes / pixels that trace down to the same leaf belong to the same segment
  if (verbose) cout << "." << flush;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      LidarSegment *s = frame->segments.getAdr(col,row);
      double minSegScoreL = minSegScore;
      double minSegScoreU = minSegScore;
      if (useTrackIDs) { // adjust minScore depending on projected tracks
        adjustSegThresh(frame, col, row, col-1, row, minSegScoreL, linSatVelWeightThresh);
        adjustSegThresh(frame, col, row, col, row-1, minSegScoreU, linSatVelWeightThresh);
      }
//      double nc  = frame->normalConfidence.get(col,row);
//      double ncL = frame->normalConfidence.get(col-1,row);
//      double ncU = frame->normalConfidence.get(col,row-1);
//      double ncmL = (nc+ncL)/2.0; // mean, [0..1]
//      double ncmU = (nc+ncU)/2.0;
      // using normalConfidence here seems to degrade results :(
      double ncmL = 1.0;
      double ncmU = 1.0;
      double segScoreL = ((1-ncmL)*1.0 + ncmL*frame->connectivityH.get(col-1,row))*frame->segmentConnectH.get(col-1,row);
      double segScoreU = ((1-ncmU)*1.0 + ncmU*frame->connectivityV.get(col,row-1))*frame->segmentConnectV.get(col,row-1);
      bool connectLeft = (segScoreL > minSegScoreL);
      bool connectUpper = (segScoreU > minSegScoreU);
      // store for visualization:
      frame->segmentThreshH.set(col-1,row,minSegScoreL);
      frame->segmentThreshV.set(col,row-1,minSegScoreU);
      frame->segmentGrowDecisionH.set(col-1,row,connectLeft);
      frame->segmentGrowDecisionV.set(col,row-1,connectUpper);
      LidarSegment *sl = frame->segments.get(col-1,row).getPtr();
      LidarSegment *su = frame->segments.get(col,row-1).getPtr();
      if (connectLeft) {
        sl->merge(s);
      } else if (col > 0) {
        LidarSegment::makeNeighbors(s,sl,segScoreL);
      }
      if (connectUpper) {
        su->merge(s);
      } else if (row > 0) {
        LidarSegment::makeNeighbors(s,su,segScoreU);
      }
      s->addPix(ColRowPair(col,row));
    }
  } // now each pixel has a segment assigned


  // 3) third pass over image: trace all pointers down to leafs and cleanup neighbor-list.
  if (verbose) cout << "." << flush;
  unsigned int segCount = 0;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      LidarSegment &seg = frame->segments.get(col,row);
      if (seg.isProxy()) {
        seg.update();
      } else {
        if (   (seg.getPxCount() >= minSegmentPtCnt)
            && (seg.getPxCount() >= minSegSizeFac/max(1e-6,pow(getSegAvgDist(seg,frame),2))) ) {
          seg.setValid(true);
          ++segCount;
        } else {
          seg.setValid(false);
        }
      } // end: else: segment is leaf, i.e. not proxy
    } // end: loop over cols
  } // end: loop over rows

  // 4) check consistency
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      LidarSegment &s = frame->segments.get(col,row);
      if (s.isValid()) {
        if (frame->distance.get(col,row) == DBL_MAX)
          if (verbose) cout << endl << "ERROR: pixel " << col << "," << row << " is invalid but its segment-object valid!";
      }
    }
  }
  if (verbose) cout << "done. created " << segCount << " segments." << flush;
}

void recolorSegments(FrameSPtr frame)
{
  int hsize, vsize;
  frame->distance.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      LidarSegment &seg = frame->segments.get(col,row);
      if (!seg.isProxy()) {
        seg.setColor(0.2 + 0.8*(double)rand()/(double)RAND_MAX,
                     0.2 + 0.8*(double)rand()/(double)RAND_MAX,
                     0.2 + 0.8*(double)rand()/(double)RAND_MAX);
      }
    }
  }
}

// functions for creating GraphML-File from segment-segment associations of two frames
void keySpecifier(std::ofstream &of)
{
  of << "  <key id=\"d3\" for=\"node\" yfiles.type=\"nodegraphics\"/>" << endl;
  of << "  <key id=\"d5\" for=\"node\" attr.name=\"description\" attr.type=\"string\"/>" << endl;
  of << "  <key id=\"d6\" for=\"edge\" yfiles.type=\"edgegraphics\"/>" << endl;
  of << "  <key id=\"d8\" for=\"edge\" attr.name=\"description\" attr.type=\"string\"/>" << endl;
}
void segmentVisitor_base(const LidarSegment* seg, std::ofstream &of, float rowOffset)
{
  float colC, rowC;
  seg->getCenter(colC, rowC);
  float size = (sqrt((double)seg->getPxCount()/M_PI)); // area -> radius -> downscale
  int r = (int)(seg->getColorR()*255.0);
  int g = (int)(seg->getColorG()*255.0);
  int b = (int)(seg->getColorB()*255.0);
  of << "      <data key=\"d5\"/>" << endl;
  of << "      <data key=\"d3\">" << endl;
  of << "        <y:ShapeNode>" << endl;
  of << "          <y:Geometry height=\"" << size << "\" width=\"" << size << "\" x=\"" << colC*3 << "\" y=\"" << rowC+rowOffset << "\"/>" << endl;
  of << "          <y:Fill color=\"#" << boost::format("%|02X|%|02X|%|02X|") % r % g % b << "\" transparent=\"false\"/>" << endl;
  of << "          <y:BorderStyle color=\"#000000\" type=\"line\" width=\"0.5\"/>" << endl;
//  of << "          <y:NodeLabel alignment=\"center\" autoSizePolicy=\"content\" fontFamily=\"Dialog\" fontSize=\"12\" fontStyle=\"plain\" hasBackgroundColor=\"false\" hasLineColor=\"false\" height=\"17.96875\" modelName=\"internal\" modelPosition=\"c\" textColor=\"#000000\" visible=\"true\" width=\"172.533203125\" x=\"-71.2666015625\" y=\"6.015625\">A1</y:NodeLabel>" << endl;
  of << "          <y:Shape type=\"circle\"/>" << endl;
  of << "        </y:ShapeNode>" << endl;
  of << "      </data>" << endl;
}
void segmentVisitorL1(const LidarSegment* seg, std::ofstream &of)
{
  segmentVisitor_base(seg, of, 0.0);
}
void segmentVisitorL2(const LidarSegment* seg, std::ofstream &of)
{
  segmentVisitor_base(seg, of, 100.0);
}
void ConnVisitor(const LidarSegment* s1, const LidarSegment* s2, unsigned int val, std::ofstream &of)
{
  (void)s1;
  (void)s2;
  int v = 200 - min(200,(int)(val*2)); // 100 --> black 0 --> light gray
  of << "      <data key=\"d8\"/>" << endl;
  of << "      <data key=\"d6\">" << endl;
  of << "        <y:PolyLineEdge>" << endl;
  of << "          <y:Path sx=\"0.0\" sy=\"0.0\" tx=\"0.0\" ty=\"0.0\"/>" << endl;
  of << "          <y:LineStyle color=\"#" << boost::format("%|02X|%|02X|%|02X|") % v % v % v << "\" type=\"line\" width=\"" << min(7,(int)(sqrt((double)val/2.0))) << ".0\"/>" << endl;
  of << "          <y:EdgeLabel alignment=\"center\" distance=\"2.0\" fontFamily=\"Dialog\" fontSize=\"12\" fontStyle=\"plain\" hasBackgroundColor=\"false\" hasLineColor=\"false\" height=\"17.96875\" modelName=\"six_pos\" modelPosition=\"tail\" preferredPlacement=\"anywhere\" ratio=\"0.5\" textColor=\"#000000\" visible=\"true\" width=\"11.634765625\" >" << val << "</y:EdgeLabel>" << endl;
  //x=\"2.0114933270290436\" y=\"38.499270405273435\"
  of << "          <y:Arrows source=\"none\" target=\"none\"/>" << endl;
  of << "          <y:BendStyle smoothed=\"false\"/>" << endl;
  of << "        </y:PolyLineEdge>" << endl;
  of << "      </data>" << endl;
}
SegmentationLinks linkSegments(const LidarImage<LidarSegment> &segImage1, const LidarImage<LidarSegment> &segImage2)
{
//  cerr << endl << "linking segments..." << flush;
  SegmentationLinks segLinks;
  int hsize, vsize;
  segImage1.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
//    cerr << "." << flush;
    for (int col = 0; col < hsize; ++col) {
      const LidarSegment &seg1 = segImage1.get(col,row);
      const LidarSegment &seg2 = segImage2.get(col,row);
      if ((seg1.isValid()) && (seg1.getPxCount() > 0) && (seg2.isValid()) && (seg2.getPxCount() > 0)) {
        segLinks.addConnection(seg1.getPtr(), seg2.getPtr());
      } else {
        if (seg1.isValid() && (seg1.getPxCount() > 0)) segLinks.addNode1(seg1.getPtr());
        if (seg2.isValid() && (seg2.getPxCount() > 0)) segLinks.addNode2(seg2.getPtr());
      }
    }
  }
//  cerr << "done" << flush;
  return segLinks;
}
void saveLinks(const SegmentationLinks &links, string filename)
{
  cerr << "saving segment-links to file..." << flush;
  links.saveToGraphML(filename, boost::bind(&keySpecifier, _1)
                              , boost::bind(&segmentVisitorL1, _1, _2)
                              , boost::bind(&segmentVisitorL2, _1, _2)
                              , boost::bind(&ConnVisitor, _1, _2, _3, _4)
  );
  cerr << "done" << flush;
}
void linkAndSave(FrameSPtr frame1, FrameSPtr frame2)
{
  SegmentationLinks sl = linkSegments(frame1->segments, frame2->segments);
  saveLinks(sl, "segLinks.graphml");
}

