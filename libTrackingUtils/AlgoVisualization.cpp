#include "AlgoVisualization.hpp"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <QColor>
#include <GL/glut.h>
#include <Gui3DQt/graphics.hpp>
#include <Gui3DQt/passatmodel.hpp>
#include "CovarEllipsoidRendering.hpp"
#include "MatrixDefs.hpp"
#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;
using namespace matrixTools::ublas;

AlgoVisualization::AlgoVisualization()
{
  Gui3DQt::PassatModel::setColor(0, 0.5, 0.5, 0.5);
}

AlgoVisualization::~AlgoVisualization()
{
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    2D Rendering Methods   /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void AlgoVisualization::renderMono(LidarImage<bool> &li, QImage &qi, bool invert)
{
  LidarImageVisualization::renderMono(li,qi,invert);
}

void AlgoVisualization::renderColored(LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs, double valInvalid, LidarImageVisualization::ColorMode2D mode)
{
  LidarImageVisualization::renderColored(li, qi, minVal, maxVal, useAbs, valInvalid, mode);
}

void AlgoVisualization::renderGrey(LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs, bool invert)
{
  LidarImageVisualization::renderGrey(li, qi, minVal, maxVal, useAbs, invert);
}

void AlgoVisualization::renderColoredNormals(LidarImage<DVector> &lin, LidarImage<double> &conf, QImage &qi)
{
  LidarImageVisualization::renderColoredNormals(lin,qi,&conf);
}

void AlgoVisualization::renderSegments(LidarImage<LidarSegment> &li, QImage &qi)
{
  int hsize, vsize;
  li.getSize(hsize,vsize);
  // render pixel in segment-colors
  qi = QImage(hsize,vsize,QImage::Format_RGB16); //Format_Mono);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      LidarSegment &s = li.get(col,row);
      if (s.isValid()) {
        qi.setPixel( col, row, QColor(s.getColorR()*255, s.getColorG()*255, s.getColorB()*255).rgb() );
      } else {
        qi.setPixel( col, row, QColor(0,0,0).rgb() );
      }
    }
  }
}

void AlgoVisualization::renderTrackProjection(const LidarImage<LidarFrame::TrackIndex> &trackProj, const LidarFrame::TrackSPtrVec &tracks, QImage &qi)
{
  int hsize, vsize;
  trackProj.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  qi.fill(QColor(0,0,0).rgb());
  // paint projected tracks
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      LidarFrame::TrackIndex tIdx = trackProj.get(col,row);
      if (tIdx != UINT_MAX) {
        const PointCloudTrack::SPtr t = tracks[tIdx];
        uint color = QColor(t->getColorR()*255, t->getColorG()*255, t->getColorB()*255).rgb();
        qi.setPixel( col, row, color );
      }
    }
  }
}

void AlgoVisualization::renderTrackProjectionCnt(LidarImage<unsigned int> &tcnt, QImage &qi)
{
  int hsize, vsize;
  tcnt.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  qi.fill(QColor(0,0,0).rgb());
  // paint projected tracks
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      unsigned int cnt = tcnt.get(col,row);
      if (cnt == 1) {
        qi.setPixel( col, row, QColor(10, 255, 10).rgb() );
      } else if (cnt >= 2) {
        qi.setPixel( col, row, QColor(255, 10, 10).rgb() );
      }
    }
  }
}

void AlgoVisualization::renderTranslationNorm(FrameSPtr frame, QImage &qi)
{
  if (frame->projSttIdx > frame->sttValidationSteps) return;
  int hsize, vsize;
  frame->distance.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  qi.fill(QColor(0,0,0).rgb());
  int hue,sat,val;
  const double maxMPS(3.0f); // speed in m/s that is displayed with full value
//  const double meanSpeedCovar(2.0f); // speed in (m/s)^2 that is displayed with half intensity
  double egoSpeed = frame->speed;
  // paint speed of projected tracks
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      LidarFrame::TrackIndex trackIdx = frame->projTracks.get(col,row);
      if (trackIdx != UINT_MAX) {
        PointCloudTrack::SPtr track = frame->shortTermTracks[frame->projSttIdx].tracks[trackIdx];
        DVector tSpeedXYZ; track->getVelocity(tSpeedXYZ); // returns a 3-dimensional vector: vx,vy,vz
        tSpeedXYZ(0) = tSpeedXYZ(0) + egoSpeed; // compensate ego-speed (which is always in local x-direction)
        DMatrix tCovar; track->getCovar(tCovar); // get complete 12x12 covariance matrix
        DMatrix tCovarSpeedXYZ = ublas::matrix_range<DMatrix>(tCovar, ublas::range(9,12), ublas::range(9,12));
        DSMatrix tCovarSpeedXYZSymm = ublas::symmetric_adaptor<DMatrix, ublas::lower>(tCovarSpeedXYZ);
        DCMatrix invCovar = invSym(tCovarSpeedXYZSymm);
        DVector tmp = ublas::prod(invCovar,tSpeedXYZ);
        double mahalDistSq = ublas::inner_prod(tSpeedXYZ,tmp); // squared mahal. dist is chi-square distributed
        double weight = exp(-0.5*mahalDistSq); // thus a square-distance of 6 (very low probability) will get weight=0.05
        double speed_magnitude = norm_2(tSpeedXYZ); // TODO (8): use mahalanobis distance instead!!!
//        double plane_angle_RAD = atan2(tSpeed(1),tSpeed(0));
//        if (plane_angle_RAD < 0) plane_angle_RAD += M_PI*2; //normalize to 0..2PI
//        double speed_covar = (tCovar(9,9) + tCovar(10,10) + tCovar(11,11))/3; // average covariance
//        hue = plane_angle_RAD/M_PI*180.0; //scale to 0..359
//        sat = 255 * fmin(1.0f,speed_magnitude/maxMPS); // scale to 0...255, 
//        val = 255 * exp(- speed_covar/meanSpeedCovar*0.69); // exp(-0.69)=0.5
        hue = 180 - 180 * min(weight*speed_magnitude,maxMPS)/maxMPS; //scale to 0..180
        sat = 255 * min(weight*speed_magnitude,maxMPS)/maxMPS; // scale to 0...255,
        val = 255 * weight; //exp(- speed_covar/meanSpeedCovar*0.69); // exp(-0.69)=0.5
//        cout << " h:" << hue << " s:" << sat << " v:" << val << flush;
        qi.setPixel( col, row, QColor::fromHsv(hue, 255, val).rgb() );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    3D Rendering Methods   /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void AlgoVisualization::render3DCarPos(bool renderCS)
{
  if (renderCS) {
  // render own coordinate system
    glPushMatrix();
    glTranslatef(0.02, 0, -0.3);
    Gui3DQt::Graphics::draw_coordinate_frame(0.8);
    glPopMatrix();
//   glLineWidth(3);
//   glBegin(GL_LINES);
//   glColor3f(1.0,0.0,0.0); // x-axis in red
//   glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(0.5, 0.0, 0.0);
//   glColor3f(0.0,1.0,0.0); // y-axis in green
//   glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(0.0, 0.5, 0.0);
//   glColor3f(0.0,0.0,1.0); // z-axis in blue
//   glVertex3d(0.0, 0.0, 0.0);
//   glVertex3d(0.0, 0.0, 0.5);
//   glEnd();
//   glLineWidth(1);
  }

  glPushMatrix();
//  glRotatef(5.0, 0.0 ,1.0 ,0.0);
  glTranslatef(0.1, 0, -1.0);
//  glTranslatef(1.25, 0, 0.9);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  Gui3DQt::PassatModel::draw(); // uint model_index = 0, double wheel_angle = 0., double velodyne_angle = 0.
  glDisable(GL_LIGHTING);
  glPopMatrix();
}

void AlgoVisualization::render3DGrid()
{
  double center_x = 0.0;
  double center_y = 0.0;
  // ++++++++ draw grid ++++++++++++
  //glEnable(GL_LINE_SMOOTH); // nur bei breiten linien von Bedeutung
  glColor3f(0.4, 0.4, 0.4);
  glLineWidth(0.5);
  glBegin(GL_LINES);
  for (int grid_x = -100; grid_x < 100; grid_x++) {
    glVertex3f(grid_x - center_x, -100 - center_y, 0);
    glVertex3f(grid_x - center_x,  100 - center_y, 0);
  }
  for (int grid_y = -100; grid_y < 100; grid_y++) {
    glVertex3f(-100 - center_x, grid_y - center_y, 0);
    glVertex3f( 100 - center_x, grid_y - center_y, 0);
  }
  glEnd();

  // ++++++++ draw circles ++++++++++++
//  glColor3f(0.6, 0.6, 0.6);
//  glPushMatrix();
//  //glRotatef(kogmo_radians_to_degrees(robot_pose.yaw), 0, 0, 1);
//  //glTranslatef(velodyne_frame->config.laser_offset.x, velodyne_frame->config.laser_offset.y, velodyne_frame->config.laser_offset.z);
//  for(int i = 10; i <= 80; i += 10) {
//    glColor3f(0.5, 0.5, 0.5);
//    glBegin(GL_LINE_LOOP);
//    for(int j = 0; j < 100; j++) {
//      double angle = j / 100.0 * M_PI * 2;
//      glVertex3f(i * cos(angle), i * sin(angle), -1.71);
//    }
//    glEnd();
//  }
//  glPopMatrix();
}

void AlgoVisualization::render3DPoints(FrameSPtr frame, bool renderUncertainty, bool invert)
{
  double r=0.1; double g=0.1; double b=0.1;
  LidarImageVisualization::ColorMode3D mode=LidarImageVisualization::ColorFixed;
  if (!invert)
    mode = LidarImageVisualization::ColorHeight;
  ParameterHeap *params = ParameterHeap::get();
  LidarImageVisualization::render3DPoints(frame->point3D,  mode, params->vis3DPointSize, r, g, b);
  if (renderUncertainty)
    LidarImageVisualization::render3DPointCovar(frame->point3D, frame->pointVariance3D);
}

void AlgoVisualization::render3DNormals(FrameSPtr frame, bool renderUncertainty, bool invert)
{
  (void)invert;
  ParameterHeap *params = ParameterHeap::get();
  LidarImage<double> *normalStdDevRAD = NULL;
  LidarImage<matrixTools::DMatrix> *normalCovar = NULL;
  if (renderUncertainty) {
    normalStdDevRAD = &frame->normalStdDevRAD;
    normalCovar = &frame->normalVariance3D;
  }
  LidarImageVisualization::render3DNormals(frame->point3D, frame->normal3D, &frame->normalConfidence, normalStdDevRAD, normalCovar, params->vis3DPointSize);
}

void AlgoVisualization::render3DConnections(LidarImage<DVector> &liP, LidarImage<double> &connR, LidarImage<double> &connD, LidarImage<double> *connRU, LidarImage<double> *connRD)
{
  LidarImageVisualization::render3DConnections(connR, 0, 0, +1, 0, liP); // hoff1, voff1, hoff2, voff2
  LidarImageVisualization::render3DConnections(connD, 0, 0, 0, +1, liP);
  if (connRD) LidarImageVisualization::render3DConnections(*connRD, 0, 0, +1, +1, liP);
  if (connRU) LidarImageVisualization::render3DConnections(*connRU, 0, 0, +1, -1, liP);
}

void AlgoVisualization::render3DSegmentationLinks(FrameSPtr frame)
{
  glBegin(GL_LINES);
  glColor4f(1.0, 1.0, 1.0, 0.5);
  int hsize, vsize;
  frame->segmentGrowDecisionH.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      if (frame->segmentGrowDecisionH.get(col,row)) {
        DVector &p1 = frame->point3D.get(col, row);
        DVector &p2 = frame->point3D.get(col+1, row);
        if ((p1(0) != DBL_MAX) && (p2(0) != DBL_MAX)) {
          glVertex3f(p1(0), p1(1), p1(2));
          glVertex3f(p2(0), p2(1), p2(2));
        }
      }
    }
  }
  frame->segmentGrowDecisionV.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      if (frame->segmentGrowDecisionV.get(col,row)) {
        DVector &p1 = frame->point3D.get(col, row);
        DVector &p2 = frame->point3D.get(col, row+1);
        if ((p1(0) != DBL_MAX) && (p2(0) != DBL_MAX)) {
          glVertex3f(p1(0), p1(1), p1(2));
          glVertex3f(p2(0), p2(1), p2(2));
        }
      }
    }
  }
  glEnd();
}

void AlgoVisualization::render3DSegments(FrameSPtr frame, bool invert)
{
  (void)invert;
  ParameterHeap *params = ParameterHeap::get();
  glPointSize(params->vis3DPointSize);
  int hsize, vsize;
  frame->segments.getSize(hsize,vsize);
  glBegin(GL_POINTS);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      LidarSegment &s = frame->segments.get(col,row);
      if (s.isValid()) {
        glColor3f(s.getColorR(), s.getColorG(), s.getColorB());
        DVector &p1 = frame->point3D.get(col, row);
        glVertex3f(p1(0), p1(1), p1(2));
      }
    }
  }
  glEnd();
}

void AlgoVisualization::render3DCorrespondences(FrameSPtr frame, FrameSPtr nframe, bool invert)
{
  (void)invert;
  glLineWidth(1);
  glBegin(GL_LINES);
  glColor4f(0.7, 0.7, 0.2, 0.3);
  int hsize, vsize;
  frame->point3D.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
//      if ((onlySegId > 0) && (s->segId != onlySegId)) continue;
      int ncol = frame->matchedColNextFrame.get(col, row);
      int nrow = frame->matchedRowNextFrame.get(col, row);
      if (ncol >= 0) {
        DVector &p1 = frame->point3D.get(col, row);
        double x1 = p1(0);
        double y1 = p1(1);
        double z1 = p1(2);
        DVector &p2 = nframe->point3D.get(ncol, nrow);
        double x2 = p2(0);
        double y2 = p2(1);
        double z2 = p2(2);
        glVertex3f(x1,y1,z1);
        glVertex3f(x2,y2,z2);
      }
    }
  }
  glEnd();
}

mdefs::DVector getCorresp(std::list< std::list<DVector> >::const_iterator oi, unsigned int idx) {
  std::list<DVector>::const_iterator ii = (*oi).begin();
  unsigned int ci = 0;
  while (ci < idx) {
    ++ci; ++ii;
  }
  return *ii;
};

void AlgoVisualization::render3DRegistrationResult(TrackRegistration::SubTrackLst::iterator strackit)
{
  TrackRegistration::SubTrackSPtr subtrack = *strackit; // = struct{uint trackPtCnt; DCMatrix points;} - predicted in EgoCS
  PointCloudTrack::SPtr track = subtrack->track; // points are relative, but state is already updated: registration/translation vector can be derived from speed vector

  // render predicted, subsampled points (which are already in Ego-coordinates)
  glPointSize(5);
  glBegin(GL_POINTS);
  glColor3f(0.3, 0.3, 1.0); // predicted points in blueish
  unsigned int pcount = subtrack->pointsTrackCS.size2();
  for (unsigned int i = 0; i < pcount; ++i) {
    mdefs::DVector pTrans = *DCMatrixColConstIterator(subtrack->pointsTrackCS, i);
    pTrans = mdefs::ublas::prod(subtrack->RInit,pTrans) + subtrack->tInit;
    glVertex3f(pTrans(0), pTrans(1), pTrans(2));
  }
  glEnd();

  /*
  DMatrix R;
  DVector t,v;
  track->getRt2EgoCS(R,t);
  track->getVelocity(v);
  v *= 0.1; // velocity->translation

  PointCloudTrack::RStateIterator  si1       = track->histAllStates.rbegin();
  PointCloudTrack::RStateIterator  si2       = track->histAllStates.rbegin();
  PointCloudTrack::RStateIterator  siend     = track->histAllStates.rend();
  if (si1 == siend) return;
  if (++si2 == siend) return;
  DMatrix tCovar; track->getCovar(tCovar);
  double tSpeedCovar = (tCovar(9,9) + tCovar(10,10) + tCovar(11,11))/3; // average covariance
  
  // render registration as lines between old position and new position
  PointCloudTrack::PointIterator pibegin, piend;
  track->getPoints(pibegin, piend);
  glLineWidth(1);
  glBegin(GL_LINES);
  const double meanSpeedCovar(2.0f); // speed in (m/s)^2 that is displayed with half intensity
  double alpha = exp(- tSpeedCovar/meanSpeedCovar*0.69); // exp(-0.69)=0.5
  glColor4f(alpha*0.2, alpha*1.0, alpha*0.4, alpha); // lines in greenish
  DVector p;
  for (PointCloudTrack::PointIterator pi = pibegin; pi != piend; ++pi) {
    p = pi->s->position; // relative to Track-CS
    p = prod(si1->R,p) + si1->t; // make relative to Ego-CS
    glVertex3f(p(0), p(1), p(2));
    p = pi->s->position; // relative to Track-CS
    p = prod(si2->R,p) + si2->t; // make relative to Ego-CS
    glVertex3f(p(0), p(1), p(2));
  }
  glEnd();
  */

  //cout << " #icpruns=" << subtrack->icpCorrespondenceBuffer.size() << flush;
  unsigned int nbIter = subtrack->icpCorrespondenceBuffer.size();
  if (nbIter == 0) return; // nothing to render
  double inc = 0.8/(double)nbIter;
  unsigned int nbElem = subtrack->icpCorrespondenceBuffer.front().size(); // assume same size of each list
  unsigned int nbCorresp = nbElem/2/pcount;
  //cout << endl << "nbCorresp = " << nbCorresp << ", nbElem = " << nbElem << flush;

  // render icp correspondences as lines (which are already in Ego-coordinates)
  glLineWidth(0.1);
  glBegin(GL_LINES);
  double col = 0.1;
  for (BOOST_AUTO(cl,subtrack->icpCorrespondenceBuffer.begin()); cl!=subtrack->icpCorrespondenceBuffer.end(); ++cl) {
    glColor3f(1-col, col, 0.0); // lines turning from red to green
    assert(cl->size()%2 == 0 && "RenderICPRegistration: uneven point count");
    BOOST_FOREACH(DVector p, *cl) {
      glVertex3f(p(0), p(1), p(2));
    }
    col += inc;
  }
  glEnd();
  // render icp search points (which are already in Ego-coordinates)
  glPointSize(3);
  glBegin(GL_POINTS);
  col = 0.1;
  for (BOOST_AUTO(cl,subtrack->icpCorrespondenceBuffer.begin()); cl!=subtrack->icpCorrespondenceBuffer.end(); ++cl) {
    glColor3f(1-col, col, 0.0); // points turning from red to green
    bool draw = true;
    BOOST_FOREACH(DVector p, *cl) { // draw every second point, i.e. searchPoint
      if (draw) glVertex3f(p(0), p(1), p(2));
      draw = !draw;
    }
    col += inc;
  }
  glEnd();
  // render icp registration trace as continuous line
  glLineWidth(1);
  for (unsigned int ei=0; ei<nbElem; ei+=2*nbCorresp) { // render line for each second point i.e. search point
    glColor3f(0.8, 0.8, 0.8);
    glBegin(GL_LINE_STRIP);
    for (BOOST_AUTO(cl,subtrack->icpCorrespondenceBuffer.begin()); cl!=subtrack->icpCorrespondenceBuffer.end(); ++cl) {
      unsigned int clSize = cl->size();
      if (clSize%2 != 0) {
        cerr << endl << "RenderICPRegistration: uneven point count";
        break;
      }
      if (clSize <= ei) {
        cerr << endl << "RenderICPRegistration: list too small: " << clSize << " but expected " << nbElem;
        break;
      }
      DVector p = getCorresp(cl, ei);
      glVertex3f(p(0), p(1), p(2));
    }
    glEnd();
  }

}

void AlgoVisualization::render3DTrackProps(PointCloudTrack *track)
{
  ParameterHeap *params = ParameterHeap::get();
  // change coordinate system to Track-CS
  glPushMatrix();
  DMatrix htm;
  track->getHTM2EgoCS(htm);
  DCMatrix htm2 = htm; // OpenGL needs column-major
  double *htm2p = htm2.data().begin();
  glMultMatrixd(htm2p);

  // draw point cloud in track-specific color
  double geom = min(1.0,4*track->pointNormalDistributionRatio()); // it is already perfect if 1/4th of the possible-normal-sphere is used
  double r = 1.0-0.8*geom;
  double g = 0.2+0.8*geom;
  double b = 0.5;
  if (track->isMovingObject())
    glColor3f(0.3*r, 0.3*g, 0.3*b);
  glColor3f(r, g, b);
  glPointSize(params->vis3DPointSize+1);
  PointCloudTrack::PointIterator pcurr, pend;
  track->getPoints(pcurr, pend);
  glBegin(GL_POINTS);
  for (; pcurr != pend; ++pcurr) {
    const DVector &v = pcurr->s->position;
    glVertex3f(v(0),v(1),v(2));
  }
  glEnd();

  glPopMatrix();  // change coordinate system back to Ego-CS
}

void AlgoVisualization::render3DTrack(PointCloudTrack *track, bool full, double R, double G, double B, DMatrix *htm2Ego)
{
  if (track == NULL) return;

  ParameterHeap *params = ParameterHeap::get();

  // change coordinate system to Track-CS
  glPushMatrix();
  DMatrix htm;
  if (htm2Ego) htm = *htm2Ego;
  else track->getHTM2EgoCS(htm);
  DCMatrix htm2 = htm; // OpenGL needs column-major
  double *htm2p = htm2.data().begin();
  glMultMatrixd(htm2p);

  // draw coordinate frame
  if (full) {
    const double axisLength = 0.2;
    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0); // x-axis in red
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(axisLength,0.0,0.0);
    glColor3f(0.0,1.0,0.0); // y-axis in reen
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,axisLength,0.0);
    glColor3f(0.0,0.0,1.0); // z-axis in blue
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,0.0,axisLength);
    glEnd();
  }
  // draw point cloud in track-specific color
  double r = R < 0.0 ? track->getColorR() : R;
  double g = G < 0.0 ? track->getColorG() : G;
  double b = B < 0.0 ? track->getColorB() : B;
  if (full)
    glColor3f(r, g, b);
  else
    glColor3f(0.5*r, 0.5*g, 0.5*b);
  PointCloudTrack::PointIterator pcurr, pend;
  // first draw VIPs then others
  track->getPoints(pcurr, pend);
  glPointSize(params->vis3DPointSize+1);
  glBegin(GL_POINTS);
  for (; pcurr != pend; ++pcurr) {
    const DVector &v = pcurr->s->position;
    if (pcurr->s->isVIP)
      glVertex3f(v(0),v(1),v(2));
  }
  glEnd();
  track->getPoints(pcurr, pend);
  glPointSize(params->vis3DPointSize);
  glBegin(GL_POINTS);
  for (; pcurr != pend; ++pcurr) {
    const DVector &v = pcurr->s->position;
    if (!pcurr->s->isVIP)
      glVertex3f(v(0),v(1),v(2));
  }
  glEnd();

  glPopMatrix();  // change coordinate system back to Ego-CS
}

void AlgoVisualization::render3DTrackVariance(PointCloudTrack *track)
{
  // move to track CS (translatory part only)
  glPushMatrix();
  DVector pos;
  track->getPosition(pos);
  glTranslatef(pos(0), pos(1), pos(2));

  // draw variance (which is, like the state, Ego-CS based)
  DMatrix covar; // 12x12 (rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot)
  track->getCovar(covar);

  DMatrixRange velCovar(covar, ublas::range(6,12), ublas::range(6,12));
  glColor4f(0.8,0.2,0.2,0.5); // velocity-variance in red
  render3DCovar66Ellipsoid(velCovar, 0.1, 0.5); // factor 0.1 to turn speed into pos
  DMatrixRange posCovar(covar, ublas::range(0, 6), ublas::range(0, 6));
  glColor4f(0.2,0.8,0.2,0.5); // position-variance in green
  render3DCovar66Ellipsoid(posCovar, 1.0, -0.5);

  // change coordinate system back to Ego-CS
  glPopMatrix();
}

void AlgoVisualization::render3DTrackHistory(PointCloudTrack *track)
{
  if (track == NULL) return;
  const double maxTrans(2.0f); // translation in m that is displayed in red
  ParameterHeap *params = ParameterHeap::get();
  PointCloudTrack::PointIterator pibegin, piend;
  track->getPoints(pibegin, piend);
  PointCloudTrack::RStateIterator sibegin = track->histAllStates.rbegin();
  PointCloudTrack::RStateIterator siend = track->histAllStates.rend();
  if (sibegin == siend) return; // should not happen as track always should have at least the current state in the history queue


  cout << endl << "history:";
  for (PointCloudTrack::RStateIterator sitmp = sibegin; sitmp != siend; ++sitmp) {
    cout << "  (" << sitmp->t << ")";
  }
  cout << flush;

  // render state-trace as lines between track-appearance points
  if (params->visMaxTrackVelVarAng == 0) {
    glLineWidth(2);
    glBegin(GL_LINES); //GL_LINESTRIP
    int toggle = -1;
    for (PointCloudTrack::PointIterator pi = pibegin; pi != piend; ++pi) {
      ++toggle;
      if (toggle%3 != 0) continue; // draw a line for every 3rd point only
      double alpha = 0.8; // let lines vanish by setting alpha*=0.8 in each iteration (thus diminishing exponentially)
      PointCloudTrack::RStateIterator si = sibegin;
      unsigned int stateCnt = 0;
      DVector p = pi->s->position; // relative to Track-CS
      p = prod(si->R,p) + si->t; // make relative to Ego-CS (according to current state)
      ++si; ++stateCnt; // previous state
      while ((si != siend) && (stateCnt <= params->visTrackHistory)) {
        DVector pOld = p;
        p = pi->s->position; // relative to Track-CS
        p = prod(si->R,p) + si->t; // make relative to Ego-CS (according to current old state)
        double dist = ublas::norm_2(pOld-p);
        //cout << dist << "  " << flush;
        int H = 270 - 270 * min(dist,maxTrans)/maxTrans;
        LidarImageVisualization::glColorHSV(H, 255, 255);
        //glColor4f(alpha*0.2, alpha*1.0, alpha*0.4, alpha); // lines in greenish
        glVertex3f(pOld(0), pOld(1), pOld(2));
        glVertex3f(p(0), p(1), p(2));
        alpha *= 0.8;
        ++si; ++stateCnt; // previous state
      }
    }
    glEnd();
  }
  // render probabilistic state-trace as lines between track-appearance points
  if (params->visMaxTrackVelVarAng != 0) {
    glLineWidth(3);
    glBegin(GL_LINES); //GL_LINESTRIP
    DMatrix R,V;
    DVector t,v;
    int toggle = -1;
    for (PointCloudTrack::PointIterator pi = pibegin; pi != piend; ++pi) {
      ++toggle;
      if ((toggle+3)%6 != 0) continue; // draw a line for every 3rd point only
      double alpha = 0.8; // let lines vanish by setting alpha*=0.8 in each iteration (thus diminishing exponentially)
      PointCloudTrack::RStateIterator si = sibegin;
      unsigned int stateCnt = 0;
      const DVector p0TrackCS = pi->s->position; // relative to Track-CS
      DVector state = DVectorRangeConst(si->x, ublas::range(0,6)); // positional-part of state
      PointCloudTrack::state2Rt_T2E(state, R, t);
      DVector pt = prod(R,p0TrackCS) + t; // make relative to Ego-CS

      ++si; ++stateCnt; // previous state
      while ((si != siend) && (stateCnt <= params->visTrackHistory)) {
        glColor4f(alpha*1.0, alpha*0.5, alpha*0.0, alpha); // lines in brownish
        glVertex3f(pt(0), pt(1), pt(2));

        DVector velocity = DVectorRangeConst(si->x, ublas::range(6,12)); // velocity-part of state
        DMatrix tCovarFull = si->P;
        DMatrix tCovarVelocity = ublas::matrix_range<DMatrix>(tCovarFull, ublas::range(6,12), ublas::range(6,12));
        DSMatrix tCovarVelocitySymm = ublas::symmetric_adaptor<DMatrix, ublas::lower>(tCovarVelocity);
        DCMatrix invCovar = invSym(tCovarVelocitySymm);
        DVector tmp = ublas::prod(invCovar,velocity);
        double mahalDistSq = ublas::inner_prod(velocity,tmp); // squared mahal. dist is chi-square distributed
        double weight = exp(-0.3*mahalDistSq); // thus a square-distance of 6 (very low probability) will get weight=0.16
        state = state - params->visMaxTrackVelVarAng * weight * velocity * params->defaultDeltaT; // calculate previous state
        PointCloudTrack::state2Rt_T2E(state, R, t);
        pt = prod(R,p0TrackCS) + t; // make relative to Ego-CS
        glVertex3f(pt(0), pt(1), pt(2));
        alpha *= 0.8;
        ++si; ++stateCnt; // previous state
      }
    }
    glEnd();
  }

}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////        Helper Methods       ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void AlgoVisualization::render3DNormalPlane(const DVector &p, const DVector &n, double radius, GLUquadricObj *quadric)
{
  bool createQuadric = (quadric == NULL);
  if (createQuadric)
    quadric = gluNewQuadric();
  gluQuadricDrawStyle(quadric, GLU_FILL);     // rendering style of the quadric (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
  double length = norm_2(n);
  if ((0.98 < length) && (length < 1.02)) { // only draw if normal vector has approx. length 1
    glPushMatrix();
    glTranslatef(p[0],p[1],p[2]); // move coord sys to point location
    glRotatef(acos(n[2])*180/M_PI, -n[1], n[0], 0.0);//rotate so that normal direction is (0,0,1).
    glColor4f(0.9, 0.9, 0.9, 0.5);
    gluPartialDisk(quadric, 0., radius, 6, 1, 0, 360);//innerRadius, outerRadius, slices, loops, startAngle, sweepAngle
    glPopMatrix();
  }
  if (createQuadric)
    gluDeleteQuadric(quadric);
}
