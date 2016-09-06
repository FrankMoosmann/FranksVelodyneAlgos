#include "Visualizer3DMapper.hpp"

#include <sys/time.h>
#include <cfloat>
#include <iostream>
//#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <QFileDialog>
#include <QInputDialog>

#include <HomogeneousTransformationMatrix.hpp>
#include <kogmo_time.h>
#include <FMUtils.hpp>
#include <PngDistanceImage.hpp>
#include <LidarImageProjector.hpp>
#include <LidarImageProjectorPNG.hpp>
#include <LidarImageFeatures.hpp>
#include <LidarImageVisualization.hpp>

#include <ParameterHeap.hpp>

using namespace std;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////        Visualizer3DMapper        ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

const double VISU_DIST_MAX = 50.0; // in meter
const double VISU_3D_SHIFTUP = 2.0;


Visualizer3DMapper::Visualizer3DMapper(Mapper &mapper_, QWidget *parent)
  : Gui3DQt::Visualizer(parent)
    ,mapper(mapper_)
    ,visAlg()
    ,ptSize(1)
{
  continReader = new QTimer(this);
  continReader->setSingleShot(false);
  connect(continReader, SIGNAL(timeout()), this, SLOT(processSelectedImage()) );

  ui.setupUi(this);
  ui.wFrameControl->setVisible(ui.gbFrameControl->isChecked());
  ui.wFrameProcessing->setVisible(ui.gbFrameProcessing->isChecked());
  ui.wVisualization->setVisible(ui.gbVisualization->isChecked());
  ui.hsImageNb->setMinimum(0);
  ui.hsImageNb->setMaximum(mapper.getImgCount()-1);
  ui.hsImageNb->setValue(0);

  connect(ui.hsImageNb, SIGNAL(valueChanged(int)), this, SLOT(imgSelect(int)) );
  connect(ui.pbProcessImage, SIGNAL(pressed()), this, SLOT(processSelectedImage()) );
  connect(ui.pbProcessNumber, SIGNAL(pressed()), this, SLOT(processCertainNumber()) );
  connect(ui.pbProcessContinuous, SIGNAL(toggled(bool)), this, SLOT(setTimer(bool)) );
  connect(ui.dsbSleep, SIGNAL(valueChanged(double)), this, SLOT(stopTimer()) );
  connect(ui.pbSaveMap, SIGNAL(pressed()), this, SLOT(saveMap()) );
  connect(ui.pbLoadMap, SIGNAL(pressed()), this, SLOT(loadMap()) );
  connect(ui.pbReset, SIGNAL(pressed()), this, SLOT(resetMapping()) );
  connect(ui.pbExport, SIGNAL(pressed()), this, SLOT(exportMap()) );
  connect(ui.pbExport3D, SIGNAL(pressed()), this, SLOT(exportMap3D()) );
  connect(ui.pbImport, SIGNAL(pressed()), this, SLOT(importMap()) );
  connect(ui.pbImport3D, SIGNAL(pressed()), this, SLOT(importMap3D()) );
  connect(ui.pbExportTraj, SIGNAL(pressed()), this, SLOT(exportTraj()) );
  connect(ui.pbImportTraj, SIGNAL(pressed()), this, SLOT(importTraj()) );

  connect(ui.pbReadImage, SIGNAL(pressed()), this, SLOT(readScan()) );
  connect(ui.pbRegister, SIGNAL(pressed()), this, SLOT(registerScan()) );
  connect(ui.pbAddRegisteredScan, SIGNAL(pressed()), this, SLOT(addScan()) );
  connect(ui.pbFilterMap, SIGNAL(pressed()), this, SLOT(filterMap()) );
  connect(ui.pbFilterMap_init, SIGNAL(pressed()), this, SLOT(filterMap_init()) );
  connect(ui.pbFilterMap_step, SIGNAL(pressed()), this, SLOT(filterMap_step()) );
  connect(ui.pbFilterMap_finish, SIGNAL(pressed()), this, SLOT(filterMap_finish()) );

  connect(ui.pbReg1, SIGNAL(pressed()), this, SLOT(reg1()) );
  connect(ui.pbReg2, SIGNAL(pressed()), this, SLOT(reg2()) );
  connect(ui.pbReg3, SIGNAL(pressed()), this, SLOT(reg3()) );

  connect(ui.pbPtPlus, SIGNAL(pressed()), this, SLOT(incPtSize()) );
  connect(ui.pbPtMinus, SIGNAL(pressed()), this, SLOT(decPtSize()) );
  connect(ui.sbMinHitCnt, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.cbVisuEstTraject, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuGPSTraject, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuMap, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.rbMapHeight, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapHeightSlices, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapDist, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapNormal, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapIntensity, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapNConf, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.cbOrigMap, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuCurrPoints, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuCurrNormals, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuLastPoints, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbVisuSubPoints, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbAdaption, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFiltering, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.sbFilterCol, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.sbFilterRow, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );

  imgSelect(ui.hsImageNb->value());
  glListIndex = glGenLists(1); // generate a display list
  update3D();
}


Visualizer3DMapper::~Visualizer3DMapper()
{
  continReader->stop();
  delete continReader;
  glDeleteLists(glListIndex, 1);
}

void Visualizer3DMapper::paintGLOpaque()
{
  glCallList(glListIndex+0);
}

void Visualizer3DMapper::paintGLTranslucent()
{
}

void Visualizer3DMapper::update3D()
{
  Map *map = mapper.getMap();
  LFrameSPtr currentFrame = mapper.getCurrentFrame();
  LFrameSPtr lastFrame = mapper.getLastFrame();
  glDeleteLists(glListIndex, 1);
  glListIndex = glGenLists(1); // generate a display list
  glNewList(glListIndex, GL_COMPILE);
  GLdouble clr[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, &(clr[0]));
  bool invert = (clr[0]+clr[1]+clr[2])/3 > 0.5; // bright background -> 1, dark background -> 0

  glPushMatrix();
  glTranslatef(0.0, 0.0, VISU_3D_SHIFTUP); // Point-shift-up

  // render trajectory
  if (ui.cbVisuEstTraject->isChecked()) {
    GLUquadricObj *quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL);     // rendering style of the quadric (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
    for (Map::trajectory_const_iterator ti=map->beginTrajectory(); ti!=map->endTrajectory(); ++ti) {
      matrixTools::DVector vehiclePos = ti->getWorldPos();
      glPushMatrix();
      glTranslatef(vehiclePos[0],vehiclePos[1],vehiclePos[2]); // move coord sys to point location
      glColor4f(1.0, 0.2, 0.2, 0.6);
      gluSphere(quadric, 0.1, 20, 20); //quadObject, radius, slices (around z-axis), stacks (along z-axis)
      glPopMatrix();
    }
    gluDeleteQuadric(quadric);
  }
  if (ui.cbVisuGPSTraject->isChecked()) {
    GLUquadricObj *quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL);     // rendering style of the quadric (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
    for (Map::trajectory_const_iterator ti=map->beginInsTrajectory(); ti!=map->endInsTrajectory(); ++ti) {
      matrixTools::DVector vehiclePos = ti->getWorldPos();
      glPushMatrix();
      glTranslatef(vehiclePos[0],vehiclePos[1],vehiclePos[2]); // move coord sys to point location
      glColor4f(0.2, 1.0, 0.2, 0.6);
      gluSphere(quadric, 0.1, 20, 20); //quadObject, radius, slices (around z-axis), stacks (along z-axis)
      glPopMatrix();
    }
    gluDeleteQuadric(quadric);
  }

  // debug rendering
  if (ui.cbAdaption->isChecked()) {
    int col = ui.sbFilterCol->value();
    int row = ui.sbFilterRow->value();
    if ((col < currentFrame->point3D.getHorizSize()) && (row < currentFrame->point3D.getVertSize())) {
      const unsigned int nbNSearched(6);
      const unsigned int nbNMinFound(3);
      const double nMaxDist(5*0.05);
      const double minNormalCoincidence(0.985); // cos(10°)=0,985
      DMatrix R2w(3,3); DVector t2w(3);
      HomogeneousTransformationMatrix::HTM_2_Rt(currentFrame->positionHTM2w, R2w, t2w);
      DVector p = ublas::prod(R2w,currentFrame->point3D.get(col,row))+t2w;
      DVector n = ublas::prod(R2w,currentFrame->normal3D.get(col,row));
      if ((norm_2(p) < 10000) && (norm_2(n) < 10000)) { // DBL_MAX
        glPointSize(6);
        glBegin(GL_POINTS);
        glColor3f(1.0,0.1,0.1); // point in red
        glVertex3f(p[0],p[1],p[2]);
        glEnd();
        glLineWidth(4);
        glBegin(GL_LINES);
        glColor3f(1.0,0.5,0.5); // point normal in redish
        glVertex3f(p[0],p[1],p[2]);
        glVertex3f(p[0]+n[0],p[1]+n[1],p[2]+n[2]);
        glEnd();
        Map::SearchResults neighborIt = map->findClosestNeighbors(p, n, nbNSearched);
        unsigned int nbNValid = 0;
        double aSum = 0.0; // LS-solution for adaptation
        double bSum = 0.0; // LS-solution for adaptation
        double avgNormalCoincidence = 0.0;
        glBegin(GL_LINES);
        for (unsigned int i=0; i<nbNSearched; ++i) { // looping over neighbors, sorted by distance
          const DVector &q = neighborIt.getCorrespPoint();
          DVector di = p-q;
          double dist = norm_2(di);
          if (dist > nMaxDist)
            break;
          const DVector &m = neighborIt.getCorrespNormal();
          double a = ublas::inner_prod(n,m);
          double b = ublas::inner_prod(di,m);
          double w = a*exp(-dist);
          if ((avgNormalCoincidence+a)/(double)(nbNValid+1) < minNormalCoincidence) // drops below threshold using this correspondence
            break;
          glColor3f(1.0,1.0,1.0); // connection in white
          glVertex3f(p[0],p[1],p[2]);
          glVertex3f(q[0],q[1],q[2]);
          glColor3f(0.7,0.7,1.0); // connections normal in blueish
          glVertex3f(q[0],q[1],q[2]);
          glVertex3f(q[0]+m[0],q[1]+m[1],q[2]+m[2]);
          avgNormalCoincidence += a;
          aSum += a*w*a;
          bSum += a*w*b;
          nbNValid++;
          if (!neighborIt.next())
            break;
        }
        glEnd();
        if (nbNValid >= nbNMinFound) {
          avgNormalCoincidence /= (double)nbNValid;
          double c = - bSum / aSum; // LS-adaption-factor
          if (c < 0.5) {
//            cout << " c=" << c << flush;
            p += n * c;
            glPointSize(6);
            glBegin(GL_POINTS);
            glColor3f(0.1,1.0,0.1); // point in green
            glVertex3f(p[0],p[1],p[2]);
            glEnd();
          }
        }
      }
    }
  }
  if ((ui.cbFiltering->isChecked()) && (map->csi != map->end())) {
    BOOST_AUTO(mcsi,map->csi);
    double fac = 1.0;
    GLUquadricObj *quadric = gluNewQuadric();
//    cout << endl;
//    for (unsigned int i=0; i<3; ++i) {
      if (mcsi != map->end()) {
        const double nMaxDist(5*0.05);//mapResol);
        SurfaceProxy &s = *mcsi;
        const DVector &p = s->origPosition;
//        cout << " rendering " << p << flush;
        const DVector &n = s->normal;
        Map::SearchResults neighborIt = map->findClosestNeighbors(s->position, s->normal, 15);
        glPointSize(6);
        glBegin(GL_POINTS);
        glColor3f(fac*1.0,fac*0.1,fac*0.1); // point in red
        glVertex3f(p[0],p[1],p[2]);
        glEnd();
        glLineWidth(4);
        glBegin(GL_LINES);
        glColor3f(fac*1.0,fac*0.5,fac*0.5); // point normal in redish
        glVertex3f(p[0],p[1],p[2]);
        glVertex3f(p[0]+n[0],p[1]+n[1],p[2]+n[2]);
        glEnd();
        neighborIt.next();
        DVector avgM = DZeroVector(3);
        double avgNc = 0.0;
        double sumW = 0.0;
        double sumWd = 0.0;
        double aSum = 0.0;
        double bSum = 0.0;
        unsigned int nbNValid = 0;
        for (unsigned int i=0; i<15; ++i) { // looping over neighbors, sorted by distance
          const DVector &q = neighborIt.getCorrespOrigPoint(); // getCorrespPoint() ? in that case cannot overwrite "position" below
          const DVector &m = neighborIt.getCorrespNormal();
          double        nci = neighborIt.getNormConf();
          DVector di = p-q;
          double dist = norm_2(di);
          if (dist > nMaxDist)
            break;
          if (dist > 0.0) { // skip search point as (probable) first correspondence
            double ai = ublas::inner_prod(n,m); // =1 if normals are in same direction
            double bi = ublas::inner_prod(di,m); // =dist if other point lies in normal direction
            double wid = exp(-dist);
            double wi = (0.5+0.5*fabs(bi/dist))*max(0.0,ai)*exp(-dist);
            sumW += wi;
            sumWd += wid;
            avgNc += wid*nci;
            avgM += wid*m;
      //      avgNormalCoincidence += wi*ai;
            aSum += ai*wi*ai;
            bSum += ai*wi*bi;
            nbNValid++;
          }
          glBegin(GL_LINES);
          glColor3f(fac*1.0,fac*1.0,fac*1.0); // connection in white
          glVertex3f(p[0],p[1],p[2]);
          glVertex3f(q[0],q[1],q[2]);
          glColor3f(fac*0.7,fac*0.7,fac*1.0); // connections normal in blueish
          glVertex3f(q[0],q[1],q[2]);
          glVertex3f(q[0]+m[0],q[1]+m[1],q[2]+m[2]);
          glEnd();
//          if (ui.cbVisuSubPoints->isChecked()) {
//            glColor4f(0.9, 0.9, 0.9, 0.5);
//            visAlg.render3DNormalPlane(q, m, 0.5, quadric);
//          }
          if (!neighborIt.next())
            break;
        }
        const double minNormalCoincidence(0.985); // cos(10°)=0,985
        const double minNormalConfidence(0.9);
        const unsigned int nbNMinFound(3);
        const double minWSum(1.5); // dist=0.5->exp(-dist)=0.6
        if ((sumW >= 0.1) && (sumWd >= minWSum) && (nbNValid >= nbNMinFound)) {
          avgNc /= sumWd;
          double avgMLength = norm_2(avgM);
          if (avgMLength > 0.0) // avoid div by zero
            avgM /= norm_2(avgM); // norm to length 1
          if ((avgNc <= minNormalConfidence) || (ublas::inner_prod(avgM,n) >= minNormalCoincidence)) {
            double c = - bSum / aSum; // LS-adaption-factor
            if (c < 0.5) {
//              cout << " c=" << c << flush;
              glPointSize(6);
              glBegin(GL_POINTS);
              glColor3f(0.1,1.0,0.1); // point in green
              glVertex3f(p[0]+n[0]*c,p[1]+n[1]*c,p[2]+n[2]*c);
              glEnd();
            }
          }
        }
        fac *= 0.7;
        ++mcsi;
      }
//    }
    gluDeleteQuadric(quadric);
  }

  // render map
  if (ui.cbOrigMap->isChecked()) {
    unsigned int minHitCnt = ui.sbMinHitCnt->value();
    bool filter = ((ui.cbFiltering->isChecked()) && (map->csi != map->end()));
    DVector drawCenter(3);
    if (filter)
      drawCenter = (**map->csi).origPosition;

    glPointSize(ptSize);
    glBegin(GL_POINTS);
    Map::const_iterator mapIt = map->begin();
    Map::const_iterator mapEnd = map->end();
    while (mapIt != mapEnd) {
//      const DVector &p = (**mapIt).position;
      const DVector &po = (**mapIt).origPosition;
//      if (norm_2(p-po) > 0.01) { // if point was adapted draw also original position
      if (((!filter) || (norm_2(drawCenter-po) < 1.0)) && (mapIt->hitCount >= minHitCnt)) {
        glColor3f(1.0, 0.5, 0.0); // orange
        glVertex3f(po[0], po[1], po[2]);
      }
      ++mapIt;
    }
    glEnd();
  }
  if (ui.cbVisuMap->isChecked()) {
    unsigned int minHitCnt = ui.sbMinHitCnt->value();
    glPointSize(ptSize);
    glBegin(GL_POINTS);
    Map::const_iterator mapIt = map->begin();
    Map::const_iterator mapEnd = map->end();
    while (mapIt != mapEnd) {
      if (mapIt->hitCount >= minHitCnt) {
        const DVector &p = (**mapIt).position;
        if (ui.rbMapNormal->isChecked()) { // color by normal vector
          const DVector &n = (**mapIt).normal;
          float r = fabs(n[0]);
          float g = fabs(n[1]);
          float b = fabs(n[2]);
//          int h,s,v; // max: 360/255/255
//          RGB2HSV((int)(r*255.0), (int)(g*255.0), (int)(b*255.0), h, s, v);
//          cout << endl << "h/s/v: " << h << "/" << s << "/" << v << flush;
//          int ri,gi,bi;
//          HSV2RGB(h, 255, 255, ri, gi, bi); // normalize
//          r = (float)ri/255.0;
//          g = (float)gi/255.0;
//          b = (float)bi/255.0;
          glColor3f(r, g, b);
        }
        if (ui.rbMapNConf->isChecked()) { // color by normal vector confidence
          const double nc = (**mapIt).normalConfidence;
          float r = 1.0-nc;
          float g = nc;
          float b = 0.5 - fabs(nc-0.5);
          glColor3f(r, g, b);
        }
        if (ui.rbMapHeight->isChecked()) { //  color by height
          float height = p(2);
          float r, g, b;
          if (height < 0) {
            float ratio = min(1.0,height/(-2.5));
            r = invert  ?  0.2        :  1.0-ratio;
            g = invert  ?  (1-ratio)  :  1.0-ratio;
            b = invert  ?  ratio      :  1.0;
          } else {
            float ratio = min(1.0,height/0.5);
            r = invert  ?  ratio      :  1.0;
            g = invert  ?  (1-ratio)  :  1.0-ratio; // normalize to range 0..1
            b = invert  ?  0.2        :  1.0-ratio;
          }
          glColor3f(r, g, b);
        }
        if (ui.rbMapHeightSlices->isChecked()) { //  color by height
          const float slice = 3.0;  // height in meter
          float relHeight = fmod(p(2)/slice, 1.0); // map to slices -1..1
          relHeight = fmod(relHeight+1.0, 1.0); // map to 0..1
          glColorHSV((int)(relHeight*360.0f), 255, 255);
        }
        if (ui.rbMapDist->isChecked()) { // color by distance from scanner
          float v = 0.2 + 0.8*max(0.0, 1.0 - (**mapIt).distFromScanner / VISU_DIST_MAX);
          glColor3f(v, v, v);
        }
        if (ui.rbMapIntensity->isChecked()) { // color by intensity
          unsigned char intens = (**mapIt).intensity;
          glColor3ub(intens, intens, intens);
        }
        glVertex3f(p[0], p[1], p[2]);
      }
      ++mapIt;
    }
    glEnd();
  }

  // render current scan
  if ((ui.cbVisuCurrPoints->isChecked()) && (currentFrame.use_count()>0)) {
    glPushMatrix();
    DCMatrix htm = currentFrame->positionHTM2w; // convert to column-major matrix
//    cout << " hmt_scan:" << htm << endl;
    glMultMatrixd(htm.data().begin());
    LidarImageVisualization::ColorMode3D mode = invert ? LidarImageVisualization::ColorFixed : LidarImageVisualization::ColorHeight;
    LidarImageVisualization::render3DPoints(currentFrame->point3D, mode, ptSize+1, 0.1, 0.1, 0.1);
    //LidarImageVisualization::render3DPoints(currentFrame->tmpVec3D, ptSize+1, 0.2, 1.0, 0.0); // original
    //LidarImageVisualization::render3DPoints(currentFrame->point3D, ptSize+2, 1.0, 0.5, 0.1); // adapted
    glPopMatrix();
  }
  if ((ui.cbVisuCurrNormals->isChecked()) && (currentFrame.use_count()>0)) {
    glPushMatrix();
    DCMatrix htm = currentFrame->positionHTM2w; // convert to column-major matrix
    glMultMatrixd(htm.data().begin());
    //LidarImageVisualization::render3DNormals(currentFrame->point3D, currentFrame->normal3D, &currentFrame->normalConfidence, &currentFrame->normalStdDevRAD, ptSize); // this renders cones
    LidarImageVisualization::render3DNormals(currentFrame->point3D, currentFrame->normal3D, &currentFrame->normalConfidence, NULL, NULL, ptSize, 0.1);
    glPopMatrix();
  }
  if ((ui.cbVisuLastPoints->isChecked()) && (lastFrame.use_count()>0)) {
    glPushMatrix();
    DCMatrix htm = lastFrame->positionHTM2w; // convert to column-major matrix
//    cout << " htm_last_scan:" << htm << endl;
    glMultMatrixd(htm.data().begin());
    LidarImageVisualization::render3DPoints(lastFrame->point3D, LidarImageVisualization::ColorFixed, ptSize+1, 0.5, 0.5, 0.5);
    glPopMatrix();
  }
  if ((ui.cbVisuSubPoints->isChecked()) && (lastFrame.use_count()>0)) {
    glPushMatrix();
    DCMatrix htm = map->getSubPtsHTM(); // convert to column-major matrix
//    cout << " htm_subsampled_points:" << htm << endl;
    glMultMatrixd(htm.data().begin());
    Map::const_iterator_samplepoints pB = map->beginSubsampledPoints();
    Map::const_iterator_samplepoints pE = map->endSubsampledPoints();
    glColor3f(1.0, 0.2, 0.2);
    glPointSize(ptSize+2);
    glBegin(GL_POINTS);
    while (pB != pE) {
      glVertex3f((*pB)[0], (*pB)[1], (*pB)[2]);
      ++pB;
    }
    glEnd();
    glPopMatrix();
  }

  glPopMatrix(); // shiftup
  glEndList();
  emit stateChanged();
}

void Visualizer3DMapper::imgSelect(int val)
{
  path fullfile(mapper.getImgName(val));
  ui.lImageName->setText(QString(fullfile.leaf().c_str()));
}

void Visualizer3DMapper::processSelectedImage()
{
  continReader->stop();
  readScan();
  if (!ui.cbGPS->isChecked())
    registerScan();
  addScan();
  setTimer(ui.pbProcessContinuous->isChecked());
}

void Visualizer3DMapper::processCertainNumber()
{
  unsigned int nb = ui.sbReadNumber->value();
  while (nb > 0) {
    processSelectedImage();
    --nb;
  }
}

void Visualizer3DMapper::setTimer(bool pressed)
{
//  ui.pbReadNumber->setEnabled(!pressed);
//  ui.pbReadNumber->setEnabled(!pressed);
  if ((pressed) && (ui.hsImageNb->tickPosition() < ui.hsImageNb->maximum()))
    continReader->start(ui.dsbSleep->value()*1000); // interval in ms
  else
    continReader->stop();
}

void Visualizer3DMapper::saveMap()
{
  QString file = QFileDialog::getSaveFileName(this, "Save map as ...", lastMapPath, "Map (*.map);;TextMap (*.txt);;XmlMap (*.xml)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  cout << endl << "saving map to " << file.toStdString() << "..." << flush;
  mapper.saveMap(file.toStdString());
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::loadMap()
{
  QString file = QFileDialog::getOpenFileName(this, "Open an existing map", lastMapPath, "All Maps (*.map *.txt *.xml);;Map (*.map);;TextMap (*.txt);;XmlMap (*.xml)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  cout << endl << "loading map from " << file.toStdString() << "..." << flush;
  mapper.loadMap(file.toStdString());
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::exportMap()
{
  QString file = QFileDialog::getSaveFileName(this, "Export map as ...", lastMapPath, "Text file (*.txt)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  cout << endl << "exporting map to " << file.toStdString() << "..." << flush;
  mapper.exportMap(file.toStdString());
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::exportMap3D()
{
  unsigned int resolFac = 1;
  QString file = QFileDialog::getSaveFileName(this, "Export map as ...", lastMapPath, "All Point Clouds (*.p3d *.p3di *.txt);;Point Cloud (*.p3d);;Point Cloud with Intensity (*.p3di)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  resolFac = QInputDialog::getInt( this, tr("Specify map resolution"), tr("resolution factor:"), resolFac, 1 );
  cout << endl << "exporting map to " << file.toStdString() << " with resolution factor " << resolFac << "..." << flush;
  mapper.exportMapPt3D(file.toStdString(), resolFac, (file[file.size()-1] == QChar('i')) ? Map::PointIntens : Map::Points);
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::importMap()
{
  QString file = QFileDialog::getOpenFileName(this, "Import an existing map", lastMapPath, "Text file (*.txt)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  cout << endl << "importing map from " << file.toStdString() << "..." << flush;
  mapper.getMap()->importMap(file.toStdString());
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::importMap3D()
{
  QString file = QFileDialog::getOpenFileName(this, "Import an existing map", lastMapPath, "All Point Clouds (*.p3d *.p3di *.txt);;Point Cloud (*.p3d);;Point Cloud with Intensity (*.p3di)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  cout << endl << "importing map from " << file.toStdString() << "..." << flush;
  bool withIntensity = (file[file.size()-1] == QChar('i'));
  mapper.getMap()->importPointCloud(file.toStdString(), withIntensity);
  cout << "done" << flush;
  update3D();
}

void Visualizer3DMapper::exportTraj()
{
  QString file = QFileDialog::getSaveFileName(this, "Export trajectory as ...", lastMapPath, "Trajectory file (*.traj)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  QString basename = file.section('.', 0, -2); // strip file extension
  cout << endl <<"stripped file " << file.toStdString() << " --> " << basename.toStdString() << endl;
  mapper.getMap()->exportTrajectories(basename.toStdString());
  update3D();
}

void Visualizer3DMapper::importTraj()
{
  QString file = QFileDialog::getOpenFileName(this, "Import an existing trajectory", lastMapPath, "Trajectory file (*.traj)");
  if (file.isEmpty()) return;
  lastMapPath = file.section('/', 0, -2);
  QString basename = file.section('.', 0, -2); // strip file extension
  cout << endl <<"stripped file " << file.toStdString() << " --> " << basename.toStdString() << endl;
  mapper.getMap()->importTrajectories(basename.toStdString());
  update3D();
}

void Visualizer3DMapper::resetMapping()
{
  mapper.reset();
  update3D();
}

void Visualizer3DMapper::readScan()
{
  unsigned int fileIdx = ui.hsImageNb->value();
  mapper.readScan(fileIdx);

  // increment slider-position, if already last position, stop timer
  if (ui.hsImageNb->tickPosition() < ui.hsImageNb->maximum()) {
    ui.hsImageNb->setValue(fileIdx+1);
  }

  QImage img;
  LFrameSPtr currentFrame = mapper.getCurrentFrame();
  visAlg.renderColored(currentFrame->distance, img, 1.5, VISU_DIST_MAX);
  emit redraw2D(img);
  update3D();
}

void Visualizer3DMapper::registerScan()
{
  mapper.registerScan(ui.hsSubsampling->value(),ui.cbUnwarp->isChecked());
//  if (ui.cbUnwarp->isChecked())
//    map->registerScan(currentFrame, lastFrame, ui.hsSubsampling->value(),ui.cbUnwarp->isChecked());
  update3D();
}
void Visualizer3DMapper::reg1()
{
  mapper.reg1(ui.hsSubsampling->value(), ui.cbUnwarp->isChecked());
  update3D();
}
void Visualizer3DMapper::reg2()
{
  mapper.reg2();
  update3D();
}
void Visualizer3DMapper::reg3()
{
  mapper.reg3();
  update3D();
}

void Visualizer3DMapper::filterMap()
{
//  map->filter();
  mapper.filterMap_init();
  mapper.filterMap_finish();
  update3D();
}
void Visualizer3DMapper::filterMap_init()
{
  mapper.filterMap_init();
  update3D();
}
void Visualizer3DMapper::filterMap_step()
{
  mapper.filterMap_step();
  update3D();
}
void Visualizer3DMapper::filterMap_finish()
{
  mapper.filterMap_finish();
  update3D();
}

void Visualizer3DMapper::addScan()
{
  mapper.addScan(ui.cbAdept->isChecked(), ui.cbUnwarp->isChecked(), ui.cbRemOld->isChecked(), ui.cbRemFar->isChecked(), ui.cbGPS->isChecked());
  update3D();
}

