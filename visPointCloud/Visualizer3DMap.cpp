#include "Visualizer3DMap.hpp"

#include <sys/time.h>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <QFileDialog>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "ConsoleProgressBar.hpp"
#include <Gui3DQt/graphics.hpp>


using namespace std;
using namespace boost::filesystem;
using namespace Gui3DQt;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////        Visualizer3DMap        ///////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

const double VISU_DIST_MAX = 80.0; // in meter

Visualizer3DMap::Visualizer3DMap(std::string mapfile, int dummy, int rangeFrom_, int rangeTo_, bool loadOnlyTraj, QWidget *parent)
  : Gui3DQt::Visualizer(parent)
    ,filename(mapfile)
    ,rangeFrom(rangeFrom_)
    ,rangeTo(rangeTo_)
    ,ptSize(1)
{
  glListIndex = glGenLists(5); // generate a display list
  for (unsigned int i=0;i<5;++i) {
    glNewList(glListIndex+i, GL_COMPILE);
    glEndList();
  }

  ui.setupUi(this);

  Qt::WindowStates flags = windowState();
  flags = Qt::WindowMaximized;
  setWindowState(flags);

  connect(ui.pbLoadMap, SIGNAL(pressed()), this, SLOT(loadMap()) );
  connect(ui.pbReloadMap, SIGNAL(pressed()), this, SLOT(reloadMap()) );
  connect(ui.pbPtPlus, SIGNAL(pressed()), this, SLOT(incPtSize()) );
  connect(ui.pbPtMinus, SIGNAL(pressed()), this, SLOT(decPtSize()) );
  connect(ui.dsbShiftup, SIGNAL(valueChanged(double)), this, SLOT(emitStateChanged()) );
  connect(ui.cbVisuCS, SIGNAL(toggled(bool)), this, SLOT(emitStateChanged()) );
  connect(ui.cbVisuEstTraject, SIGNAL(toggled(bool)), this, SLOT(emitStateChanged()) );
  connect(ui.cbVisuGPSTraject, SIGNAL(toggled(bool)), this, SLOT(emitStateChanged()) );
  connect(ui.cbVisuMap, SIGNAL(toggled(bool)), this, SLOT(updateMap3D()) );
  connect(ui.hsColorShift, SIGNAL(valueChanged(int)), this, SLOT(updateMap3D()));
  connect(ui.cbColorInvert, SIGNAL(toggled(bool)), this, SLOT(updateMap3D()));
  connect(ui.dsbSliceHeight, SIGNAL(valueChanged(double)), this, SLOT(sliceHeightChanged()) );
  connect(ui.dsbHeightPitch, SIGNAL(valueChanged(double)), this, SLOT(heightAngleChanged()) );
  connect(ui.dsbHeightRoll, SIGNAL(valueChanged(double)), this, SLOT(heightAngleChanged()) );
  connect(ui.rbMapHeight, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapSlices, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapNormal, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapNConf, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbMapSDist, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );

  connect(ui.gbLimit, SIGNAL(toggled(bool)), this, SLOT(updateLimit3D()) );
  connect(ui.cbLimitApply, SIGNAL(toggled(bool)), this, SLOT(updateMap3D()) );
  connect(ui.pbLimitShow, SIGNAL(toggled(bool)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMinX, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMinY, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMinZ, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMaxX, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMaxY, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );
  connect(ui.dsbMaxZ, SIGNAL(valueChanged(double)), this, SLOT(updateLimit3D()) );

  checkFileType();
  reloadTraj();
  if (!loadOnlyTraj)
    reloadMap();

}

Visualizer3DMap::~Visualizer3DMap()
{
  glDeleteLists(glListIndex, 5);
}

void Visualizer3DMap::paintGLOpaque()
{
  glPushMatrix();
  glTranslatef(0.0, 0.0, -tminZ+ui.dsbShiftup->value()); // Point-shift-up

  if (ui.cbVisuMap->isChecked()) {
//    glCallList(glListIndex+0); // map
    renderer.render(ptSize, &indices);
  }
  if (ui.cbVisuCS->isChecked())
    glCallList(glListIndex+2); // cs
  if (ui.cbVisuEstTraject->isChecked())
    glCallList(glListIndex+3); // traj 1
  if (ui.cbVisuGPSTraject->isChecked())
    glCallList(glListIndex+4); // traj 2

  glPopMatrix();
}

void Visualizer3DMap::paintGLTranslucent()
{
  glPushMatrix();
  glTranslatef(0.0, 0.0, -tminZ+ui.dsbShiftup->value()); // Point-shift-up
  glCallList(glListIndex+1); // box
  glPopMatrix();
}

void Visualizer3DMap::checkFileType()
{
  mapIsP3d = !(filename.substr(filename.size()-3, 3) == "txt");
  mapHasIntensity = mapIsP3d && (filename[filename.size()-1] == 'i' || filename[filename.size()-2] == 'i'); // .p3di .p3dib
  mapIsBinary = mapIsP3d && (filename[filename.size()-1] == 'b'); // .p3db .p3dib
}

void Visualizer3DMap::renderTraj(const deque<GlVec3> &trajectory, bool colorByHeight, bool renderSpheres)
{
  float trMinZ = DBL_MAX;
  float trMaxZ = -DBL_MAX;
  if (colorByHeight) {
    BOOST_FOREACH(const GlVec3 &pt, trajectory) {
      trMinZ = min(trMinZ,pt.z);
      trMaxZ = max(trMaxZ,pt.z);
    }
  }
  float zRange = trMaxZ-trMinZ;

  if (renderSpheres) {
    GLUquadricObj *quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL);     // rendering style of the quadric (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
    BOOST_FOREACH(const GlVec3 &pt, trajectory) {
      glPushMatrix();
      glTranslatef(pt.x,pt.y,pt.z); // move coord sys to point location
      gluSphere(quadric, 0.1, 20, 20); //quadObject, radius, slices (around z-axis), stacks (along z-axis)
      glPopMatrix();
    }
    gluDeleteQuadric(quadric);
  }

  glLineWidth(ptSize+1);
  glBegin(GL_LINES);
  unsigned int i=0;
  unsigned int s=trajectory.size();
  BOOST_FOREACH(const GlVec3 &pt, trajectory) {
    if (colorByHeight)
      glColorHSV((int)(120.0f-120.0f*(pt.z-trMinZ)/zRange), 255, 255); // H = 0..259; S,V=0..255
    if ((i!=0) && (i != s-1)) // don't double-draw at start and end
      glVertex3f(pt.x,pt.y,pt.z);
    glVertex3f(pt.x,pt.y,pt.z);
    ++i;
  }
  glEnd();
}

void Visualizer3DMap::updateMap3D()
{
  // update the colors and indices vector
  indices.clear();
  indices.reserve(renderer.size());
  if (ui.cbVisuMap->isChecked()) {
    recalcWorldMinMaxZ();
    unsigned int mode = 0; // (ui.rbMapHeight->isChecked())
    if (ui.rbMapSlices->isChecked()) mode = 1;
    if (ui.rbMapNormal->isChecked()) mode = 2;
    if (ui.rbMapNConf->isChecked()) mode = 3;
    if (ui.rbMapSDist->isChecked()) mode = 4;
    const bool  limitBox = ui.cbLimitApply->isChecked();
    const float limitMinX = ui.dsbMinX->value();
    const float limitMaxX = ui.dsbMaxX->value();
    const float limitMinY = ui.dsbMinY->value();
    const float limitMaxY = ui.dsbMaxY->value();
    const float limitMinZ = ui.dsbMinZ->value();
    const float limitMaxZ = ui.dsbMaxZ->value();
    const float sliceHeight = min((float)ui.dsbSliceHeight->value(),pmaxZ-pminZ);
    const float roll = ui.dsbHeightRoll->value()*M_PI/180.0f;
    const float pitch = ui.dsbHeightPitch->value()*M_PI/180.0f;
    const int   colorIdxShift = ui.hsColorShift->value()*360/100;
    const bool  colorIdxInv = ui.cbColorInvert->isChecked();
    cout << endl << "rendering map:" << endl;
    // drawing indices are just all the points (0, 1, ..., n-1)

    int r,g,b;
    unsigned int div = max(renderer.size()/1000,(size_t)1);
    ConsoleProgressBar progress(max(renderer.size()/div, (size_t)1));

    for (unsigned int i=0; i<renderer.size(); ++i) {
      if (i%div == 0)
        progress.increment();
      const GlVec3 &pt = renderer.pointAt(i);
      const Attrib &at = attribs[i];
      const float rotZ = rotateZ(pt.x, pt.y, pt.z, roll, pitch);
      GlCol3 &cl = renderer.colorAt(i);
      if ((limitBox) && ( (pt.x < limitMinX) || (pt.x > limitMaxX)
                          || (pt.y < limitMinY) || (pt.y > limitMaxY)
                          || (rotZ < limitMinZ) || (rotZ > limitMaxZ)))
        continue;
      indices.push_back(i);
      switch (mode) {
        case 0 : { // color complete range as one color-slice
          int ci = (int)((max(min(rotZ,pmaxZ+0.3f),pminZ-2.0f)-(pminZ-2.0f))*360.0/((pmaxZ+0.3f)-(pminZ-2.0f)));
          //int ci = (int)((max(min(rotZ,pmaxZ),pminZ)-pminZ)*255.0/(pmaxZ-pminZ));
          ci = (ci + colorIdxShift) % 360;
          if (colorIdxInv) ci = 360 - ci - 1;
          HSV2RGB(ci*240/360, 255, 255, r,g,b); // don't use full spectrum, but only red to blue
          cl = GlCol3(r, g, b);
        } break;
        case 1 : { // color slices of 2m repeatedly
          int ci = ((int)((rotZ-pminZ)*360.0/sliceHeight)) % 360; // 0..259
          ci = (ci + colorIdxShift) % 360; // 0..259
          if (colorIdxInv) ci = 360 - ci - 1;
          HSV2RGB(ci, 255, 255, r,g,b);
          cl = GlCol3(r, g, b);
        } break;
        case 2 : {
          int ci = at.normalDirectionColorIndex;
          HSV2RGB(ci*360/255, 255, 255, r,g,b);
          cl = GlCol3(r, g, b);
        } break;
        case 3 : {
          int ci = at.normalConfidenceColorIndex;
          cl = GlCol3(255-ci, ci, 128-abs(ci-128));
        } break;
        case 4 : {
          int ci = at.scannerDistColorIndex;
          cl = GlCol3(ci, ci, ci);
        } break;
      }
    } // end: loop over points
    progress.finish();
  } // end: visu map activated


//  glNewList(glListIndex+0, GL_COMPILE);
//  glPointSize(ptSize);
//  glDrawElements(GL_POINTS, /* Primitivtyp */
//      indices.size(), /* Anzahl Indizes */
//      GL_UNSIGNED_INT, /* Typ der Indizes */
//      &indices[0]); /* Index-Array */
//
//  glEndList();
  emit stateChanged();
}

void Visualizer3DMap::updateLimit3D()
{
  glNewList(glListIndex+1, GL_COMPILE);
  if ((ui.pbLimitShow->isChecked()) && (ui.gbLimit->isChecked())) {
    glColor4f(0.5, 0.5, 0.5, 0.5);
    Graphics::draw_cube_cage( ui.dsbMinX->value(), ui.dsbMaxX->value(),
                              ui.dsbMinY->value(), ui.dsbMaxY->value(),
                              ui.dsbMinZ->value(), ui.dsbMaxZ->value());
  }
  glEndList();
  if (ui.cbLimitApply->isChecked())
    updateMap3D();
  emit stateChanged();
}

void Visualizer3DMap::updateTraj3D()
{
  glNewList(glListIndex+3, GL_COMPILE);
  glColor3f(1.0, 0.2, 0.2);
  renderTraj(traj, true);
  glEndList();

  glNewList(glListIndex+4, GL_COMPILE);
  glColor3f(0.2, 1.0, 0.2);
  renderTraj(trajINS);
  glEndList();

  glNewList(glListIndex+2, GL_COMPILE);
//  graphics_draw_coordinate_frame(1.0); // produces X errors
  glEndList();

  emit stateChanged();
}

void Visualizer3DMap::importTrajectory(deque<GlVec3> &traj, string filename)
{
  cout << endl << "- loading trajectory " << filename << "..." << flush;
  ifstream in(filename.c_str());
  if (!in.good())
    throw runtime_error("Map::importTrajectory: problem opening the file for reading");
  traj.clear();
  string line;
  int64_t timestamp;
  double htm2w[4][4];
  bool globErr = false;
  double oldPt[3];
  double currPt[3];
  double trajLength = 0.0;
  while (getline(in, line)) { // read one line into string
    istringstream istr;   // streaming the string containing a line of the calib_file
    istr.str(line);
    istr >> timestamp;
    bool err = false;
    for (unsigned int row=0; row<4;++row) {
      for (unsigned int col=0; col<4;++col) {
        if (istr.eof())
          err=true;
        istr >> htm2w[row][col];
      }
    }
    if (err)
      globErr = true;
    else {
      memcpy(oldPt, currPt, sizeof(double)*3);
      currPt[0]=htm2w[0][3];
      currPt[1]=htm2w[1][3];
      currPt[2]=htm2w[2][3];
      if (!traj.empty())
        trajLength += sqrt(pow(oldPt[0]-currPt[0],2)+pow(oldPt[1]-currPt[1],2)+pow(oldPt[2]-currPt[2],2));
      traj.push_back(GlVec3(currPt));
    }
  }
  if (globErr)
    cerr << "ERROR" << flush;
  else
    cout << "imported with length " << trajLength << flush;
}

void Visualizer3DMap::loadMap()
{
  QString file;
  file.fromStdString(filename);
  file = QFileDialog::getOpenFileName(this, "Import an existing map", file.section('/', 0, -2), "All valid files (*.p3d *.p3di *.txt);;Maps (*.txt);;Point Clouds (*.p3d);;Point Clouds with Intensity (*.p3di)");
  if (file.isEmpty()) return;
  filename = file.toStdString();
  checkFileType();
  reloadTraj();
  reloadMap();
}

template <typename T>
T getVal(ifstream &ifs, bool binary) {
  T t;
  if (binary)
    ifs.read(reinterpret_cast<char *>(&t),sizeof(T));
  else
    ifs >> t;
  return t;
}

template <>
unsigned char getVal(ifstream &ifs, bool binary) {
  unsigned char t;
  if (binary) {
    ifs.read(reinterpret_cast<char *>(&t),sizeof(unsigned char));
  } else {
    unsigned int tt;
    ifs >> tt;
    t = tt;
  }
  return t;
}


void Visualizer3DMap::reloadTraj()
{
  if (mapIsP3d)
    return;

  // might throw (use try-catch?):
  importTrajectory(traj, filename + ".traj");
  importTrajectory(trajINS, filename + ".ins.traj");

  if (traj.size() != trajINS.size())
    cerr << endl << "WARNING! estimated trajectory has size " << traj.size() << " but INS has size " << trajINS.size() << flush;

  //set min/maxZ of trajectory
  tminZ = DBL_MAX;
  tmaxZ = -DBL_MAX;
  BOOST_FOREACH(const GlVec3 &pt, traj) {
    tminZ = min(tminZ,pt.z);
    tmaxZ = max(tmaxZ,pt.z);
  }

  updateTraj3D();
}

void Visualizer3DMap::reloadMap()
{
  // estimate number of points in file
  const unsigned int bytesPerPointTxt = 100; // 100 characters/bytes per point
  const unsigned int bytesPerPointP3D = 25; // 25 characters/bytes per point +4 if intensity
  const unsigned int bytesPerPointP3DI = bytesPerPointP3D+4;
  const unsigned int bytesPerPointP3DB = 3*sizeof(double); // 3*sizeof(double) bytes per point +1 if intensity
  const unsigned int bytesPerPointP3DIB = bytesPerPointP3DB+1;
  unsigned int bytesPerPoint;
  if (mapIsP3d) {
    if (mapIsBinary) {
      bytesPerPoint = mapHasIntensity ? bytesPerPointP3DIB : bytesPerPointP3DB;
    } else {
      bytesPerPoint = mapHasIntensity ? bytesPerPointP3DI : bytesPerPointP3D;
    }
  } else {
    bytesPerPoint = bytesPerPointTxt;
  }
  unsigned long filesize = file_size(filename);
  unsigned long expectedFilePts = filesize / bytesPerPoint;
  long expectedLoadPts = expectedFilePts;
  if (rangeTo > 0) expectedLoadPts = min(expectedLoadPts,(long)rangeTo);
  if (rangeFrom > 0) expectedLoadPts -= min(expectedLoadPts,(long)rangeFrom);
  long expectedReadLimit = expectedFilePts;
  if (rangeTo > 0) expectedReadLimit = min(expectedReadLimit,(long)rangeTo);

  // output information
  cout << endl << "- loading map" << flush;
  string type = (mapIsBinary) ? "binary" : "text";
  if (mapIsP3d) cout << " as " << type << " p3d" << flush;
  if (mapHasIntensity) cout << " with intensities" << flush;
  //cout << endl << "size of " << filename << " is " << filesize << endl;

  // load map
  pminX = FLT_MAX;
  pminY = FLT_MAX;
  pmaxX = -FLT_MAX;
  pmaxY = -FLT_MAX;
  ifstream in;
  if (mapIsBinary)
    in.open(filename.c_str(), ios::binary);
  else
    in.open(filename.c_str());
  if (!in.good())
    throw runtime_error("Map::importMap: problem opening the file for reading");
  // clear all scans stored until now
  renderer.clear(); renderer.reserve(expectedLoadPts*1.1);
  attribs.clear(); attribs.reserve(expectedLoadPts*1.1);
  string line;
  Attrib newAttr;
  GlVec3 newP3D;
  GlCol3 newCol;
  unsigned int dummyUI;
  double dummyD;
  int hitCnt = 1;
  double nx = 0.0;
  double ny = 0.0;
  double nz = 0.0;
  double scanD = 0.0;
  double nConf = 0.0;
  int h,s,v;
  //unsigned int hitCnt = 1;
  unsigned int i=0;
  unsigned int skipped=0;
  const int minHitCnt = ui.sbMinHitCnt->value();
  cout << endl << "about to load ~" << expectedLoadPts << " points" << endl;
  ConsoleProgressBar progress(max(expectedReadLimit/10000, 1l));
  while ((mapIsP3d && (!in.eof())) || (getline(in, line))) { // if p3d test for EOF, otherwise read one line into string
    if ((++i) % 10000 == 0)
      progress.increment();
    if ((rangeFrom >= 0) && (i < rangeFrom)) continue;
    if ((rangeTo > 0) && (i > rangeTo)) break;
    // load point data into temporary variables
    if (mapIsP3d) {
      newP3D.x = getVal<double>(in, mapIsBinary);
      newP3D.y = getVal<double>(in, mapIsBinary);
      newP3D.z = getVal<double>(in, mapIsBinary);
      if (mapHasIntensity)
        scanD = (double)(getVal<unsigned char>(in, mapIsBinary))/255.; // will be used to color map
    } else {
      istringstream istr;   // streaming the string containing a line of the calib_file
      istr.str(line);
      istr >> dummyUI; // pose-idx
      istr >> newP3D.x;
      istr >> newP3D.y;
      istr >> newP3D.z;
      istr >> dummyD; // orig-pos
      istr >> dummyD; // orig-pos
      istr >> dummyD; // orig-pos
      istr >> nx; // [-1..1]
      istr >> ny; // [-1..1]
      istr >> nz; // [-1..1]
      istr >> nConf; // tmp.normalConfidence; // [0..1]
      istr >> scanD; // s.distFromScanner;
      if (!istr.eof()) {
        istr >> hitCnt;
        if (hitCnt < minHitCnt) {
          ++skipped;
          continue;
        }
      }
      scanD = 0.2 + 0.8*max(0.0, 1.0 - scanD / VISU_DIST_MAX);
    }
    // generate color information and store
    RGB2HSV(abs(nx*255), abs(ny*255), abs(nz*255), h, s, v);
    newAttr.normalDirectionColorIndex = (unsigned char)(h*255/360); // map from 0..360 to 0..255
    newAttr.normalConfidenceColorIndex = (unsigned char)(nConf*255); // map to 0..255
    newAttr.scannerDistColorIndex = (unsigned char)(scanD*255); // map to 0..255
    pminX = min(pminX,newP3D.x);
    pminY = min(pminY,newP3D.y);
    pmaxX = max(pmaxX,newP3D.x);
    pmaxY = max(pmaxY,newP3D.y);
    renderer.push_back(newP3D,newCol);
    attribs.push_back(newAttr);
  } // end: loop over points in file
  progress.finish();
  cout << endl << "loaded " << renderer.size() << " points, skipped " << skipped << " points" << flush;
  recalcWorldMinMaxZ();
  ui.dsbMinX->setMaximum(pmaxX);
  ui.dsbMinX->setMinimum(pminX);
  ui.dsbMaxX->setMaximum(pmaxX);
  ui.dsbMaxX->setMinimum(pminX);
  ui.dsbMinY->setMaximum(pmaxY);
  ui.dsbMinY->setMinimum(pminY);
  ui.dsbMaxY->setMaximum(pmaxY);
  ui.dsbMaxY->setMinimum(pminY);
  ui.dsbMinZ->setMaximum(pmaxZ);
  ui.dsbMinZ->setMinimum(pminZ);
  ui.dsbMaxZ->setMaximum(pmaxZ);
  ui.dsbMaxZ->setMinimum(pminZ);
  if ((ui.dsbMinX->value() == 0.0) && (ui.dsbMaxX->value() == 0.0)) {
    ui.dsbMinX->setValue(pminX);
    ui.dsbMaxX->setValue(pmaxX);
  }
  if ((ui.dsbMinY->value() == 0.0) && (ui.dsbMaxY->value() == 0.0)) {
    ui.dsbMinY->setValue(pminY);
    ui.dsbMaxY->setValue(pmaxY);
  }
  if ((ui.dsbMinZ->value() == 0.0) && (ui.dsbMaxZ->value() == 0.0)) {
    ui.dsbMinZ->setValue(pminZ);
    ui.dsbMaxZ->setValue(pmaxZ);
  }
  ui.dsbShiftup->setMinimum(-pmaxZ);
  ui.dsbShiftup->setMaximum(-pminZ);
  if (mapIsP3d)
    ui.rbMapSDist->setText("Intensity");
  else
    ui.rbMapSDist->setText("Scan-Dist");
  ui.rbMapSDist->setEnabled((!mapIsP3d) || mapHasIntensity);
  ui.rbMapNormal->setEnabled(!mapIsP3d);
  ui.rbMapNConf->setEnabled(!mapIsP3d);

  // call render
  if (ui.rbMapHeight->isChecked())
    updateMap3D();
  else
    ui.rbMapHeight->setChecked(true);
}

void Visualizer3DMap::recalcWorldMinMaxZ()
{
  const float roll = ui.dsbHeightRoll->value()*M_PI/180.0f;
  const float pitch = ui.dsbHeightPitch->value()*M_PI/180.0f;
  pminZ = FLT_MAX;
  pmaxZ = -FLT_MAX;
  float z;
  for (size_t i=0; i<renderer.size(); ++i) {
    const Gui3DQt::PointCloudRenderer::GlVec3& pt = renderer.pointAt(i);
    z = rotateZ(pt.x, pt.y, pt.z, roll, pitch);
    pminZ = min(pminZ,z);
    pmaxZ = max(pmaxZ,z);
  }
}

