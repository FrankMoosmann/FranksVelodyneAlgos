#include "VisualizerTracks.hpp"

#include <sys/time.h>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <QFileDialog>
#include <QInputDialog>

#include <Gui3DQt/graphics.hpp>
#include <LidarImageVisualization.hpp>
#include <HomogeneousTransformationMatrix.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////        VisualizerTracks        ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

const float VISU_DIST_MAX = 80.0; // in meter
const float VISU_DIST_MIN = 5.0; // in meter
const float PRECACHE_WAIT = 50; // in msec
const float DEFAULT_SLICE = 2.0;
const float PASSAT_HEIGHT = 1.9;
const unsigned long long WORLD_TRACK_ID = 1;

inline float sqr(float x) {return x*x;};

unsigned long long VisualizerTracks::PrecacheThread::precache(int idx)
{
  if (vis.loadedTracks[idx] == NULL) {
    // data not yet loaded --> load into temporary buffer
    FrameData *newdata = vis.loadFrame(idx);
    vis.trackMutex.lock();
    if (vis.loadedTracks[idx] == NULL) // check again, data could have been loaded by main loop in the meantime
      vis.loadedTracks[idx] = newdata;
    else
      delete newdata;
    vis.trackMutex.unlock();
  }
  return vis.loadedTracks[idx]->getMemUsage();
}

void VisualizerTracks::PrecacheThread::freeMem(int idx)
{
  vis.trackMutex.lock();
  if (vis.loadedTracks[idx] != NULL) {
    //cout << "freeing frame " << vis.images[idx] << " (" << vis.loadedTracks[idx]->getMemUsage()/1024 << " kB)" << endl;
    delete vis.loadedTracks[idx];
    vis.loadedTracks[idx] = NULL;
  }
  vis.trackMutex.unlock();
}

void VisualizerTracks::PrecacheThread::run()
{
  cout << endl << "precaching..." << flush;
  int currIdx = vis.loadedTracks.size(); // this value can never be selected -> 1st if in while-loop fires
  int idxOffset = 0;
  unsigned long long trackMemUsed = 0;
  // loop around current index and load data
  while (trackMemUsed < bufferMem) {
    if (currIdx != (int)vis.ui.hsImageNb->value()) { // selected image has changed
      currIdx = vis.ui.hsImageNb->value();
      idxOffset = 0;
      vis.assertFrameData(currIdx); // locks mutex
      trackMemUsed = vis.loadedTracks[currIdx]->getMemUsage();
    }

    if ((currIdx < idxOffset) && (currIdx+idxOffset >= (int)vis.loadedTracks.size()))
      break; // index+-offset out of bounds -> precashing finished

    if (currIdx+idxOffset < (int)vis.loadedTracks.size()) // load next frame
      trackMemUsed += precache(currIdx + idxOffset);
    if (currIdx >= idxOffset) // load previous frame
      trackMemUsed += precache(currIdx - idxOffset);

    idxOffset++;
    usleep(PRECACHE_WAIT);
  }
  // delete all other loaded data
  for (int idx = currIdx+idxOffset; idx < (int)vis.loadedTracks.size(); ++idx) {
    freeMem(idx);
  }
  for (int idx = currIdx-idxOffset; idx >= 0; --idx) {
    freeMem(idx);
  }
  cout << "precached indices " << currIdx << "+-" << idxOffset << " eating " << trackMemUsed/1024 << "/" << bufferMem/1024 << " kB" << endl;
}

VisualizerTracks::VisualizerTracks(string trkDir, unsigned int bufferMemMB, bool loadPoints_, bool loadWorld_, float shiftup, VisuMode visMode, bool exitAfterLastFrame_, Gui3DQt::VisualizerPassat *passatVis_, QWidget *parent)
  : Gui3DQt::Visualizer(parent)
  , passatVis(passatVis_)
  , listIsUpdating(false)
  , currWorldTrack(NULL)
  , worldModelRenderer(NULL)
  , loadPoints(loadPoints_)
  , exitAfterLastFrame(exitAfterLastFrame_)
  , ptSizeScan(2)
  , ptSizeMovingTrk(2)
  , ptSizeStaticTrk(1)
  , ptSizeWorld(1)
{
  glListIndex = glGenLists(2); // generate a display list

  playTimer = new QTimer(this);
  playTimer->setSingleShot(false);
  connect(playTimer, SIGNAL(timeout()), this, SLOT(imgNext()) );

  precacheThread = new PrecacheThread(*this, bufferMemMB);

  // determine directory where image files are located
  string imgDir;
  string sourcedirfile = path(trkDir / path("sourcedir")).string();
  ifstream imgSourceFile(sourcedirfile.c_str());
  getline(imgSourceFile,imgDir);
  imgSourceFile.close();
  if (imgDir.length() == 0) {
    string errstr = "Source-dir-file " + sourcedirfile + " is empty";
    throw logic_error(errstr);
  }
  try {
    proj = new PNGImageProjector(imgDir+"/img.cfg");
  } catch (std::invalid_argument) {
    proj = new PNGImageProjector(imgDir+"/../velodyne_config.txt");
  }

  // determine track files and corresponding image files
  //cout << path(trkDir / path("sourcedir")).string() << " --> " << imgDir << endl;
  const path trackbase(trkDir);
  const path imagebase(imgDir);
  const boost::regex textfilter( ".*\\.txt" ); // "\\" is transformed into "\" at compile-time
  boost::smatch what; // match result
  list<string> tmpTrackList;
  if (exists(trackbase)) {
    directory_iterator end;
    for (directory_iterator iter(trackbase); iter != end; ++iter) {
      path currFile = *iter;
      if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
      if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
      if (!boost::regex_match(currFile.filename().string(), what, textfilter)) continue; // Skip if no match
      tmpTrackList.push_back(currFile.filename().string()); // File matches, store it
      //cout << "pushing " << currFile.filename() << endl;
    }
  }
  tmpTrackList.sort(); // sort alphabetically
  BOOST_FOREACH(string tf, tmpTrackList) {
    //cout << "stem is " << path(tf).stem() << endl;
    path imgFile = imagebase / path(tf).stem(); // without ".txt"
    path trkFile = trackbase / tf;
    //cout << "testing " << imgFile << " AND " << trkFile << endl;
    if (is_regular_file(imgFile) && is_regular_file(trkFile)) {
      //cout << "--> success" << endl;
      images.push_back(imgFile.string());
      tracks.push_back(trkFile.string());
    }
  }
  loadedTracks.resize(tracks.size());
  if (loadedTracks.size() == 0) {
    string errstr = "No corresponding files found in imgDir " + imgDir + " and trkDir " + trkDir;
    throw logic_error(errstr);
  }
  cout << "found " << tracks.size() << " corresponding images/trackfiles" << endl;

  // load world model
  if (loadWorld_) {
    string p3dFile = path(trkDir / path("world.p3d")).string();
    loadWorldModelRenderer(p3dFile);
  }

  // setup GUI
  ui.setupUi(this);
  ui.hsImageNb->setMinimum(0);
  ui.hsImageNb->setMaximum(loadedTracks.size()-1);
  ui.hsImageNb->setValue(0);
  ui.hsImageNb->setTracking(false);
  ui.dsbShiftUp->setValue(shiftup);
  ui.gbDrawWorld->setEnabled(worldModelRenderer);

  switch(visMode) {
    case VM_TBB : {
      // do nothing, default settings from UI model are used
    } break;
    case VM_TPC : {
      ui.rbImgColorConst->setChecked(true);
      ui.sbColorIntens->setValue(3);
      ui.cbTrkBB->setChecked(false);
      ptSizeScan = 2;
      ptSizeMovingTrk = 3;
    } break;
  }

  // Image control
  connect(ui.pbNext, SIGNAL(pressed()), this, SLOT(imgNext()) );
  connect(ui.pbPrev, SIGNAL(pressed()), this, SLOT(imgPrev()) );
  connect(ui.pbImgPlay, SIGNAL(toggled(bool)), this, SLOT(imgPlay(bool)));
  connect(ui.hsImageNb, SIGNAL(valueChanged(int)), this, SLOT(imgSelect(int)) );
  connect(ui.hsImageNb, SIGNAL(sliderMoved(int)), this, SLOT(updateSelectedImg(int)) );
  connect(ui.lwTrackIDs, SIGNAL(itemSelectionChanged()), this, SLOT(update3D()) );
  connect(ui.cbSelectAll, SIGNAL(toggled(bool)), this, SLOT(selectAll(bool)));
  connect(ui.pbSelectAll, SIGNAL(pressed()), this, SLOT(selectAll()) );
  connect(ui.pbSelectNone, SIGNAL(pressed()), this, SLOT(selectNone()) );
  connect(ui.pbSelectInvert, SIGNAL(pressed()), this, SLOT(selectInvert()) );

  // Draw Scan
  connect(ui.gbDrawImage, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.rbImgColorConst, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbImgColorHeight, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbImgColorDistance, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.sbColorIntens, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.pbImgPtPlus, SIGNAL(pressed()), this, SLOT(incImgPtSize()) );
  connect(ui.pbImgPtMinus, SIGNAL(pressed()), this, SLOT(decImgPtSize()) );

  // Draw World Model
  connect(ui.gbDrawWorld, SIGNAL(toggled(bool)), this, SLOT(emitStateChanged()) );
  connect(ui.pbWrldPtPlus, SIGNAL(pressed()), this, SLOT(incWrldPtSize()) );
  connect(ui.pbWrldPtMinus, SIGNAL(pressed()), this, SLOT(decWrldPtSize()) );
  connect(ui.dsbWorldHeight, SIGNAL(valueChanged(double)), this, SLOT(recolorWorldModel()) );
  connect(ui.dsbWorldRoll, SIGNAL(valueChanged(double)), this, SLOT(recolorWorldModel()) );
  connect(ui.dsbWorldPitch, SIGNAL(valueChanged(double)), this, SLOT(recolorWorldModel()) );
  connect(ui.rbWorldHeight, SIGNAL(toggled(bool)), this, SLOT(recolorWorldModel(bool)) );
  connect(ui.rbWorldGray, SIGNAL(toggled(bool)), this, SLOT(recolorWorldModel(bool)) );
  connect(ui.dsbShiftUp, SIGNAL(valueChanged(double)), this, SLOT(shiftupChanged()) );

  // Draw Tracks
  connect(ui.gbDrawTracks, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbTrkPointcloud, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbTrkBB, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbTrkCS, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.pbTrkPtPlus, SIGNAL(pressed()), this, SLOT(incTrkPtSize()) );
  connect(ui.pbTrkPtMinus, SIGNAL(pressed()), this, SLOT(decTrkPtSize()) );
  //connect(ui.pbTrkRecolor, SIGNAL(pressed()), this, SLOT() );
  connect(ui.rbTrkID, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkAge, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkUp, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkHeight, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkSpeed, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkNormal, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbTrkNConf, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );

  connect(ui.cbStatic, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.pbStatPtPlus, SIGNAL(pressed()), this, SLOT(incStatPtSize()) );
  connect(ui.pbStatPtMinus, SIGNAL(pressed()), this, SLOT(decStatPtSize()) );
  connect(ui.rbStatGray, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbStatHeight, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbStatNormal, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );
  connect(ui.rbStatNConf, SIGNAL(toggled(bool)), this, SLOT(update3D(bool)) );

  connect(ui.cbFilterAge, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFilterUp, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFilterPtCnt, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFilterDist, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFilterHeight, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.cbFilterVolume, SIGNAL(toggled(bool)), this, SLOT(update3D()) );
  connect(ui.sbFilterAge, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.sbFilterUp, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.sbFilterPtCnt, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.sbFilterDist, SIGNAL(valueChanged(int)), this, SLOT(update3D()) );
  connect(ui.dsbFilterHeight, SIGNAL(valueChanged(double)), this, SLOT(update3D()) );
  connect(ui.dsbFilterVolume, SIGNAL(valueChanged(double)), this, SLOT(update3D()) );

  imgSelect(ui.hsImageNb->value());

  recolorWorldModel();
  shiftupChanged();
  update3D();
}


VisualizerTracks::~VisualizerTracks()
{
  precacheThread->terminate();
  playTimer->stop();

  glDeleteLists(glListIndex, 2);

  BOOST_FOREACH(FrameData* fd,  loadedTracks)
    delete fd; // safe if NULL
  delete proj;
  delete worldModelRenderer;
}

void VisualizerTracks::paintGLOpaque()
{
  glPushMatrix();
  glTranslatef(0.0, 0.0, ui.dsbShiftUp->value());
  if (worldModelRenderer && (ui.gbDrawWorld->isChecked())) {
    glPushMatrix();
    if (currWorldTrack) {
      DCMatrix htm = glHtm2TrackCS(*currWorldTrack); // OpenGL needs column-major
      double *htmptr = htm.data().begin();
      glMultMatrixd(htmptr);
    }
    worldModelRenderer->render(ptSizeWorld);
    glPopMatrix();
  }
  glCallList(glListIndex+0);
  glPopMatrix();
}

void VisualizerTracks::paintGLTranslucent()
{
  glPushMatrix();
  glTranslatef(0.0, 0.0, ui.dsbShiftUp->value());
  glCallList(glListIndex+1);
  glPopMatrix();
}

void VisualizerTracks::loadWorldModelRenderer(string p3dFile)
{
  const unsigned int bytesPerPointP3D = 25; // 25 characters/bytes per point
  if (exists(path(p3dFile))) {
    unsigned long filesize = file_size(p3dFile);
    unsigned long expectedFilePts = filesize / bytesPerPointP3D;
    cout << "loading ~" << expectedFilePts << " points from file " << p3dFile << "..." << flush;
    delete worldModelRenderer;
    worldModelRenderer = new Gui3DQt::PointCloudRenderer();
    worldModelRenderer->reserve(expectedFilePts*1.1);
    ifstream in(p3dFile.c_str());
    if (!in.good())
      throw runtime_error("VisualizerTracks::loadWorldModelRenderer: problem opening the file for reading");
    Gui3DQt::PointCloudRenderer::GlVec3 point;
    Gui3DQt::PointCloudRenderer::GlCol3 color;
    while (!in.eof()) {
      in >> point.x;
      in >> point.y;
      in >> point.z;
      worldModelRenderer->push_back(point, color);
    }
    cout << "done" << endl;
  }
}

void VisualizerTracks::recalcWorldMinMaxZ()
{
  const float roll = ui.dsbWorldRoll->value()*M_PI/180.0f;
  const float pitch = ui.dsbWorldPitch->value()*M_PI/180.0f;
  worldMinZ = FLT_MAX;
  worldMaxZ = -FLT_MAX;
  float z;
  for (size_t i=0; i<worldModelRenderer->size(); ++i) {
    const Gui3DQt::PointCloudRenderer::GlVec3& pt = worldModelRenderer->pointAt(i);
    z = rotateZ(pt.x, pt.y, pt.z, roll, pitch);
    worldMinZ = min(worldMinZ,z);
    worldMaxZ = max(worldMaxZ,z);
  }
}

void VisualizerTracks::recolorWorldModel()
{
  if (!worldModelRenderer)
    return;
  recalcWorldMinMaxZ();
  const float sliceHeight = min((float)ui.dsbWorldHeight->value(),worldMaxZ-worldMinZ);
  const float roll = ui.dsbWorldRoll->value()*M_PI/180.0f;
  const float pitch = ui.dsbWorldPitch->value()*M_PI/180.0f;
  int r,g,b;
  for (unsigned int i=0; i<worldModelRenderer->size(); ++i) {
    const Gui3DQt::PointCloudRenderer::GlVec3 &pt = worldModelRenderer->pointAt(i);
    Gui3DQt::PointCloudRenderer::GlCol3 &cl = worldModelRenderer->colorAt(i);
    if (ui.rbWorldHeight->isChecked()) { // color slices of 2m repeatedly
      float z = rotateZ(pt.x, pt.y, pt.z, roll, pitch);
      int relZ = ((int)((z-worldMinZ)*255.0/sliceHeight)) % 255;
      HSV2RGB(relZ*360/255, 255, 255, r,g,b);
      cl = Gui3DQt::PointCloudRenderer::GlCol3(r, g, b);
    }
    if (ui.rbWorldGray->isChecked())
      cl = Gui3DQt::PointCloudRenderer::GlCol3(100, 100, 100);
  }
  emit stateChanged();
}

bool VisualizerTracks::skipTrack(const TrackData &tdata)
{
  QList<QListWidgetItem*> listItems;
  QListWidgetItem*        listWidget = NULL;
  if (!exitAfterLastFrame) {
    listItems = ui.lwTrackIDs->findItems( trackIDString(tdata), Qt::MatchExactly );
    //cout << "found matching " << listItems.size() << " items" << endl;
    assert( listItems.size() == 1 );
    listWidget = listItems[0];
  }
  if (listWidget && !listWidget->isSelected()) return true;
  if ((!ui.cbStatic->isChecked())      && (tdata.uid == 1)) return true;
  if ((ui.cbFilterAge->isChecked())    && ((int)tdata.age <= ui.sbFilterAge->value())) return true;
  if ((ui.cbFilterUp->isChecked())     && ((int)tdata.lastUp >= ui.sbFilterUp->value())) return true;
  if ((ui.cbFilterPtCnt->isChecked())  && ((int)tdata.nbPts <= ui.sbFilterPtCnt->value())) return true;
  if ((ui.cbFilterDist->isChecked())   && (norm(tdata.center) >= ui.sbFilterDist->value())) return true;
  if ((ui.cbFilterHeight->isChecked()) && ((tdata.bbExtensions[2] == 0)
                                           || (log(max(tdata.bbExtensions[0],tdata.bbExtensions[1])/tdata.bbExtensions[2]) >= ui.dsbFilterHeight->value()))) return true;
  if ((ui.cbFilterVolume->isChecked()) && ((tdata.bbExtensions[0]*tdata.bbExtensions[1]*tdata.bbExtensions[2]) >= ui.dsbFilterVolume->value())) return true;


  return false;
}

void VisualizerTracks::applyColor(const TrackData &tdata)
{
  if (tdata.uid == 1) {
    glColor3f(0.5, 0.5, 0.5); // gray
  } else {
    if (ui.rbTrkID->isChecked())
      glColor3f(tdata.colorR, tdata.colorG, tdata.colorB);
    if (ui.rbTrkAge->isChecked()) {
      float mult = max(0.0f,min(1.0f,(float)(tdata.age-3)/10.0f));
      //glColor3f(mult*tdata.colorR, mult*tdata.colorG, mult*tdata.colorB);
      glColorHSV(mult*120, 255, 255); // H = 0..360; S,V=0..255
    }
    if (ui.rbTrkUp->isChecked()) {
      float mult = 1.0f-min(1.0f,(float)tdata.lastUp/10.0f);
      //glColor3f(mult*tdata.colorR, mult*tdata.colorG, mult*tdata.colorB);
      glColorHSV(mult*120, 255, 255); // H = 0..360; S,V=0..255
    }
    if (ui.rbTrkSpeed->isChecked()) {
      float speed = sqrt(sqr(tdata.state[9])+sqr(tdata.state[10])+sqr(tdata.state[11])); // in m/s
      float mult = min(1.0f,speed/10.0f);
      glColorHSV(241-mult*240, 255, 255); // H = 0..360; S,V=0..255
    }
  }
}

void VisualizerTracks::update3D()
{
  if (listIsUpdating) return;

//  GLdouble clr[4];
//  glGetDoublev(GL_COLOR_CLEAR_VALUE, &(clr[0]));
//  bool invert = (clr[0]+clr[1]+clr[2])/3 > 0.5; // bright background -> 1, dark background -> 0

  trackMutex.lock();
  const FrameData *frame = loadedTracks[ui.hsImageNb->value()]; // assume that pointer is valid
  trackMutex.unlock();

  if ((frame->tracks.size() > 0) && (frame->tracks[0].uid == WORLD_TRACK_ID)) // tracks are sorted, so world track must be at first position
    currWorldTrack = &(frame->tracks[0]);
  else
    currWorldTrack = NULL;

  glNewList(glListIndex+0, GL_COMPILE); // replaces list (i.e. deletes old list before creating a new one)

  if (ui.gbDrawImage->isChecked()) {
    const LidarImage<double> *distance = &(frame->distanceImg);
    double i = (double)ui.sbColorIntens->value()/9.0;
    LidarImageVisualization::ColorMode3D mode = LidarImageVisualization::ColorFixed;
    if (ui.rbImgColorHeight->isChecked()) mode = LidarImageVisualization::ColorHeight;
    if (ui.rbImgColorDistance->isChecked()) mode = LidarImageVisualization::ColorDist;
    LidarImageVisualization::render3DPoints(frame->pointcloud, mode, ptSizeScan, i, i, i, distance, VISU_DIST_MIN, VISU_DIST_MAX);
  }

  if (ui.gbDrawTracks->isChecked()) {
    //QList<QListWidgetItem*> selTracks = ui.lwTrackIDs->selectedItems();
    BOOST_FOREACH(const TrackData &tdata, frame->tracks) { // loop over tracks
      if (skipTrack(tdata))
        continue;

      applyColor(tdata);
      if (ui.cbTrkPointcloud->isChecked()) {
        if (tdata.uid == 1) {
          glPointSize(ptSizeStaticTrk);
          glBegin(GL_POINTS);
          BOOST_FOREACH(const TrackData::Surface &s, tdata.ptsVIP)
            drawStaticSurface(s);
          BOOST_FOREACH(const TrackData::Surface &s, tdata.ptsADD)
            drawStaticSurface(s);
          glEnd();
        } else {
          glPointSize(ptSizeMovingTrk+1);
          glBegin(GL_POINTS);
          BOOST_FOREACH(const TrackData::Surface &s, tdata.ptsVIP)
            drawMovingSurface(s);
          glEnd();
          glPointSize(ptSizeMovingTrk);
          glBegin(GL_POINTS);
          BOOST_FOREACH(const TrackData::Surface &s, tdata.ptsADD)
            drawMovingSurface(s);
          glEnd();
        }
      }
      if (ui.cbTrkCS->isChecked()) {
        // change coordinate system to Track-CS
        glPushMatrix();
        DCMatrix htm2 = glHtm2TrackCS(tdata); // OpenGL needs column-major
        double *htm2p = htm2.data().begin();
        glMultMatrixd(htm2p);

        // draw coordinate frame
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
        glPopMatrix();
      }

    } // loop tracks
  } // draw tracks checked

  glEndList();

  glNewList(glListIndex+1, GL_COMPILE); // replaces list

  if ((ui.gbDrawTracks->isChecked()) && (ui.cbTrkBB->isChecked())) {
    GLUquadric *quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_LINE); //GLU_LINE,GLU_FILL
    glLineWidth(ptSizeMovingTrk);
    BOOST_FOREACH(const TrackData &tdata, frame->tracks) { // loop over tracks
      if (skipTrack(tdata))
        continue;
      if ((tdata.bbExtensions[0] == 0) || (tdata.bbExtensions[1] == 0) || (tdata.bbExtensions[2] == 0))
        continue;
      if (tdata.uid == 1) // don't render bounding box for world track
        continue;
      
      applyColor(tdata);

      glPushMatrix();
      glTranslatef(tdata.center[0], tdata.center[1], tdata.center[2]);
      glRotatef(tdata.bbAngle*180.0/M_PI, 0, 0, 1);
      glBegin(GL_LINES);
      glVertex3f(0,0,0);
      glVertex3f(1,0,0);
      glEnd();
      glScalef(max(0.01f,tdata.bbExtensions[0]), max(0.01f,tdata.bbExtensions[1]), max(0.01f,tdata.bbExtensions[2]));
      //glColor4f( 0.8, 0.8, 0.8, alpha);
      //Gui3DQt::Graphics::draw_cube_solid(-1,1,-1,1,-1,1);
      //glColor4f( 0.8, 0.8, 0.8, alpha/2); //const grey color for lines
      Gui3DQt::Graphics::draw_cube_cage(-1,1,-1,1,-1,1);
      glPopMatrix();
    }
    gluDeleteQuadric(quad);
  }

  glEndList();
  emit stateChanged();
}

DCMatrix VisualizerTracks::glHtm2TrackCS(const TrackData &tdata)
{
  DMatrix htm = YawPitchRollXYZ_2_HTM(tdata.state[2], tdata.state[1], tdata.state[0], tdata.state[3], tdata.state[4], tdata.state[5]);
//  cout << " state";
//  for (int i=0;i<6;++i)
//    cout << " " << tdata.state[i];
//  cout << " htm";
//  for (int i=0;i<4;++i)
//    for (int j=0;j<4;++j)
//      cout << " " << htm(i,j);
  return DCMatrix(htm); // OpenGL needs column-major
}

void VisualizerTracks::drawMovingSurface(const TrackData::Surface &s)
{
  drawSurface(s, ui.rbTrkHeight->isChecked(), ui.rbTrkNormal->isChecked(), ui.rbTrkNConf->isChecked());
}
void VisualizerTracks::drawStaticSurface(const TrackData::Surface &s)
{
  drawSurface(s, ui.rbStatHeight->isChecked(), ui.rbStatNormal->isChecked(), ui.rbStatNConf->isChecked());
}
void VisualizerTracks::drawSurface(const TrackData::Surface &s, bool height, bool normal, bool nConf)
{
  if (height) {
    float hrel = min(1.0f,max(0.0f,(s.pt[2]+DEFAULT_SLICE)/DEFAULT_SLICE));
    glColor3f(hrel, 0.2f, 1.0f-hrel); //0.5f - fabs(hrel-0.5f)
  }
  if (normal)
    glColor3f(fabs(s.nr[0]), fabs(s.nr[1]), fabs(s.nr[2]));
  if (nConf)
    glColor3f(1.0f-s.nc, s.nc, 0.5f - fabs(s.nc-0.5f));
  glVertex3f(s.pt[0], s.pt[1], s.pt[2]);
}

QString VisualizerTracks::trackIDString(const TrackData& track)
{
  QString str;
  str.setNum(track.uid);
  return str;
}

void VisualizerTracks::imgPlay(bool state)
{
  playTimer->setInterval(ui.sbPlayInterval->value());
  if (state)
    playTimer->start();
  else
    playTimer->stop();
}

void VisualizerTracks::imgNext()
{
  if (ui.hsImageNb->value() == ui.hsImageNb->maximum()) {
    playTimer->stop();
    if (exitAfterLastFrame)
      exit(0);
  } else {
    ui.hsImageNb->setValue(ui.hsImageNb->value()+1);
  }
}

void VisualizerTracks::imgPrev()
{
  ui.hsImageNb->setValue(ui.hsImageNb->value()-1);
}

void VisualizerTracks::imgSelect(int val)
{
  // load files if not yet done
  assertFrameData(val);
  updateSelectedImg(val);
  
  // update list-widget, keeping (de-)selected items (de-)selected
  trackMutex.lock();
  const FrameData *frame = loadedTracks[val];
  trackMutex.unlock();
  
  if (!exitAfterLastFrame) {
    listIsUpdating = true;
    // 1) save which items were selected
    QList<QListWidgetItem*> selTracks = ui.lwTrackIDs->selectedItems();
    QList<unsigned long long> selectedUids;
    BOOST_FOREACH(QListWidgetItem* wid, selTracks) {
      selectedUids.push_back(wid->text().toLongLong());
    }
    // 2) loop over tracks creating new list items if they don't exist;
    ui.lwTrackIDs->clear(); // should free up / delete all memory of the items
    BOOST_FOREACH(const TrackData &tdata, frame->tracks) { // loop over tracks
      QListWidgetItem *ni = new QListWidgetItem(trackIDString(tdata), ui.lwTrackIDs);
      if ((ui.cbSelectAll->isChecked()) || (selectedUids.contains(tdata.uid)))
        ni->setSelected(true);
    }
    listIsUpdating = false;
  }

  if (exitAfterLastFrame) {
    if (val>0) precacheThread->freeMem(val-1);
    if (val<(int)loadedTracks.size()-1) precacheThread->freeMem(val+1);
  } else {
    cout << endl << "starting precache thread" << flush;
    precacheThread->start();
  }

  update3D();
}

void VisualizerTracks::updateSelectedImg(int val)
{
  path fullfile(images[val]);
  ui.lImageName->setText(QString(fullfile.leaf().c_str()));
}

FrameData* VisualizerTracks::loadFrame(int idx)
{
  FrameData *newdata = new FrameData(images[idx], tracks[idx], tracks[idx]+".bb", *proj, loadPoints);
  //cout << "loaded frame " << images[idx] << " using " << newdata->getMemUsage()/1024 << " kB " << endl;
  return newdata;
}


void VisualizerTracks::assertFrameData(int val)
{
  trackMutex.lock();
  if (loadedTracks[val] == NULL) {
    loadedTracks[val] = loadFrame(val);
  }
  trackMutex.unlock();
}

void VisualizerTracks::selectAll(bool newState)
{
  for (int i=0; i<ui.lwTrackIDs->count(); ++i)
    ui.lwTrackIDs->item(i)->setSelected(newState);
}

void VisualizerTracks::selectInvert()
{
  for (int i=0; i<ui.lwTrackIDs->count(); ++i)
    ui.lwTrackIDs->item(i)->setSelected(!(ui.lwTrackIDs->item(i)->isSelected()));
}

void VisualizerTracks::shiftupChanged()
{
  if (passatVis) {
    cout << endl << "setting pose" << flush;
    passatVis->setPose(Gui3DQt::VisualizerPassat::translate_to_velobase_x,
                  Gui3DQt::VisualizerPassat::translate_to_velobase_y,
                  Gui3DQt::VisualizerPassat::translate_to_groundcenter_z - PASSAT_HEIGHT + ui.dsbShiftUp->value(), //translate_to_groundcenter_z
                  0);
  }
  emitStateChanged();
}
