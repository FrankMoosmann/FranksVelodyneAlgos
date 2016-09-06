#include "VisualizerObstacleTracking.hpp"

#include <sys/time.h>
#include <cfloat>
#include <iostream>
#include <sstream>
#include <fstream>
#include <QInputDialog>
#include <boost/bind.hpp>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "AlgoFrameFeatures.hpp"
#include "AlgoSegmentation.hpp"
#include "ParameterHeap.hpp"
#include "AlgoRegistration.hpp"

using namespace std;
using namespace matrixTools;
using namespace Gui3DQt;
namespace fs = boost::filesystem;

const double VISU_2D_DIST_MAX = 80.0; // in meter
const double VISU_2D_DDERIV_MAX = 10.0; // in meter
const double VISU_3D_SHIFTUP = 1.9f;


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////         Con-/Destructor           ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

VisualizerObstacleTracking::VisualizerObstacleTracking(ObstacleTracking &tracker_, int dummy, MNavWidget* glWid_, QWidget *parent)
  : Visualizer(parent)
  , tracker(tracker_)
  , glWid(glWid_)
  , vis()
  , trackMergeIdx(-1)
  , blockParamHeapWrite(false)
{
  (void)dummy;

  ui.setupUi(this);
  setParamsDefault();
  ui.sbBuffSize->setValue(tracker.bufferSize());
  ui.hsFrameSelect->setMinimum(-(tracker.bufferSize()-1));
  FrameSPtr frame = tracker.getFrame(0);
  ui.sbTrackBufferIdx->setMaximum(frame->sttValidationSteps); //sttValidationSteps refering to LTT

  glListIndex = glGenLists(1); // generate a display list
  connect(ui.hsFrameSelect, SIGNAL(valueChanged(int)), this, SLOT(render2D()) );
  connect(ui.hsFrameSelect, SIGNAL(valueChanged(int)), this, SLOT(render3D()) );

  paramsChanged();
  
  // 2D-Visualization:
  // unfortunately the following will execute render2D() twice when selection changes (deselect + select)
  connect(ui.gbVisu2DSelect, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuDistance, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuDistDiff, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuIntensity, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuHorizDeriv, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuVertDeriv, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuHorizConn, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuVertConn, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuHorizSegConn, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuVertSegConn, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuHorizSegThresh, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuVertSegThresh, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuHorizSegDecision, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuVertSegDecision, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuNormals, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuNormalConfidence, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuSegments, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuProjTracks, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuProjTrackCnt, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuProjTracksExt, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuSpeed, SIGNAL(toggled(bool)), this, SLOT(render2D()) );
  connect(ui.rbVisuPixMask, SIGNAL(toggled(bool)), this, SLOT(render2D()) );

  // 3D-Visualization:
  connect(ui.pbCenter3DPos, SIGNAL(pressed()), this, SLOT(center3DCamPos()) );
  connect(ui.cbVisuPassat3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuPassatCS, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuGrid3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuPoints3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuPointCovar3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuConn3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuNormals3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuNormalCovar3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuSegCrit3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuSegConn3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuSegDecision3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuSegments3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuFeatCorresp3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuTracks3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuTrackHist3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuTrackVariance3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuTrackProps3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuTrackMergeDec3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbVisuRegistration3D, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.cbRestrictVisu, SIGNAL(toggled(bool)), this, SLOT(render3D()) );
  connect(ui.pbPointPlus, SIGNAL(pressed()), this, SLOT(incPointSize()) );
  connect(ui.pbPointMinus, SIGNAL(pressed()), this, SLOT(decPointSize()) );
  connect(ui.spTrackHist, SIGNAL(valueChanged(int)), this, SLOT(paramsChanged()) );
  connect(ui.dsbHistAng, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.dsbHistPos, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.spTrackHist, SIGNAL(valueChanged(int)), this, SLOT(render3D()) );
  connect(ui.dsbHistAng, SIGNAL(valueChanged(double)), this, SLOT(render3D()) );
  connect(ui.dsbHistPos, SIGNAL(valueChanged(double)), this, SLOT(render3D()) );
  connect(ui.sbTrackBufferIdx, SIGNAL(valueChanged(int)), this, SLOT(render3D()) );
  connect(ui.sbTMMerge, SIGNAL(valueChanged(int)), this, SLOT(render3D()) );

  // changes of parameters:
  connect(ui.dsbConnMaxDist, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.dsbConnMaxDistRatio, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.dsbConvexThresh, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.dsbNormalThresh, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.dsbMinSegScore, SIGNAL(valueChanged(double)), this, SLOT(paramsChanged()) );
  connect(ui.sbMinSegSizePt, SIGNAL(valueChanged(int)), this, SLOT(paramsChanged()) );
  connect(ui.sbMinSegSizeCm, SIGNAL(valueChanged(int)), this, SLOT(paramsChanged()) );
  connect(ui.cbSegUseTrackIDs, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbDynamicMapping, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbSimpleMapping, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbOverwriteMapping, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbICPDistDiffWeight, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbICPOcclWeight, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );
  connect(ui.cbRegWithFeatMatch, SIGNAL(toggled(bool)), this, SLOT(paramsChanged()) );

  // buttons:
  connect(ui.pbResetAll, SIGNAL(pressed()), this, SLOT(resetBuffer()) );
  connect(ui.pbProcFrame, SIGNAL(pressed()), this, SLOT(processFrame()) );
  connect(ui.pbProcAllFrames, SIGNAL(pressed()), this, SLOT(processAllFrames()) );
  connect(ui.pbSetPVec1, SIGNAL(pressed()), this, SLOT(setPVec1()) );
  connect(ui.pbWriteDebugImages, SIGNAL(pressed()), this, SLOT(writeDebugImages()) );
  connect(ui.pbCalculateFeatures, SIGNAL(pressed()), this, SLOT(generateFeatures()) );
  connect(ui.pbCalculateSegments, SIGNAL(pressed()), this, SLOT(generateSegments()) );
  connect(ui.pbRecolorSegments, SIGNAL(pressed()), this, SLOT(generateSegmentColors()) );
  connect(ui.pbRecolorAllSegments, SIGNAL(pressed()), this, SLOT(generateAllSegmentColors()) );
  connect(ui.pbCalculateRegistration, SIGNAL(pressed()), this, SLOT(generateRegistration()) );
  connect(ui.pbCreateTracks, SIGNAL(pressed()), this, SLOT(generateTracks()) );
  connect(ui.pbTMInit, SIGNAL(pressed()), this, SLOT(trackMergeInit()) );
  connect(ui.pbTMMerge, SIGNAL(pressed()), this, SLOT(trackMergeMerge()) );
  connect(ui.pbTMKeep, SIGNAL(pressed()), this, SLOT(trackMergeKeep()) );
  connect(ui.pbTMIgnore, SIGNAL(pressed()), this, SLOT(trackMergeIgnore()) );
  connect(ui.pbTMBack, SIGNAL(pressed()), this, SLOT(trackMergeBack()) );
  connect(ui.pbTMPlus, SIGNAL(pressed()), this, SLOT(trackMergePlus()) );

  connect(ui.pbRestartReg, SIGNAL(pressed()), this, SLOT(resetSegReg()) );
  connect(ui.pbRegisterCurSeg, SIGNAL(pressed()), this, SLOT(registerCurrSeg()) );
  connect(ui.pbRegFinishLayer, SIGNAL(pressed()), this, SLOT(registerFinishLayer()) );
  connect(ui.pbRegisterFinish, SIGNAL(pressed()), this, SLOT(registerFinish()) );

  connect(ui.pbSaveICPConfig, SIGNAL(pressed()), this, SLOT(saveRegistrConfig()) );
  connect(ui.pbSegSave, SIGNAL(pressed()), this, SLOT(saveSegments()) );
  connect(ui.pbSegLoad, SIGNAL(pressed()), this, SLOT(loadSegments()) );

  ui.pbTMMerge->setShortcut(tr("Ctrl+M"));
  ui.pbTMKeep->setShortcut(tr("Ctrl+K"));
  ui.pbTMIgnore->setShortcut(tr("Ctrl+I"));
  ui.pbTMPlus->setShortcut(tr("Ctrl+U"));
  ui.pbTMBack->setShortcut(tr("Ctrl+J"));
  //ui.sbTMMerge->setShortcut();
  QAction *sbInc = new QAction(this);
  sbInc->setShortcut(tr("Ctrl+O"));
  connect(sbInc, SIGNAL(triggered(bool)), this, SLOT(trackMergeSBInc()) );
  QAction *sbDec = new QAction(this);
  sbDec->setShortcut(tr("Ctrl+L"));
  connect(sbDec, SIGNAL(triggered(bool)), this, SLOT(trackMergeSBDec()) );
  ui.gbTrackCreation->addAction(sbInc);
  ui.gbTrackCreation->addAction(sbDec);


  render2D();
  render3D();
}

VisualizerObstacleTracking::~VisualizerObstacleTracking()
{
  glDeleteLists(glListIndex, 1);
}


void write(const DVector &v, ostream &out) {
  out << v(0) << "," << v(1) << "," << v(2) << endl;
}
void write(const DCMatrix &m, unsigned int col, ostream &out) {
  out << m(0,col) << "," << m(1,col) << "," << m(2,col) << endl;
}
void writeTrans(const DMatrix &R, const DVector &t, ostream &out) {
  double yawRAD, pitchRAD, rollRAD, tx, ty, tz;
  HomogeneousTransformationMatrix::Rt_2_YawPitchRollXYZ( R, t, yawRAD, pitchRAD, rollRAD, tx, ty, tz);
  out << yawRAD<< "," << pitchRAD<< "," << rollRAD << "," << tx << "," << ty << "," << tz << endl;
}
void VisualizerObstacleTracking::saveRegistrConfig() {
  FrameSPtr world = tracker.registrator->getFrame();
  int hsize,vsize;
  world->point3D.getSize(hsize,vsize);
  TrackRegistration::SubTrackSPtr model = *(tracker.registrator->nextTrackToRegister);
  unsigned int nbPts = model->pointsTrackCS.size2();

  cout << endl << "saving model of size " << nbPts << "..." << flush;
  DVector mCenter = model->tInit;
//  DVector mCenter = ublas::zero_vector<double>(3);
//  for (unsigned int row=0; row<nbPts; ++row) {
//    mCenter(0) += model->pointsEgoCS(0,row);
//    mCenter(1) += model->pointsEgoCS(1,row);
//    mCenter(2) += model->pointsEgoCS(2,row);
//  }
//  mCenter /= (double)nbPts;

  ostringstream ossWPt, ossWNr, ossWNC, ossWNS;
  ostringstream ossMPt, ossMNr, ossMNC, ossMNS;
  // world:
  for (int h=0; h<hsize; ++h)
    for (int v=0; v<vsize; ++v)
      if (norm_1(world->normal3D.get(h,v)) > 0.1) // normal valid -> point valid
        if (norm_2(world->point3D.get(h,v) - mCenter) < 40) {
          //cout << "." << flush;
          write(world->point3D.get(h,v), ossWPt);
          write(world->normal3D.get(h,v), ossWNr);
          ossWNC << world->normalConfidence.get(h,v) << endl;
          ossWNS << world->normalStdDevRAD.get(h,v) << endl;
        }
  string s1 = ossWPt.str();
  cout << "world points converted into " << s1.length() << " characters" << flush;

  // model:
  for (unsigned int row=0; row<nbPts; ++row) {
    write(model->pointsTrackCS, row, ossMPt);
    write(model->pointNormalsTrackCS, row, ossMNr);
    ossMNC << model->normalConfidence(row) << endl;
    ossMNS << model->normalStdDevRAD(row) << endl;
  }

  ofstream outfile("ICPTestSetting.txt", ios::out);
  outfile << ossWPt.str() << endl; // world points
  outfile << ossWNr.str() << endl; // world normals
  outfile << ossMPt.str() << endl; // model points
  outfile << ossMNr.str() << endl; // model normals

  outfile << ossWNC.str() << endl; // world normal confidence
  outfile << ossWNS.str() << endl; // normal std-dev
  outfile << ossMNC.str() << endl; // normal normal confidence
  outfile << ossMNS.str() << endl; // normal normal std-dev

  // model transformation:
  writeTrans(model->RInit, model->tInit, outfile);
  outfile << endl;
  cout << endl << "successfully saved";
}

void VisualizerObstacleTracking::saveSegments()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) currRow = 0; // if no row is selected, use current frame
  FrameSPtr frame = tracker.getFrame(currRow);
  saveToFile(frame->segments, "segments.png");
  PngDistanceImage dstImg(frame->distance.getHorizSize(), frame->distance.getVertSize());
  dstImg.setDistances(frame->distance.getAdr(0,0), 2);
  dstImg.save("distance.png");
}

void VisualizerObstacleTracking::loadSegments()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) currRow = 0; // if no row is selected, use current frame
  FrameSPtr frame = tracker.getFrame(currRow);
  loadFromFile(frame->segments, "segments.png");
  recolorSegments(frame);
  render2D();
  render3D();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////          Visualization            ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int VisualizerObstacleTracking::getSelectedFrameNumber() const {
  return -ui.hsFrameSelect->value();
}

void VisualizerObstacleTracking::setSelectedFrameNumber(int row) {
  setMergingActive(false);
  ui.hsFrameSelect->setValue(-row);
}

void VisualizerObstacleTracking::paintGLOpaque()
{
  glCallList(glListIndex+0);
}

void VisualizerObstacleTracking::paintGLTranslucent()
{
}


// functions for creating GraphML-File from TrackIndex-TrackIndex associations
void trackindexKeySpecifier(std::ofstream &of)
{
  of << "  <key id=\"d3\" for=\"node\" yfiles.type=\"nodegraphics\"/>" << endl;
  of << "  <key id=\"d5\" for=\"node\" attr.name=\"description\" attr.type=\"string\"/>" << endl;
  of << "  <key id=\"d6\" for=\"edge\" yfiles.type=\"edgegraphics\"/>" << endl;
  of << "  <key id=\"d8\" for=\"edge\" attr.name=\"description\" attr.type=\"string\"/>" << endl;
}
void trackindexVisitor_base(const PointCloudTrack::SPtr track, std::ofstream &of, float rowOffset)
{
  const float colmult = 1;
  const float rowmult = 1;
  float colC, rowC;
  track->getCenter(colC, rowC);
  float size = (sqrt((double)track->getPointCount()/M_PI)); // area -> radius -> downscale
  int r = (int)(track->getColorR()*255.0);
  int g = (int)(track->getColorG()*255.0);
  int b = (int)(track->getColorB()*255.0);
  of << "      <data key=\"d5\"/>" << endl;
  of << "      <data key=\"d3\">" << endl;
  of << "        <y:ShapeNode>" << endl;
  of << "          <y:Geometry height=\"" << size << "\" width=\"" << size << "\" x=\"" << colC*colmult << "\" y=\"" << rowC*rowmult+rowOffset << "\"/>" << endl;
  of << "          <y:Fill color=\"#" << boost::format("%|02X|%|02X|%|02X|") % r % g % b << "\" transparent=\"false\"/>" << endl;
  of << "          <y:BorderStyle color=\"#000000\" type=\"line\" width=\"0.5\"/>" << endl;
//  of << "          <y:NodeLabel alignment=\"center\" autoSizePolicy=\"content\" fontFamily=\"Dialog\" fontSize=\"12\" fontStyle=\"plain\" hasBackgroundColor=\"false\" hasLineColor=\"false\" height=\"17.96875\" modelName=\"internal\" modelPosition=\"c\" textColor=\"#000000\" visible=\"true\" width=\"172.533203125\" x=\"-71.2666015625\" y=\"6.015625\">A1</y:NodeLabel>" << endl;
  of << "          <y:Shape type=\"circle\"/>" << endl;
  of << "        </y:ShapeNode>" << endl;
  of << "      </data>" << endl;
}
void trackindexVisitorL1(const LidarFrame::TrackIndex tIdx, std::ofstream &of, LidarFrame::TrackSPtrVec trackVec)
{
  PointCloudTrack::SPtr track = trackVec[tIdx];
  trackindexVisitor_base(track, of, 0.0);
}
void trackindexVisitorL2(const LidarFrame::TrackIndex tIdx, std::ofstream &of, LidarFrame::TrackSPtrVec trackVec)
{
  PointCloudTrack::SPtr track = trackVec[tIdx];
  trackindexVisitor_base(track, of, 200.0);
}
void trackindexVisitorL2L(const PointCloudTrack::IdT tId, std::ofstream &of, LidarFrame::TrackSPtrList trackList)
{
  BOOST_FOREACH(PointCloudTrack::SPtr track, trackList) {
    if (track->getUID() == tId) {
      trackindexVisitor_base(track, of, 200.0);
      break;
    }
  }
}
void trackindexConnVisitor(const LidarFrame::TrackIndex tIdx1, const LidarFrame::TrackIndex tIdx2, unsigned int val, std::ofstream &of)
{
  (void)tIdx1;
  (void)tIdx2;
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
void saveTrackLinkageGraph(LidarFrame::TrackLinkageGraphSPtr &graph, LidarFrame::TrackSPtrVec tracks1, LidarFrame::TrackSPtrVec tracks2, string filename)
{
  cerr << "saving track-links to file..." << flush;
  graph->saveToGraphML(filename, boost::bind(&trackindexKeySpecifier, _1)
                               , boost::bind(&trackindexVisitorL1, _1, _2, tracks1)
                               , boost::bind(&trackindexVisitorL2, _1, _2, tracks2)
                               , boost::bind(&trackindexConnVisitor, _1, _2, _3, _4)
  );
  cerr << "done" << flush;
}
void saveTrackLinkageGraph(LidarFrame::LTTrackLinkageGraphSPtr &graph, LidarFrame::TrackSPtrVec tracks, LidarFrame::TrackSPtrList lttTracks, string filename)
{
  cerr << "saving track-links to file..." << flush;
  graph->saveToGraphML(filename, boost::bind(&trackindexKeySpecifier, _1)
                               , boost::bind(&trackindexVisitorL1, _1, _2, tracks)
                               , boost::bind(&trackindexVisitorL2L, _1, _2, lttTracks)
                               , boost::bind(&trackindexConnVisitor, _1, _2, _3, _4)
  );
  cerr << "done" << flush;
}

void VisualizerObstacleTracking::writeDebugImages() {
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) currRow = 0; // if no row is selected, use current frame
  FrameSPtr frame = tracker.getFrame(currRow);

  QString imgBaseName = ui.leImgBaseName->text();
  imgBaseName.append('/');
  QString imgNbBaseName;
  QString currImgName;
  unsigned int frameCounter = 1;
  do {
    string number = (boost::format("%1$05d") % frameCounter).str();
    imgNbBaseName = imgBaseName;
    imgNbBaseName.append(number.c_str());
    currImgName = imgNbBaseName;
    currImgName.append("_dst.png");
    frameCounter++;
  } while (fs::exists(currImgName.toStdString()));
  cout << "store frame as " << imgNbBaseName.toStdString() << endl;

  QImage img;

  vis.renderColored(frame->distance, img, 1.5, VISU_2D_DIST_MAX); // Distance
  currImgName = imgNbBaseName; currImgName.append("_dst.png");
  img.save(currImgName);
  if (glWid) {
    img = glWid->grabFrameBuffer(); // bool withAlpha = false
    currImgName = imgNbBaseName; currImgName.append("_3d.png");
    img.save(currImgName);
  }
//  renderTrackProjection(frame->projTracks, frame, frame->projSttIdx, img);
//  currImgName = imgNbBaseName; currImgName.append("_prj.png");
//  img.save(currImgName);
//  vis.renderTrackProjectionCnt(frame, frame->projTrackCnt, img);
//  currImgName = imgNbBaseName; currImgName.append("_prjcnt.png");
//  img.save(currImgName);
//  renderTrackProjection(frame->extProjTracks, frame, frame->projSttIdx, img);
//  currImgName = imgNbBaseName; currImgName.append("_prj_ext.png");
//  img.save(currImgName);
//  renderTrackProjection(frame->projLttTracks, frame, frame->sttValidationSteps, img);
//  currImgName = imgNbBaseName; currImgName.append("_prj_ltt.png");
//  img.save(currImgName);
//  bool useSegPrj = ui.cbSegUseTrackIDs->isChecked();
//  ui.cbSegUseTrackIDs->setChecked(false);
//  generateSegments();
//  vis.renderSegments(frame->segments, img);
//  currImgName = imgNbBaseName; currImgName.append("_seg.png");
//  img.save(currImgName);
//  vis.renderMono(frame->segmentGrowDecisionH, img);
//  currImgName = imgNbBaseName; currImgName.append("_seg_decH.png");
//  img.save(currImgName);
//  vis.renderMono(frame->segmentGrowDecisionV, img);
//  currImgName = imgNbBaseName; currImgName.append("_seg_decV.png");
//  img.save(currImgName);
//  ui.cbSegUseTrackIDs->setChecked(true);
//  generateSegments();
//  vis.renderSegments(frame->segments, img);
//  currImgName = imgNbBaseName; currImgName.append("_segprj.png");
//  img.save(currImgName);
//  vis.renderMono(frame->segmentGrowDecisionH, img);
//  currImgName = imgNbBaseName; currImgName.append("_segprj_decH.png");
//  img.save(currImgName);
//  vis.renderMono(frame->segmentGrowDecisionV, img);
//  currImgName = imgNbBaseName; currImgName.append("_segprj_decV.png");
//  img.save(currImgName);
//  ui.cbSegUseTrackIDs->setChecked(useSegPrj);
//    vis.renderTranslationNorm(frame, img);
//    currImgName = imgNbBaseName; currImgName.append("_speed.png");
//    img.save(currImgName);
//  currImgName = imgNbBaseName; currImgName.append("_lnk2t-1.graphml");
//  saveTrackLinkageGraph(frame->shortTermTracks[0].linksT1, frame->shortTermTracks[0].tracks, frame->shortTermTracks[1].tracks, currImgName.toStdString());
//  currImgName = imgNbBaseName; currImgName.append("_lnk2ltt.graphml");
//  saveTrackLinkageGraph(frame->shortTermTracks[0].lttLinks, frame->shortTermTracks[0].tracks, frame->longTermTracks, currImgName.toStdString());
}

void VisualizerObstacleTracking::render2D()
{
  if (!ui.gbVisu2DSelect->isChecked())
    return;

  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  FrameSPtr frame = tracker.getFrame(currRow);
//  FrameSPtr nframe;
  QImage img;
  if (ui.rbVisuDistance->isChecked())
    vis.renderColored(frame->distance, img, 1.5, VISU_2D_DIST_MAX);
  if (ui.rbVisuDistDiff->isChecked())
    vis.renderGrey(frame->distanceDiffLastFrame, img, 0.0, 1.0, true);
  if (ui.rbVisuIntensity->isChecked())
    vis.renderGrey(frame->intensity, img, -0.2, 1.0); // -0.2 will scale rendering from 0.2..1. Invalid values will be 0.0 (thus, recognizable)
  if (ui.rbVisuHorizDeriv->isChecked())
    vis.renderGrey(frame->distDerivativeH, img, -VISU_2D_DDERIV_MAX, VISU_2D_DDERIV_MAX);
  if (ui.rbVisuVertDeriv->isChecked())
    vis.renderGrey(frame->distDerivativeV, img, -VISU_2D_DDERIV_MAX, VISU_2D_DDERIV_MAX);
  if (ui.rbVisuHorizConn->isChecked())
    vis.renderGrey(frame->connectivityH, img, -0.2, 1.0);
  if (ui.rbVisuVertConn->isChecked())
    vis.renderGrey(frame->connectivityV, img, -0.2, 1.0);
  if (ui.rbVisuHorizSegConn->isChecked())
    vis.renderGrey(frame->segmentConnectH, img, -0.2, 1.0);
  if (ui.rbVisuVertSegConn->isChecked())
    vis.renderGrey(frame->segmentConnectV, img, -0.2, 1.0);
  if (ui.rbVisuHorizSegThresh->isChecked())
    vis.renderGrey(frame->segmentThreshH, img, -0.2, 1.0);
  if (ui.rbVisuVertSegThresh->isChecked())
    vis.renderGrey(frame->segmentThreshV, img, -0.2, 1.0);
  if (ui.rbVisuHorizSegDecision->isChecked())
    vis.renderMono(frame->segmentGrowDecisionH, img);
  if (ui.rbVisuVertSegDecision->isChecked())
    vis.renderMono(frame->segmentGrowDecisionV, img);
  if (ui.rbVisuNormals->isChecked())
    vis.renderColoredNormals(frame->normal3D, frame->normalConfidence, img);
  if (ui.rbVisuNormalConfidence->isChecked())
    vis.renderColored(frame->normalConfidence, img, 0.000001, 1.0, false, 0.0, LidarImageVisualization::GREENRED);
  if (ui.rbVisuSegments->isChecked())
    vis.renderSegments(frame->segments, img);
  if (ui.rbVisuProjTracks->isChecked())
    renderTrackProjection(frame->projTracks, frame, frame->projSttIdx, img);
  if (ui.rbVisuProjTrackCnt->isChecked())
    vis.renderTrackProjectionCnt(frame->projTrackCnt, img);
  if ((ui.rbVisuProjTracksExt->isChecked()) && (frame->projSttIdx < frame->sttValidationSteps))
    renderTrackProjection(frame->extProjTracks, frame, frame->projSttIdx, img);
  if (ui.rbVisuPixMask->isChecked())
    vis.renderMono(frame->regAvailPixels, img);
  if (ui.rbVisuSpeed->isChecked())
    vis.renderTranslationNorm(frame, img);
  emit redraw2D(img);
}

void VisualizerObstacleTracking::render3D()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected

  FrameSPtr frame = tracker.getFrame(currRow);

  glDeleteLists(glListIndex, 1);
  glListIndex = glGenLists(1); // generate a display list
  glNewList(glListIndex, GL_COMPILE);
  GLdouble clr[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, &(clr[0]));
  bool invert = (clr[0]+clr[1]+clr[2])/3 > 0.5; // bright background -> 1, dark background -> 0

  if (ui.cbVisuGrid3D->isChecked())                                       // render grid (don't shift up)
    vis.render3DGrid();

  glPushMatrix();
  glTranslatef(0.0, 0.0, VISU_3D_SHIFTUP); // Point-shift-up

  if (ui.cbVisuPassat3D->isChecked())                                      // render Passat
    vis.render3DCarPos(ui.cbVisuPassatCS->isChecked());
  if (ui.cbVisuPoints3D->isChecked())                                     // render 3D Points
    vis.render3DPoints(frame, ui.cbVisuPointCovar3D->isChecked(), invert);
  if (ui.cbVisuNormals3D->isChecked())                                    // render Normals
    vis.render3DNormals(frame, ui.cbVisuNormalCovar3D->isChecked(), invert);
  if (ui.cbVisuConn3D->isChecked())                                       // render ConnectivityGraph
    vis.render3DConnections(frame->point3D, frame->connectivityH, frame->connectivityV);
  if (ui.cbVisuSegCrit3D->isChecked())                                    // render SegmentConnectivityGraph
    vis.render3DConnections(frame->point3D, frame->segmentCritR, frame->segmentCritD, &frame->segmentCritRU, &frame->segmentCritRD);
  if (ui.cbVisuFeatCorresp3D->isChecked()) // HACK, remove!
    vis.render3DConnections(frame->point3D, frame->segmentCritR, frame->segmentCritD);
  if (ui.cbVisuSegConn3D->isChecked())                                    // render SegmentConnectivityDecision
    vis.render3DConnections(frame->point3D, frame->segmentConnectH, frame->segmentConnectV);
  if (ui.cbVisuSegDecision3D->isChecked())                                    // render SegmentConnectivityDecision
    vis.render3DSegmentationLinks(frame);
  if (ui.cbVisuSegments3D->isChecked())                                   // render segments as 3D Points
    vis.render3DSegments(frame, invert);
//  if (ui.cbVisuFeatCorresp3D->isChecked()) {                              // render Feature Correspondences
//    try {
//      FrameSPtr lframe = tracker.getFrame(currRow+1);
//      vis.render3DCorrespondences(lframe, frame, invert);
//    } catch (std::range_error) {
//    }
//  }
  if (trackMergeIdx >= 0) {   // render track and link-based registration results to decide about merging
    // calculate and cout errors of various transformations (induced by linked tracks and itself)
    // use stt[STTVS-1].links to determine associated tracks to stt[STTVS]
    LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
    PointCloudTrack::SPtr track = buffer.tracks[trackMergeIdx];

    // 1st) render track in current state
    vis.render3DTrack(track.get(), true, 1.0, 0.2, 0.2); // render track in red

    // 2nd) render associated track
    LidarFrame::TrackLinks links = frame->getMergedSTTLinks(frame->sttValidationSteps-1, trackMergeIdx, true);
    int selectedLink = ui.sbTMMerge->value();
    ui.pbTMKeep->setEnabled(track->getLastUpdateCounter() == 0);
    if ((selectedLink >= 0) && (selectedLink < (int)links.size()) && (links[selectedLink].track.get() != NULL)) {
      PointCloudTrack::SPtr linkedTrack = links[selectedLink].track;
      unsigned int          linkStrength = links[selectedLink].linkStrength;
      if (track->getAge()+1 != track->histAllStates.size())
        throw std::runtime_error("render3D(): track-history-size doesn't match age");
      if (linkedTrack->getAge()+1 != linkedTrack->histAllStates.size())
        throw std::runtime_error("render3D(): tc-history-size doesn't match age");
      if (linkedTrack->getAge() <= track->getAge())
        throw std::runtime_error("render3D(): tc-history-size < track-history-size");

      // print current state of track
      DVector v1 = track->getState();
      DVector v2 = track->histAllStates.back().x;
      if (fabs(ublas::norm_2(v1-v2)) > 0.00001)
        throw std::runtime_error("render3D(): states don't match");
      tracker.registrator->usePixMask(false);
      pair<double,double> err = tracker.registrator->calculateError(track, track->getState());
      cout << endl << " Err [" << track->getUID() << "]:       \t" << err.first << "\t" << err.second
           << " up:" << track->getLastUpdateCounter()
           << " age:" << track->getAge()
           << " state: " << track->getState();

      // print properties of link
      ui.pbTMMerge->setEnabled(track->getMergeTrackIcpSuccess(linkedTrack.get()));
      // assumption: linked tracks from higher buffer-level were already registered
      // calculate move of "tc" within the same period of time, thus -p0+pcurr
      DVector origstate = DVectorRange(track->histAllStates.front().x, ublas::range(0,6));
      DVector tLinkedNow, t2linked;
      DMatrix RLinkedNow, R2linked;
      tLinkedNow = linkedTrack->histAllStates.back().t;
      RLinkedNow = linkedTrack->histAllStates.back().R;
      track->getRt2OtherTrkCS(R2linked, t2linked, *linkedTrack, track->getAge(), track->getAge());
      DVector tcomp; DMatrix Rcomp, Hcomp;
      Rcomp = ublas::prod(RLinkedNow,R2linked);
      tcomp = ublas::prod(RLinkedNow,t2linked) + tLinkedNow;
      Hcomp = Rt_2_HTM(Rcomp, tcomp);
      DVector xcomp; PointCloudTrack::Rt2state_T2E(Rcomp, tcomp, xcomp);
      bool linkICPSuccess = track->getMergeTrackIcpSuccess(linkedTrack.get());
      bool mergeworld = false; // associated track is rendered anyway
      float rg = linkICPSuccess ? .1 : .9;
      float b  = mergeworld     ? .1 : .9;
      vis.render3DTrack(linkedTrack.get(), false, 0.5, 0.5, 0.5); // render linked track in gray
      vis.render3DTrack(track.get(), true, rg, rg, b, &Hcomp); // render track with new transform in blue/dark/bright
      pair<double,double> errLink = tracker.registrator->calculateError(track, xcomp);
      cout << endl << "-Err [" << linkedTrack->getUID() << "](|" << linkStrength << "|/" << links.size() << "):\t" << errLink.first << "\t" << errLink.second
           << " up:" << linkedTrack->getLastUpdateCounter()
           << " age:" << linkedTrack->getAge()
           << " statediff: " << ublas::norm_2(origstate-xcomp) << " state: " << xcomp << " orig: " << origstate;
      tracker.registrator->usePixMask(true);
    } else {
      ui.pbTMMerge->setEnabled(false);
    }
    cout << flush;
  }
  if (ui.cbVisuTrackMergeDec3D->isChecked()) {                                // render track-merge-decision
    LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
    BOOST_FOREACH(PointCloudTrack::SPtr track, buffer.tracks) {
      PointCloudTrack::SPtr mergeInto = track->trackMaxMergeScore;
      if (track->mergeDecision == PointCloudTrack::MDKeep)
        vis.render3DTrack(track.get(), true, 1., .2, .2, (DMatrix*)NULL); // render in red in current state
      else if (track->trackMaxMergeScore.get() != NULL) { // will be merged or ignored
        DVector tLinkedNow = mergeInto->histAllStates.back().t;
        DMatrix RLinkedNow = mergeInto->histAllStates.back().R;
        DMatrix R2linked; DVector t2linked;
        track->getRt2OtherTrkCS(R2linked, t2linked, *mergeInto, track->getAge(), track->getAge());
        DVector tcomp; DMatrix Rcomp, Hcomp;
        Rcomp = ublas::prod(RLinkedNow,R2linked);
        tcomp = ublas::prod(RLinkedNow,t2linked) + tLinkedNow;
        Hcomp = Rt_2_HTM(Rcomp, tcomp);
        //bool linkICPSuccess = track->getMergeTrackIcpSuccess();
        bool linkICPSuccess = (track->mergeDecision == PointCloudTrack::MDMerge);
        bool mergeworld = (mergeInto == frame->getLttWorldTrack());
        float rg = linkICPSuccess ? .1 : .9;
        float b  = mergeworld     ? .1 : .9;
        vis.render3DTrack(track.get(), true, rg, rg, b, &Hcomp); // render in gray/blue/ in state of merge-into-track
      }
    }
  }
  if (ui.cbVisuTracks3D->isChecked()) {     // render Tracks
    if ((ui.cbRestrictVisu->isChecked()) && (tracker.registrator.get() != NULL)) { // render current track only
      if (tracker.registrator->currTrackRegistered != tracker.registrator->subsampledTracks.end()) {
        PointCloudTrack::SPtr track = (**(tracker.registrator->currTrackRegistered)).track;
        vis.render3DTrack(track.get());
        if (ui.cbCoutTrackDetails->isChecked()) {
          unsigned int pointCount;
          double minDist, maxDist;
          track->getStatistics(pointCount, minDist, maxDist);
          cout << endl << "     " << frame->getTime() << " track[" << track->getUID() << "](" << pointCount << "): minDist=" << minDist << ", maxDist=" << maxDist << flush;
        }
      }
      // also render next track as preview
      if (tracker.registrator->nextTrackToRegister != tracker.registrator->subsampledTracks.end()) {
        PointCloudTrack::SPtr track = (**(tracker.registrator->nextTrackToRegister)).track;
        vis.render3DTrack(track.get(), false);
      }
    } else {
      callForEachTrack(frame, boost::bind(&AlgoVisualization::render3DTrack, &vis, _1, true, -1., -1., -1., (DMatrix*)NULL));
    }
  }
  if ((ui.cbVisuRegistration3D->isChecked()) && (tracker.registrator.get() != NULL)) {   // render Registration Result of current "Registrator" track
    if (ui.cbRestrictVisu->isChecked()) {
      if (tracker.registrator->currTrackRegistered != tracker.registrator->subsampledTracks.end()) {
        vis.render3DRegistrationResult(tracker.registrator->currTrackRegistered);
        if (ui.cbCoutTrackDetails->isChecked()) {
          TrackRegistration::SubTrackSPtr subtrack = *tracker.registrator->currTrackRegistered;
          PointCloudTrack::SPtr track = subtrack->track;
          cout << endl << "     track[" << track->getUID() << "]@l" << subtrack->trackBufferIdx << ":";
          cout << endl << " PRED: x=" << track->histXPredict << ", P=" << track->histPPredict;
          cout << endl << " ICP: ";
          if (track->histRegUsedFeatureMatching) cout << "(featOK) ";
          BOOST_FOREACH(double err, track->histRegICPError) cout << " d=" << err;
          cout << endl << " UP: zExp=" << track->histZExp << ", z=" << track->histZ << ", R=" << track->histR << ", KGain=" << track->histKGain;
          cout << endl << "     x=" << track->getState() << ", P=" << track->getCovar();
        }
      }
    } else {
      for (BOOST_AUTO(tracki, tracker.registrator->subsampledTracks.begin()); tracki!=tracker.registrator->subsampledTracks.end(); ++tracki) {
        vis.render3DRegistrationResult(tracki);
      }
    }
  }
  if (ui.cbVisuTrackHist3D->isChecked()) {  // render Track History
    if ((ui.cbRestrictVisu->isChecked()) && (tracker.registrator.get() != NULL)) {
      if (tracker.registrator->currTrackRegistered != tracker.registrator->subsampledTracks.end()) {
        PointCloudTrack::SPtr track = (**(tracker.registrator->currTrackRegistered)).track;
        vis.render3DTrackHistory(track.get());
      }
    } else {
      callForEachTrack(frame, boost::bind(&AlgoVisualization::render3DTrackHistory, &vis, _1), true);
    }
  }
  if (ui.cbVisuTrackVariance3D->isChecked()) {                               // render Track Variance of current "Registrator" track
    if ((ui.cbRestrictVisu->isChecked()) && (tracker.registrator.get() != NULL)) {
      if (tracker.registrator->currTrackRegistered != tracker.registrator->subsampledTracks.end()) {
        PointCloudTrack::SPtr track = (**(tracker.registrator->currTrackRegistered)).track;
        vis.render3DTrackVariance(track.get());
      }
      // also render next track as preview
      if (tracker.registrator->nextTrackToRegister != tracker.registrator->subsampledTracks.end()) {
        PointCloudTrack::SPtr track = (**(tracker.registrator->nextTrackToRegister)).track;
        vis.render3DTrackVariance(track.get());
      }
    } else {
      callForEachTrack(frame, boost::bind(&AlgoVisualization::render3DTrackVariance, &vis, _1));
    }
  }
  if (ui.cbVisuTrackProps3D->isChecked()) {                                 // render Tracks
    callForEachTrack(frame, boost::bind(&AlgoVisualization::render3DTrackProps, &vis, _1));
  }

  glPopMatrix(); // Point-shift-up

  glEndList();
  emit stateChanged();
}

void VisualizerObstacleTracking::callForEachTrack(FrameSPtr frame, TrackProcessFunction func, bool skipWorldTrack)
{
  unsigned int sttIdx = ui.sbTrackBufferIdx->value();
  PointCloudTrack* wt = skipWorldTrack ? frame->getLttWorldTrack().get() : NULL;
  if (sttIdx == frame->sttValidationSteps) // render long-term-tracks
    callForEachTrack(frame->longTermTracks.begin(), frame->longTermTracks.end(), func, wt);
  if (sttIdx < frame->sttValidationSteps) // render short-term-tracks
    callForEachTrack(frame->shortTermTracks[sttIdx].tracks.begin(), frame->shortTermTracks[sttIdx].tracks.end(), func, wt);
};

void VisualizerObstacleTracking::renderTrackProjection(const LidarImage<LidarFrame::TrackIndex> &trackProj, FrameSPtr frame, unsigned int projSttIdx, QImage &img)
{
  if (projSttIdx <= frame->sttValidationSteps) {
    if (projSttIdx == frame->sttValidationSteps) {
      LidarFrame::TrackSPtrVec tmpTrkVec(frame->longTermTracks.begin(), frame->longTermTracks.end());
      vis.renderTrackProjection(trackProj, tmpTrkVec, img);
    } else {
      vis.renderTrackProjection(trackProj, frame->shortTermTracks[projSttIdx].tracks, img);
    }
  }
}

void VisualizerObstacleTracking::center3DCamPos()
{
  if (!glWid) return;
//  int currRow = getSelectedFrameNumber();
//  if (currRow < 0) return; // if no row is selected
//  FrameSPtr frame = tracker.getFrame(currRow);
//  // set camera to that position
  double pan; double tilt; double range; double x_offset; double y_offset; double z_offset;
  glWid->getCameraPos(pan, tilt, range, x_offset, y_offset, z_offset);
//  x_offset = frame->positionHTM2w(0,3);
//  y_offset = frame->positionHTM2w(1,3);
//  z_offset = frame->positionHTM2w(2,3);
//  pan = 0.0;
//  tilt = 90.0;
//  range = 25.0;
//  glWid->setCameraPos(pan, tilt, range, x_offset, y_offset, z_offset);
  glWid->setCameraPos(0.0, 90.0, 25.0, 0.0, 0.0, 0.0);
  render3D();
}
  
void VisualizerObstacleTracking::incPointSize()
{
  ParameterHeap *params = ParameterHeap::get();
  params->vis3DPointSize++;
  render3D();
}

void VisualizerObstacleTracking::decPointSize()
{
  ParameterHeap *params = ParameterHeap::get();
  if (params->vis3DPointSize > 1) params->vis3DPointSize--;
  render3D();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////         Data Processing           ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void VisualizerObstacleTracking::setParamsDefault()
{
  ParameterHeap *params = ParameterHeap::get();
  this->setUpdatesEnabled(false);
  blockParamHeapWrite = true;
  ui.dsbConnMaxDist->setValue(params->connMaxDst);
  ui.dsbConnMaxDistRatio->setValue(params->connMaxDstRatio);
  ui.dsbConvexThresh->setValue(params->segConvThresh);
  ui.dsbNormalThresh->setValue(params->segNormThresh);
  ui.dsbMinSegScore->setValue(params->segMinScore);
  ui.sbMinSegSizePt->setValue(params->segMinSizePt);
  ui.sbMinSegSizeCm->setValue(params->segMinSizeM*100);
  ui.cbSegUseTrackIDs->setChecked(params->segUseTrackIDs);
  ui.cbICPDistDiffWeight->setChecked(params->regIcpUseDistDiffWeight);
  ui.cbICPOcclWeight->setChecked(params->regIcpUseOcclWeight);
  ui.cbRegWithFeatMatch->setChecked(params->regUseFeatMatch);
  ui.cbDynamicMapping->setChecked(params->useDynamicMapping);
  ui.cbSimpleMapping->setChecked(params->useSimpleMapping);
  ui.cbOverwriteMapping->setChecked(params->useOverwriteMapping);
  ui.spTrackHist->setValue(params->visTrackHistory);
  ui.dsbHistPos->setValue(params->visMaxTrackVelVarPos);
  ui.dsbHistAng->setValue(params->visMaxTrackVelVarAng);
  blockParamHeapWrite = false;
  this->setUpdatesEnabled(true);
  this->update();
}
void VisualizerObstacleTracking::paramsChanged()
{
  if (blockParamHeapWrite) return;
  ParameterHeap *params = ParameterHeap::get();
  params->connMaxDst = ui.dsbConnMaxDist->value();
  params->connMaxDstRatio = ui.dsbConnMaxDistRatio->value();
  params->segConvThresh = ui.dsbConvexThresh->value();
  params->segNormThresh = ui.dsbNormalThresh->value();
  params->segMinScore = ui.dsbMinSegScore->value();
  params->segMinSizePt = ui.sbMinSegSizePt->value();
  params->segMinSizeM = (double)ui.sbMinSegSizeCm->value()/100.;
  params->segUseTrackIDs = ui.cbSegUseTrackIDs->isChecked();
  params->regIcpUseDistDiffWeight = ui.cbICPDistDiffWeight->isChecked();
  params->regIcpUseOcclWeight = ui.cbICPOcclWeight->isChecked();
  params->regUseFeatMatch = ui.cbRegWithFeatMatch->isChecked();
  params->useDynamicMapping = ui.cbDynamicMapping->isChecked();
  params->useSimpleMapping = ui.cbSimpleMapping->isChecked();
  params->useOverwriteMapping = ui.cbOverwriteMapping->isChecked();
  params->visTrackHistory = ui.spTrackHist->value();
  params->visMaxTrackVelVarPos = ui.dsbHistPos->value();
  params->visMaxTrackVelVarAng = ui.dsbHistAng->value();
}
void VisualizerObstacleTracking::setPVec1()
{
  bool pressedOK;
  QString pstring = QInputDialog::getText(this, "Parameter Vector Input", "Vector:", QLineEdit::Normal, QString(""), &pressedOK );
  if (pressedOK) {
    QStringList pslist = pstring.split(",", QString::SkipEmptyParts);
    cout << endl;
    if (pslist.size() == 1) // comma was not the separator -> split by spaces
      pslist = pstring.split(" ", QString::SkipEmptyParts);
    vector<double> pvec;
    for (int i=0; i<pslist.size(); ++i)
      pvec.push_back(pslist.at(i).trimmed().toFloat());
    ParameterHeap *params = ParameterHeap::get();
    params->setSelection1(pvec);
    setParamsDefault(); // update GUI elements
  }
}

void VisualizerObstacleTracking::newFrameAvailable()
{
  FrameSPtr f = tracker.getFrame(0);
  cerr << endl << "Grabbed frame " << f->sourcefile << endl << "  at time " << f->recTime << " = " << f->getTime() << flush;
  setMergingActive(false);
  if (ui.cbAutoCalculate->isChecked()) {
    processFrame();
  } else {
    render2D();
    render3D();
  }
}
  
void VisualizerObstacleTracking::resetBuffer()
{
  FrameSPtr f = tracker.getFrame(0);
  unsigned int buffSize = ui.sbBuffSize->value();
  tracker.reset(buffSize);
  setMergingActive(false);
  render2D();
  render3D();
}


void VisualizerObstacleTracking::setProcessingSelection(bool features, bool registration, bool segmentation, bool tracking)
{
  ui.gbFeatures->setChecked(features);
  ui.gbRegistration->setChecked(registration);
  ui.gbSegmentation->setChecked(segmentation);
  ui.gbTrackCreation->setChecked(tracking);
}

void VisualizerObstacleTracking::processFrame()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  //tracker.processFrame(currRow);
  if (ui.gbFeatures->isChecked())      tracker.calcFeatures(currRow);
  if (ui.gbRegistration->isChecked())  tracker.calcRegistration(currRow);
////  tracker.calcEgoMotion(currRow);
////  tracker.calcTrackEgoMotionUpdate(currRow);
  if (ui.gbSegmentation->isChecked())  tracker.calcSegmentation(currRow);
  if (ui.gbTrackCreation->isChecked()) tracker.calcTracks(currRow);
  render2D();
  render3D();
  if (ui.cbAutoStoreImg->isChecked())
    writeDebugImages();
}

void VisualizerObstacleTracking::processAllFrames()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) currRow = tracker.bufferSize()-1;
  while (currRow >= 0) {
    setSelectedFrameNumber(currRow);
    processFrame(); // will update GUI
    --currRow;
  }
}

void VisualizerObstacleTracking::generateFeatures()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  tracker.calcFeatures(currRow);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::generateSegments()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  tracker.calcSegmentation(currRow);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::generateSegmentColors()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  recolorSegments(tracker.getFrame(currRow));
  render2D();
  render3D();
}

void VisualizerObstacleTracking::generateAllSegmentColors()
{
  for (unsigned int row=0; row < tracker.bufferSize(); ++row) {
    try {
      recolorSegments(tracker.getFrame((int)row));
    } catch (range_error) {
    }
  }
  render2D();
  render3D();
}

void VisualizerObstacleTracking::generateRegistration()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  tracker.calcRegistration(currRow);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::setMergingActive(bool active) {
  ui.pbTMMerge->setEnabled(active);
  ui.sbTMMerge->setEnabled(active);
  ui.pbTMKeep->setEnabled(active);
  ui.pbTMIgnore->setEnabled(active);
  ui.pbTMBack->setEnabled(active);
  ui.pbTMPlus->setEnabled(active);
}

void VisualizerObstacleTracking::nextMergeTrack(FrameSPtr frame, bool jumpOver) // helper function called by the ones below
{
  (void)jumpOver;
  LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
  if (trackMergeIdx >= (int)buffer.tracks.size()-1) {
    trackMergeIdx = -1;
    setMergingActive(false);
    return;
  }
  PointCloudTrack::SPtr track;
  LidarFrame::TrackLinks links;
//  do {
    ++trackMergeIdx;
    track = buffer.tracks[trackMergeIdx];
    links = frame->getMergedSTTLinks(frame->sttValidationSteps-1, trackMergeIdx, true);
//  } while (jumpOver && (trackMergeIdx < (int)buffer.tracks.size()-1));// && (track->getLastUpdateCounter() > 0) && (links.size() <= 1));
  if (trackMergeIdx < (int)buffer.tracks.size()) {
    ui.sbTMMerge->setMinimum(0);
    ui.sbTMMerge->setMaximum(links.size()-1);
    ui.sbTMMerge->setValue(0);
    switch (track->mergeDecision) {
    case PointCloudTrack::MDIgnore:
      ui.pbTMIgnore->setFocus();
      break;
    case PointCloudTrack::MDKeep:
      ui.pbTMKeep->setFocus();
      break;
    case PointCloudTrack::MDMerge:
      ui.pbTMMerge->setFocus();
      break;
    }
    PointCloudTrack *linkedTrack = track->trackMaxMergeScore.get();
    if (linkedTrack) {
      for (unsigned int i=0; i<links.size(); ++i) {
        if (linkedTrack == links[i].track.get()) {
          ui.sbTMMerge->setValue(i);
          break;
        }
      }
    }
  } else {
    trackMergeIdx = -1;
    setMergingActive(false);
  }
}

void VisualizerObstacleTracking::trackMergeBack() // helper function called by the ones below
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) { // if no row is selected
    setMergingActive(false);
    return;
  }
  if (trackMergeIdx <= 5) { // will start from the beginning
    trackMergeInit();
    return;
  }
  FrameSPtr frame = tracker.getFrame(currRow);
  trackMergeIdx -= 6;
  nextMergeTrack(frame, false);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergePlus() // helper function called by the ones below
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) { // if no row is selected
    setMergingActive(false);
    return;
  }
  FrameSPtr frame = tracker.getFrame(currRow);
  trackMergeIdx += 20;
  nextMergeTrack(frame); // will set trackMergeIdx to -1 when finished
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergeInit()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) { // if no row is selected
    setMergingActive(false);
    return;
  }
  FrameSPtr frame = tracker.getFrame(currRow);
  cout << endl << "initializing track-management...current buffer:" << flush;
  for (unsigned int idx = 0; idx < frame->shortTermTracks.size(); ++idx) {
    if (idx == frame->sttValidationSteps-1)
      cout << "[" << frame->shortTermTracks[idx].tracks.size() << "]" << flush;
    else
      cout << "(" << frame->shortTermTracks[idx].tracks.size() << ")" << flush;
  }
  if (   (frame->shortTermTracks[frame->sttValidationSteps-1].tracks.size() > 0)
      && (frame->shortTermTracks[frame->sttValidationSteps].tracks.size() > 0)) {
    trackMergeIdx = -1; // next command will increase Idx by 1
    nextMergeTrack(frame, false);
    setMergingActive(true);
  } else {
    setMergingActive(false);
  }
  ui.cbVisuTrackMergeDec3D->setChecked(false);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergeMerge()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  FrameSPtr frame = tracker.getFrame(currRow);
  LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
  if ((trackMergeIdx >= 0) && (trackMergeIdx < (int)buffer.tracks.size())) {
    LidarFrame::TrackLinks links = frame->getMergedSTTLinks(frame->sttValidationSteps-1, trackMergeIdx, true);
    if ((ui.sbTMMerge->value() >= 0) && (ui.sbTMMerge->value() < (int)links.size())) {
      buffer.tracks[trackMergeIdx]->mergeDecision = PointCloudTrack::MDMerge;
      buffer.tracks[trackMergeIdx]->trackMaxMergeScore = links[ui.sbTMMerge->value()].track;
    }
    nextMergeTrack(frame);
  }
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergeKeep()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  FrameSPtr frame = tracker.getFrame(currRow);
  LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
  if ((trackMergeIdx >= 0) && (trackMergeIdx < (int)buffer.tracks.size())) {
    buffer.tracks[trackMergeIdx]->mergeDecision = PointCloudTrack::MDKeep;
    buffer.tracks[trackMergeIdx]->trackMaxMergeScore.reset(); // no merge-link
    nextMergeTrack(frame);
  }
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergeIgnore()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  FrameSPtr frame = tracker.getFrame(currRow);
  LidarFrame::ExtTrackList &buffer = frame->shortTermTracks[frame->sttValidationSteps-1];
  if ((trackMergeIdx >= 0) && (trackMergeIdx < (int)buffer.tracks.size())) {
    LidarFrame::TrackLinks links = frame->getMergedSTTLinks(frame->sttValidationSteps-1, trackMergeIdx, true);
    if ((ui.sbTMMerge->value() >= 0) && (ui.sbTMMerge->value() < (int)links.size())) {
      buffer.tracks[trackMergeIdx]->mergeDecision = PointCloudTrack::MDIgnore;
      buffer.tracks[trackMergeIdx]->trackMaxMergeScore = links[ui.sbTMMerge->value()].track;
    }
    nextMergeTrack(frame);
  }
  render2D();
  render3D();
}

void VisualizerObstacleTracking::trackMergeSBInc()
{
  ui.sbTMMerge->setValue(ui.sbTMMerge->value()+1);
}

void VisualizerObstacleTracking::trackMergeSBDec()
{
  ui.sbTMMerge->setValue(ui.sbTMMerge->value()-1);
}

void VisualizerObstacleTracking::generateTracks()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  trackMergeIdx = -1;
  tracker.calcTracks(currRow); // will switch buffers -> index is not valid any more
  setMergingActive(false);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::resetSegReg()
{
  int currRow = getSelectedFrameNumber();
  if (currRow < 0) return; // if no row is selected
  tracker.registrationInit(currRow);
  render2D();
  render3D();
}

void VisualizerObstacleTracking::registerCurrSeg()
{
  tracker.registrationNext();
  render2D();
  render3D();
}

void VisualizerObstacleTracking::registerFinish()
{
  while (tracker.registrationNext()) {};
  render2D();
  render3D();
}

void VisualizerObstacleTracking::registerFinishLayer()
{
  while (tracker.registrationNext(true)) {};
  render2D();
  render3D();
}

