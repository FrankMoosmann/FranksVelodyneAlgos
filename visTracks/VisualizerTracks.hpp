#ifndef VISUALIZERTRACKS_H_
#define VISUALIZERTRACKS_H_

#include <GL/glut.h>
#include <QTimer>
#include <QThread>
#include <QMutex>
#include <vector>

#include <Gui3DQt/Visualizer.hpp>
#include <Gui3DQt/PointCloudRenderer.hpp>
#include <Gui3DQt/VisualizerPassat.hpp>
#include <LidarImageProjectorPNG.hpp>
#include <MatrixDefs.hpp>

#include "ui_VisualizerTracks.h"
#include "FrameData.hpp"

class VisualizerTracks : public Gui3DQt::Visualizer
{
  Q_OBJECT

public:
  enum VisuMode {VM_TBB, VM_TPC}; // track-bounding-boxes, track-point-clouds

  VisualizerTracks(std::string trackDir, unsigned int bufferMemMB, bool loadPoints = true, bool loadWorld = false, float shiftup = 0.0, VisuMode mode = VM_TBB, bool exitAfterLastFrame = false, Gui3DQt::VisualizerPassat *passatVis = NULL, QWidget *parent = 0);
  virtual ~VisualizerTracks();

  virtual void paintGLOpaque();
  virtual void paintGLTranslucent();

private:
  class PrecacheThread : public QThread {
  private:
    VisualizerTracks &vis;
    unsigned long long bufferMem;
    unsigned long long precache(int idx); // returns size in bytes
  public:
    PrecacheThread(VisualizerTracks &vis_, unsigned int bufferMemMB) : QThread(&vis_), vis(vis_), bufferMem(bufferMemMB*1024*1024) {};
    virtual ~PrecacheThread() {};
    void run();
    void freeMem(int idx);
  };

  Gui3DQt::VisualizerPassat *passatVis;
  Ui::VisualizerTracksClass ui;
  bool listIsUpdating; // avoid repainting when items are added/removed
  QTimer *playTimer;
  PrecacheThread *precacheThread;

  // should all have the same size
  std::vector<std::string>     images;
  std::vector<std::string>     tracks;
  std::vector<FrameData*>      loadedTracks;
  QMutex                       trackMutex;
  const TrackData             *currWorldTrack; // will be NULL or will point to a valid TrackData object
  Gui3DQt::PointCloudRenderer *worldModelRenderer;
  float                        worldMinZ;
  float                        worldMaxZ;
  PNGImageProjector           *proj;
  bool                         loadPoints;
  bool                         exitAfterLastFrame;

  GLuint glListIndex;
  unsigned int ptSizeScan;
  unsigned int ptSizeMovingTrk;
  unsigned int ptSizeStaticTrk;
  unsigned int ptSizeWorld;

  static QString trackIDString(const TrackData&);
  static matrixTools::DCMatrix glHtm2TrackCS(const TrackData&);
  void loadWorldModelRenderer(std::string p3dFile); // safe loading from file
  void recalcWorldMinMaxZ();
  void drawMovingSurface(const TrackData::Surface &s);
  void drawStaticSurface(const TrackData::Surface &s);
  void drawSurface(const TrackData::Surface &s, bool height, bool normal, bool nConf);
  bool skipTrack(const TrackData &tdata);
  void applyColor(const TrackData &tdata);
  FrameData* loadFrame(int idx);

  inline float rotateZ(const float &x, const float &y, const float &z, const float &roll, const float &pitch) { // method for approximately calculating updated z coordinate after rotation
    return z - x*tan(pitch) + y*tan(roll);
  }

public slots:
  void imgPlay(bool);
  
private slots:
  void emitStateChanged() {emit stateChanged();};
  void recolorWorldModel(bool doUpdate) {if (doUpdate) recolorWorldModel();};
  void recolorWorldModel();
  void update3D(bool doUpdate) {if (doUpdate) update3D();};
  void update3D();

  void imgNext();
  void imgPrev();
  void imgSelect(int val);
  void updateSelectedImg(int val);
  void assertFrameData(int val);
  void selectAll(bool);
  void selectAll() {selectAll(true);};
  void selectNone() {selectAll(false);};
  void selectInvert();

  void incImgPtSize() {++ptSizeScan; update3D();};
  void decImgPtSize() {if (ptSizeScan >= 2) --ptSizeScan; update3D();};
  void incTrkPtSize() {++ptSizeMovingTrk; update3D();};
  void decTrkPtSize() {if (ptSizeMovingTrk >= 2) --ptSizeMovingTrk; update3D();};
  void incStatPtSize() {++ptSizeStaticTrk; update3D();};
  void decStatPtSize() {if (ptSizeStaticTrk >= 2) --ptSizeStaticTrk; update3D();};
  void incWrldPtSize() {++ptSizeWorld; recolorWorldModel();};
  void decWrldPtSize() {if (ptSizeWorld >= 2) --ptSizeWorld; recolorWorldModel();};
  void shiftupChanged();

};

#endif // VISUALIZERTRACKS_H_
