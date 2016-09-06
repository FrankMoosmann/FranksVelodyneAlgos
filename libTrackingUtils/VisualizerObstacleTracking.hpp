#ifndef VISUALIZEROBSTACLETRACKING_H
#define VISUALIZEROBSTACLETRACKING_H

#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <GL/glut.h>
#include <QtCore/QThread>

#include <Gui3DQt/Visualizer.hpp>
#include <Gui3DQt/MNavWidget.hpp>

#include "ui_VisualizerObstacleTracking.h"
#include "LidarFrame.hpp"
#include "ObstacleTracking.hpp"
#include "AlgoVisualization.hpp"


class VisualizerObstacleTracking : public Gui3DQt::Visualizer
{
  Q_OBJECT

public:
  VisualizerObstacleTracking(ObstacleTracking &tracker, int dummy, Gui3DQt::MNavWidget* glWid, QWidget *parent = 0);
  virtual ~VisualizerObstacleTracking();

  virtual void paintGLOpaque();
  virtual void paintGLTranslucent();

  void setProcessingSelection(bool features, bool registration, bool segmentation, bool tracking);

public slots:
  void newFrameAvailable();
  void processFrame();

public: // should be private??
  Ui::VisualizerObstacleTrackingClass ui;
  ObstacleTracking &tracker;

  Gui3DQt::MNavWidget* glWid;
  AlgoVisualization vis;

  GLuint glListIndex;
  
  int trackMergeIdx;
  void setMergingActive(bool);
  void nextMergeTrack(FrameSPtr frame, bool jumpOver = true);
  bool blockParamHeapWrite;

  int getSelectedFrameNumber() const;
  void setSelectedFrameNumber(int);

  typedef boost::function<void (PointCloudTrack *track)> TrackProcessFunction;
  void callForEachTrack(FrameSPtr frame, TrackProcessFunction func, bool skipWorldTrack = false); // calls function of each selected track to visualize
  template <class TrackSPtrInputIterator>
  void callForEachTrack(TrackSPtrInputIterator begin, TrackSPtrInputIterator end, TrackProcessFunction func, PointCloudTrack *skipTrack = NULL)
  {
    for (; begin != end; ++begin) {
      if (begin->get() != skipTrack)
        func(begin->get());
    }
  };
  void renderTrackProjection(const LidarImage<LidarFrame::TrackIndex> &trackProj, FrameSPtr frame, unsigned int projSttIdx, QImage &img);

public slots:
  void render2D();
  void render3D();
  void center3DCamPos();
  void incPointSize();
  void decPointSize();

  void setPVec1();
  void paramsChanged();
  void setParamsDefault();
  void resetBuffer();
  void writeDebugImages();
  
  void processAllFrames();
  void generateFeatures();
  void generateSegments();
  void generateSegmentColors();
  void generateAllSegmentColors();
  void generateRegistration();
  void trackMergeBack();
  void trackMergePlus();
  void trackMergeInit();
  void trackMergeMerge();
  void trackMergeKeep();
  void trackMergeIgnore();
  void trackMergeSBInc();
  void trackMergeSBDec();
  void generateTracks();

  void resetSegReg();
  void registerCurrSeg();
  void registerFinishLayer();
  void registerFinish();

  void saveRegistrConfig();
  void saveSegments();
  void loadSegments();

};

#endif // VISUALIZEROBSTACLETRACKING_H
