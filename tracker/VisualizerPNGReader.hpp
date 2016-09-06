#ifndef VISUALIZERPNGREADER_H_
#define VISUALIZERPNGREADER_H_

#include <GL/glut.h>
#include <QTimer>

#include <Gui3DQt/Visualizer.hpp>
#include <ObstacleTracking.hpp>
#include <VisualizerObstacleTracking.hpp>

#include "PNGReader.hpp"
#include "ui_VisualizerPNGReader.h"


class VisualizerPNGReader : public Gui3DQt::Visualizer
{
  Q_OBJECT

public:
  VisualizerPNGReader(PNGReader &pngReader_, ObstacleTracking &tracker, VisualizerObstacleTracking &trackerVis, bool fillBuffer = false, QWidget *parent = 0); //!< param dummy needed against segfault. dunny why. problem with qt
  virtual ~VisualizerPNGReader();

  virtual void paintGLOpaque() {};
  virtual void paintGLTranslucent() {};

  void selectImg(unsigned int nb) {ui.hsImageNb->setValue(nb);};

private:
  Ui::VisualizerPNGReaderClass ui;
  QTimer *continReader;

  ObstacleTracking &tracker;
  VisualizerObstacleTracking &trackerVis;
  PNGReader &pngReader;

//private slots:
public slots:
  void imgSelect(int val);
  void readImage();
  void readNImages();
  void fillBuffer();
  void readContinuously(bool pressed);
  void loopToNextTrackMerge();
  void iFeelLucky();

};

#endif // VISUALIZERPNGREADER_H_
