#ifndef VISUALIZER3DMAPPER_H_
#define VISUALIZER3DMAPPER_H_

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <GL/glut.h>
#include <QTimer>

#include <MatrixDefs.hpp>
#include <CsvReader.hpp>
#include <LidarImageProjector.hpp>
#include <AlgoVisualization.hpp>
#include <Gui3DQt/Visualizer.hpp>

#include "ui_Visualizer3DMapper.h"
#include "Frame.hpp"
#include "Map.hpp"
#include "Mapper.hpp"

class Visualizer3DMapper : public Gui3DQt::Visualizer
{
  Q_OBJECT

public:
  Visualizer3DMapper(Mapper &mapper, QWidget *parent = 0);
  virtual ~Visualizer3DMapper();

  virtual void paintGLOpaque();
  virtual void paintGLTranslucent();

public:
//private:
  Ui::Visualizer3DMapperClass ui;

//  const double maxDist; //!< maximum scanning distance to use points from
  Mapper &mapper;
  QTimer *continReader;
  QString lastMapPath;

  AlgoVisualization visAlg;
  GLuint glListIndex;
  unsigned int ptSize;

//private slots:
public slots:
  void update3D(bool doUpdate) {if (doUpdate) update3D();};
  void update3D();
  void incPtSize() {++ptSize; update3D();};
  void decPtSize() {if (ptSize > 2) --ptSize; update3D();};
  void imgSelect(int val);
  void processSelectedImage();
  void processCertainNumber();
  void setTimer(bool pressed);
  void stopTimer() {setTimer(false);};
  void resetMapping();
  void saveMap();
  void loadMap();
  void exportMap();
  void exportMap3D();
  void importMap();
  void importMap3D();
  void exportTraj();
  void importTraj();
  void readScan();
  void registerScan();
  void reg1();
  void reg2();
  void reg3();
  void filterMap();
  void filterMap_init();
  void filterMap_step();
  void filterMap_finish();
  void addScan();

};

#endif // VISUALIZER3DMAPPER_H_
