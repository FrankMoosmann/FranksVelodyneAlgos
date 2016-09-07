#include "VisualizerPNGReader.hpp"


#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <kogmo_time.h>
#include <FMUtils.hpp>

#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////        VisualizerPNGReader        ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

VisualizerPNGReader::VisualizerPNGReader(PNGReader &pngReader_, ObstacleTracking &tracker_, VisualizerObstacleTracking &trackerVis_, bool fillBuffer, QWidget *parent)
  : Gui3DQt::Visualizer(parent)
  , tracker(tracker_)
  , trackerVis(trackerVis_)
  , pngReader(pngReader_)
{
  continReader = new QTimer(this);

  ui.setupUi(this);
  ui.hsImageNb->setMinimum(0);
  ui.hsImageNb->setMaximum(pngReader.getImageCount()-1);
  ui.hsImageNb->setValue(0);
  //imgSelect(ui.hsImageNb->value());
  connect(ui.hsImageNb, SIGNAL(valueChanged(int)), this, SLOT(imgSelect(int)) );
  connect(ui.pbReadImage, SIGNAL(pressed()), this, SLOT(readImage()) );
  connect(ui.pbReadN, SIGNAL(pressed()), this, SLOT(readNImages()) );
  connect(ui.pbEmpty, SIGNAL(pressed()), this, SLOT(fillBuffer()) );
  connect(ui.pbReadContinuous, SIGNAL(toggled(bool)), this, SLOT(readContinuously(bool)) );
  connect(continReader, SIGNAL(timeout()), this, SLOT(readImage()) );
  connect(ui.pbLoop, SIGNAL(pressed()), this, SLOT(loopToNextTrackMerge()) );
  connect(ui.pbLucky, SIGNAL(pressed()), this, SLOT(iFeelLucky()) );
  if (fillBuffer)
    this->fillBuffer();
}


VisualizerPNGReader::~VisualizerPNGReader()
{
}

void VisualizerPNGReader::imgSelect(int val)
{
  path fullfile(pngReader.getDstImageName(val));
  ui.lImageName->setText(QString(fullfile.filename().c_str()));
}

void VisualizerPNGReader::readImage() // this method does the job
{
  // 1) determine current index
  unsigned int fileIdx = ui.hsImageNb->value();

  // 2) load png-file and transform to frame
  bool incBufferIdx = true;
  FrameSPtr lf = tracker.getFrame(incBufferIdx ?  0 : 1);
  FrameSPtr cf = tracker.getFrame(incBufferIdx ? -1 : 0); // get future buffer pointer
  pngReader.readImage(fileIdx, cf, lf);
  if (incBufferIdx) tracker.next();

  // 3) increment slider-position, if already last position, stop timer
  if (ui.hsImageNb->value() >= ui.hsImageNb->maximum()) {
    continReader->stop();
  } else {
    ui.hsImageNb->setValue(fileIdx+1);
  }

  // 4) signal obstacleTracker that new frame is available
  trackerVis.newFrameAvailable(); // update GUI and possibly process frame
}

void VisualizerPNGReader::readNImages() // this method does the job
{
  for (int n=0; (n < ui.sbReadN->value()) && (ui.hsImageNb->value() < ui.hsImageNb->maximum()); ++n)
    readImage();
}

void VisualizerPNGReader::fillBuffer()
{
  size_t limit = min(tracker.bufferSize(), (unsigned int)pngReader.getImageCount());
  for (size_t i=0; i<limit; ++i) {
    readImage();
  }
}

void VisualizerPNGReader::readContinuously(bool pressed)
{
  if (pressed)
    continReader->start(100); // interval in ms
  else
    continReader->stop();
}


void VisualizerPNGReader::loopToNextTrackMerge()
{
  for (int n=0; (n < ui.sbLoopCount->value()) && (ui.hsImageNb->value() < ui.hsImageNb->maximum()); ++n) {
    trackerVis.generateTracks(); // will skip automatically if tracks were already generated
    trackerVis.ui.cbAutoCalculate->setChecked(false);
    readImage();
    trackerVis.generateFeatures();
    trackerVis.generateRegistration();
    trackerVis.generateSegments();
    trackerVis.ui.cbVisuTrackMergeDec3D->setChecked(true);
  }
  // beep ?
}

void VisualizerPNGReader::iFeelLucky()
{
  // check if CreateTracks was already pressed and if not execute this step
  trackerVis.ui.cbAutoCalculate->setChecked(true);
  ui.sbReadN->setValue(4);
  readNImages();
  loopToNextTrackMerge();
  // beep ?
}
