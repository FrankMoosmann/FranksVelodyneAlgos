#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <EnumCreation.hpp>
#include <ParamStructCreation.hpp>
#include <LidarImageProjectorPNG.hpp>
#include <Gui3DQt/Gui.hpp>
#include <Gui3DQt/VisualizerCamControl.hpp>
#include <Timer.hpp>
#include <FMUtils.hpp>

#include "ObstacleTracking.hpp"
#include "VisualizerPNGReader.hpp"
#include <NormalSampler.hpp>
#include <ParameterHeap.hpp>

#include <boost/static_assert.hpp>
#include <kogmo_time.h>
#include <MatrixDefs.hpp>
#include <HomogeneousTransformationMatrix.hpp>


// definition of Main Parameters
#define MAIN_PARAMS_ARR (7, ( \
       (std::string,               sourcedir,        ".", \
                                                               "specify directory with png files", 0) \
      ,(unsigned int,              from,             0,        "processing starts with this image number", 0) \
      ,(unsigned int,              to,               0,        "if specified, processing immediately runs without GUI until this image number", 0) \
      ,(unsigned int,              horizRes,         870,      "horizontal resolution of range image (for UDP|RTDB)", 0) \
      ,(float,                     imgStartAng,      180,      "start angle for image acquisition (for UDP|RTDB)", 0) \
      ,(unsigned int,              bufferSize,       3,        "size of tracking buffer", 0) \
      ,(bool,                      profile,          false,    "if specified, profiling is carried out on the image range specified", 0) \
      ) )
CREATE_PARAM_STRUCT(MainParams, "Main-Parameter", MAIN_PARAMS_ARR)
CREATE_PARAM_STRUCT_STREAMOP(MainParams, MAIN_PARAMS_ARR, " ", 0)


void handler(int sig) {
  void *array[30];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 30);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

int main(int argc, char *argv[])
{
  using namespace std;
  namespace fs = boost::filesystem;
  namespace po = boost::program_options;
  
  const fs::path currPath = fs::initial_path();
  const fs::path binPath = fs::path(argv[0]).branch_path().relative_path();
  
  /*
  signal(SIGSEGV, handler);   // install our handler
  signal(SIGINT, handler);
  signal(SIGTERM, handler);
  signal(SIGABRT, handler);
  signal(SIGFPE, handler);
  signal(SIGILL, handler);
  */

  // ------------------------------------------------
  // -------------- parse arguments -----------------
  // ------------------------------------------------
  ParameterHeap* params = ParameterHeap::get();
  MainParams mainParams;
  mainParams.parameterDescription.add_options()
    ("help", "produce standard help message")
    ("full-help", "produce extended help message")
    ;
  po::options_description allOptions;
  allOptions.add(mainParams.parameterDescription);
  allOptions.add(params->fullParameterDescription);
  po::positional_options_description posOpt;
  posOpt.add("sourcedir",1);
  posOpt.add("from",1);
  posOpt.add("to",1);
  po::variables_map poVM;
  try {
    //po::store(parse_command_line(argc, argv, mainParams, po::command_line_style::unix_style ^ po::command_line_style::allow_short), poVM);
    po::store(po::command_line_parser(argc, argv).options(allOptions).positional(posOpt).run(), poVM);
    po::notify(poVM);
  } catch (exception &e) {
    cerr << "Error when parsing arguments:" << e.what() << endl;
    return 1;
  }

  if (poVM.count("help")) {
    cout << mainParams.parameterDescription << endl;
    return 0;
  }

  if (poVM.count("full-help")) {
    cout << allOptions << endl;
    return 0;
  }

  // ------------------------------------------------
  // -------------- run program ---------------------
  // ------------------------------------------------
  const boost::array<double,25> pvec = {{
/*-286520+*/  0.202108, 1.81887, 0.520373, -1.05362, 0.103174, -2.41372, 2.09533, 20.3268, 4.85689, 15.4985, 0.522013, 8, 0.456825, 1, 0, 0, 1, 1, 1, 1, 0, 4, 0, 0, 1
/*-209471*/  //0.147802, 1.90782, 0.25095, -2.04829, 0.141517, -7.67024, 1.16484, 19.1772, 5.16697, 14.8062, 0.45036, 8         , 0.450437, 1, 0, 0, 1, 1, 1, 1, 0, 4, 0, 0, 1
/*-235143*/  //0.20203, 1.81866, 0.963305, 0.769714, -0.00338368, -6.43715, 1.93201, 19.8571, 5.01758, 16.9035, 0.387119, 8   , 0.317836, 1, 0, 0, 1, 1, 1, 1, 2, 4, 0, 0, 1
/*-291012*/  //0.134167, 1.72293, 2.23169, 2.02115, 0.0044602, -14.1693, 0.977674, 19.844, 4.85896, 14.4566, 0.376234, 8      , 0.384149, 1, 1, 0, 1, 1, 1, 1, 1, 4, 0, 0, 1
  }};
  params->setSelection1(pvec);

  LidarImageProjector *proj = NULL;
  PNGReader *pngReader = NULL;
  ObstacleTracking *tracker = NULL;
  VisualizerObstacleTracking *trackerVis = NULL;
  Gui3DQt::Visualizer *readerVis = NULL;
  // read png images from disk
  pngReader = new PNGReader(mainParams.sourcedir);
  if (mainParams.to > pngReader->getImageCount()) {
    cout << endl << "reducing target image number to " << pngReader->getImageCount() << flush;
    mainParams.to = pngReader->getImageCount();
  }
  if (mainParams.from > pngReader->getImageCount()) {
    cout << endl << "reducing start image number to " << pngReader->getImageCount() << flush;
    mainParams.from = pngReader->getImageCount();
  }
  proj = new PNGImageProjector(pngReader->getImageConfigName());
  tracker = new ObstacleTracking(mainParams.bufferSize, *proj);
  if (mainParams.to > mainParams.from) {
    cout << endl << "processing from " << mainParams.from << " to " << mainParams.to << flush;
    // simply process range, then exit
    if (mainParams.profile) {
      //bool features, bool registration, bool segmentation, bool tracking)
      trackerVis->setProcessingSelection(true, true, true, true);
//        trackerVis->setProcessingSelection(true, false, true, false);
      cout << endl << "starting to profile" << flush;
    }
    Timer t;
    for (unsigned int i=mainParams.from; i<mainParams.to; ++i) {
      tracker->next(); // shift buffers
      FrameSPtr cf = tracker->getFrame(0); // get current buffer
      FrameSPtr lf = tracker->getFrame(1);
      t.start("read");
      pngReader->readImage(i, cf, lf);
      t.stop();
      t.start("process");
      tracker->processFrame();
      t.stop();
    }
    if (mainParams.profile) {
      cout << endl;
      t.plot();
    }
    delete tracker; // will write the remaining points to disk
    delete proj;
    delete pngReader;
    cout << endl;
    return 0;
  }
  Gui3DQt::Gui myGui("DATMOSLAM", argc, argv, Gui3DQt::MainWindow::GM_3D2D);
  Gui3DQt::MNavWidget* glWid = myGui.getQGlWidget();
  trackerVis = new VisualizerObstacleTracking(*tracker, 0, glWid);
  readerVis = new VisualizerPNGReader(*pngReader, *tracker, *trackerVis, false);
  ((VisualizerPNGReader*)readerVis)->selectImg(mainParams.from);
  // start GUI
  myGui.registerVisualizer(readerVis, "Reader", Gui3DQt::MainWindow::VM_Plain);
  myGui.registerVisualizer(trackerVis, "Tracker", Gui3DQt::MainWindow::VM_Plain);
  myGui.registerVisualizer(new Gui3DQt::VisualizerCamControl(*myGui.getQGlWidget()), "Pilot");
  cout << endl << endl << "bringing up GUI..." << flush; sleep(1);
  myGui.exec();
  // don't delete readerVis and trackerVis, they are owned by GUI
  delete tracker; // will write the remaining points to disk
  delete proj;
  delete pngReader;
  cout << endl;
  return 0;
}
