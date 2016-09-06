#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// includes from USR_INC:
#include <LidarImageProjectorPNG.hpp>
#include <Gui3DQt/Gui.hpp>
#include <Gui3DQt/VisualizerGrid.hpp>
#include <Gui3DQt/VisualizerPassat.hpp>
#include <Gui3DQt/VisualizerCamControl.hpp>

// includes from this project:
#include "VisualizerTracks.hpp"

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // ------------------------------------------------
  // ------------ parse arguments -------------------
  // ------------------------------------------------
  string trackdir = ".";
  string outdir = ".";
  string pilotSource = "";
  float  shiftup = 0.0;
  unsigned int bufferMemMB = 100;
  po::options_description progOpt;
  progOpt.add_options()
    ("help", "produce help message")
    ("noPilot", "specify to suppress creation of pilot")
    ("visNoBB", "specify to suppress bounding-box-visualization")
    ("trackdir", po::value(&trackdir), "Directory containing track details")
    ("bufferMemMB", po::value(&bufferMemMB)->default_value(bufferMemMB), "specify the amount of memory used for buffering neighboring frames")
    ("skipPoints", "specify in order to not load track appearance")
    ("loadWorld", "specify in order to load accumulated 3D world model")
    ("pilotBuffer", po::value(&pilotSource), "load pilot buffer from specified file")
    ("shiftup", po::value(&shiftup)->default_value(shiftup), "shift 3D scene upwards by this amount of meter")
    ("autogenScreenshots", po::value(&outdir), "automatically generated output images in specified directory")
    ;

  po::positional_options_description posOpt;
  posOpt.add("trackdir", 1);

  po::variables_map poVM;
  try {
    po::store(po::command_line_parser(argc, argv).options(progOpt).positional(posOpt).run(), poVM);
    po::notify(poVM);
  } catch (exception e) {
    cerr << "Error when parsing arguments" << endl << progOpt << endl;
    return 1;
  }

  if (poVM.count("help")) {
    cout << progOpt << endl;
    return 1;
  }

  if (!fs::exists(trackdir)) {
    cout << "Image directory \"" << trackdir << "\" is invalid" << endl;
    cout << progOpt << endl;
    return 1;
  }

  bool loadPoints = !poVM.count("skipPoints");
  bool loadWorld = poVM.count("loadWorld");
  bool autogen = (poVM.count("autogenScreenshots")) && (fs::exists(outdir));
  VisualizerTracks::VisuMode visMode = poVM.count("visNoBB") ? VisualizerTracks::VM_TPC : VisualizerTracks::VM_TBB;

  // ------------------------------------------------
  // -------------- run program ---------------------
  // ------------------------------------------------
  Gui3DQt::Gui myGui("TrackVisualizer", argc, argv, Gui3DQt::MainWindow::GM_3D, true);
  Gui3DQt::VisualizerGrid *gVis = new Gui3DQt::VisualizerGrid();
  Gui3DQt::VisualizerPassat *pVis = new Gui3DQt::VisualizerPassat(0.5, 0.5, 0.5, 0, false); // double r = 0.5, double g = 0.5, double b = 0.5, unsigned int modelindex = 0, bool rotateVelodyne = true
  pVis->setPose(Gui3DQt::VisualizerPassat::translate_to_velobase_x,
                Gui3DQt::VisualizerPassat::translate_to_velobase_y,
                Gui3DQt::VisualizerPassat::translate_to_groundcenter_z,
                0);
  VisualizerTracks *tVis = new VisualizerTracks(trackdir,bufferMemMB,loadPoints,loadWorld,shiftup,visMode,autogen,pVis);
  Gui3DQt::VisualizerCamControl *cVis = new Gui3DQt::VisualizerCamControl(*myGui.getQGlWidget());
  myGui.registerVisualizer(tVis, "Tracks", Gui3DQt::MainWindow::VM_Plain);
  myGui.registerVisualizer(pVis, "Passat");
  myGui.registerVisualizer(gVis, "Grid", Gui3DQt::MainWindow::VM_Groupbox, false);
  if (!poVM.count("noPilot"))
    myGui.registerVisualizer(cVis, "Pilot");
  if ((poVM.count("pilotBuffer")) && (fs::exists(pilotSource)))
    cVis->loadCamBuff(pilotSource);
  if (autogen) {
    if (poVM.count("pilotBuffer")) cVis->camPosSelect();
    myGui.getMainWindow()->setControlPanelVisible(false);
    myGui.getMainWindow()->setImageOutputDir(outdir);
    myGui.getMainWindow()->setGrabbingActive(true);
    tVis->imgPlay(true);
  }
  cout << endl << endl << "bringing up GUI..." << flush; sleep(1);
  myGui.exec();

  cout << endl;
  return 0;
}
