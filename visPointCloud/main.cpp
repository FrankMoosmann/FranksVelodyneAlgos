#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <Gui3DQt/Gui.hpp>
#include <Gui3DQt/VisualizerCamControl.hpp>
#include <Gui3DQt/VisualizerGrid.hpp>

#include "Visualizer3DMap.hpp"

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // ------------------------------------------------
  // ------------ parse arguments -------------------
  // ------------------------------------------------
  string mapfile;
  int from, to;
  po::options_description progOpt;
  progOpt.add_options()
    ("help", "produce help message")
    ("inputmap,i", po::value(&mapfile), "REQUIRED: Map-mainfile")
    ("from", po::value(&from)->default_value(-1), "Skip this number of points in the beginning")
    ("to", po::value(&to)->default_value(-1), "Abort at this number points index")
    ("traj,t", "set to only show trajectories for .txt maps")
    ;
  po::positional_options_description posOpt;
  posOpt.add("inputmap", 1);

  po::variables_map poVM;
  try {
    po::store(po::command_line_parser(argc, argv).options(progOpt).positional(posOpt).run(), poVM);
    po::notify(poVM);
  } catch (exception e) {
    cout << "Error when parsing arguments" << endl << progOpt << endl;
    return 1;
  }

  if (poVM.count("help")) {
    cout << "This is a small program to efficiently visualize large point clouds" << endl;
    cout << "Several file-types can be loaded:" << endl;
    cout << " .txt:     one point per line: pose-idx  x  y  z  x_orig  y_o  z_o  nx  ny  nz  nConf  scanD" << endl;
    cout << " .p3d:     text file with  x z y  all in one line" << endl;
    cout << " .p3di:    text file with  x z y intens  all in one line" << endl;
    cout << " .p3d(i)b: the two above but in binary format (x/y/z: double, i:uchar)" << endl;
    cout << "Program arguments:" << endl;
    cout << progOpt << endl;
    return 0;
  }
  if (!fs::exists(mapfile)) {
    cout << "mapfile \"" << mapfile << "\" does not exist" << endl;
    cout << progOpt << endl;
    return 1;
  }

  // ------------------------------------------------
  // -------------- run program ---------------------
  // ------------------------------------------------
  Gui3DQt::Gui myGui("MapVisualizer", argc, argv, Gui3DQt::MainWindow::GM_3D);
  myGui.registerVisualizer(new Visualizer3DMap(mapfile, 0, from, to, poVM.count("traj")), "3D Map", Gui3DQt::MainWindow::VM_Plain);
  myGui.registerVisualizer(new Gui3DQt::VisualizerCamControl(*myGui.getQGlWidget()), "Cam Control", Gui3DQt::MainWindow::VM_Plain);
  cout << endl << endl << "bringing up GUI..." << flush; sleep(1);
  myGui.exec();

  cout << endl;
  return 0;
}
