#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <Gui3DQt/Gui.hpp>
#include <Gui3DQt/VisualizerCamControl.hpp>
#include <GridND.hpp>

// includes from this project:
#include "DatasetReaderVeloslam.hpp"
#include "DatasetReaderVelotracking.hpp"
#include "Visualizer3DMapper.hpp"

using namespace std;
using namespace matrixTools;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
//  GridNDtest(false);  cout << endl;
  //testOrthonormalization(); //  cout << endl;

  // ------------------------------------------------
  // ------------ parse arguments -------------------
  // ------------------------------------------------
  string imgdir = ".";
  string savefile = "";
  string exportfile = "";
  vector<unsigned int> ranges;
  bool adept, deskew, filter, remOldMap, remFarMap, gps;
  unsigned int mapResolCM;
  double maxDist;
  const unsigned int nbRegSubsamples = 2000;
  po::options_description progOpt;
  progOpt.add_options()
    ("help", "produce help message")
    ("directory", po::value(&imgdir), "Directory containing range images")
    ("range",po::value<vector<unsigned int>>(&ranges)->multitoken(),"if specified all these image ranges {a,b} are processed straight away. Indexing starts at 0")
    ("deskew,d", po::value<bool>(&deskew)->default_value(true), "specify in order to deskew scans")
    ("adept,a", po::value<bool>(&adept)->default_value(true), "specify in order to use scan adaption")
    ("gps,g", "use GPS for localization instead of scan matching")
    ("incremental,i", "specify in order to achieve incremental scan matching only")
    ("local,l", "specify in order to keep only local map")
    ("filter,f", "specify in order to filter map afterwards")
    ("resolution,r", po::value<unsigned int>(&mapResolCM)->default_value(5), "map resolution in cm")
    ("maxScanDist,m",po::value(&maxDist)->default_value(50.0),"maximum range values used for mapping (in meter)")
    ("save", po::value<string>(&savefile), "if this and range is specified, map is saved to the specified file. disables GUI")
    ("export", po::value<string>(&exportfile), "if this and range is specified, map is exported to the specified file. disables GUI")
    ("gui", "start GUI although processed scans were already stored into a file")
    ;

  po::positional_options_description posOpt;
  posOpt.add("directory", 1);

  po::variables_map poVM;
  try {
    po::store(po::command_line_parser(argc, argv).options(progOpt).positional(posOpt).run(), poVM);
    po::notify(poVM);
  } catch (exception e) {
    cerr << "Error when parsing arguments" << endl << progOpt << endl;
    return 1;
  }
  gps = poVM.count("gps");
  remOldMap = poVM.count("incremental");
  remFarMap = poVM.count("local");
  filter = poVM.count("filter");

  if (poVM.count("help")) {
    cout << progOpt << endl;
    return 1;
  }
  if (!fs::exists(imgdir)) {
    cout << "Image directory \"" << imgdir << "\" is invalid" << endl;
    cout << progOpt << endl;
    return 1;
  }
  if (ranges.size() % 2 != 0) {
    cout << "Ranges must be specified as pairs. E.g. \"--range 1 10 12 15\"" << endl;
    cout << progOpt << endl;
    return 1;
  }

  bool showGUI = (poVM.count("gui") || ranges.empty() || (!poVM.count("save") && !poVM.count("export")));
  cout << endl << "running program with the following configuration:";
  cout << endl << " - directory: " << imgdir;
  cout << endl << " - ranges:    ";
  for (unsigned int i=0; i<ranges.size(); i+=2) cout << ranges[i] << "-" << ranges[i+1] << " ";
  cout << endl << " - GPS-mapping  " << (gps ? "true" : "false");
  cout << endl << " - adeption:    " << (adept ? "true" : "false");
  cout << endl << " - deskewing:   " << (deskew ? "true" : "false");
  cout << endl << " - incremental: " << (remOldMap ? "true" : "false");
  cout << endl << " - localMap:    " << (remFarMap ? "true" : "false");
  cout << endl << " - filtering:   " << (filter ? "true" : "false");
  cout << endl << " - resolution:  " << (double(mapResolCM))/100.0;
  cout << endl << endl;

  // ------------------------------------------------
  // -------------- run program ---------------------
  // ------------------------------------------------
  DataSetReaderVeloSlam dataSetReaderSlam(imgdir);
  DataSetReaderVeloTracking dataSetReaderTracking(imgdir);
  cout << endl << "SlamReader detected " << dataSetReaderSlam.getFrameCount() << " frames";
  cout << endl << "TrackingReader detected " << dataSetReaderTracking.getFrameCount() << " frames";
  cout << endl;
  DataSetReader *dataSetReader = dataSetReaderSlam.getFrameCount() > 0 ? (DataSetReader*)&dataSetReaderSlam : (DataSetReader*)&dataSetReaderTracking;
  Mapper mapper(dataSetReader, (double(mapResolCM))/100.0, maxDist, exportfile);

  // process requested frames
  if (!ranges.empty()) {
    unsigned int maxIdx = dataSetReader->getFrameCount()-1;
    for (unsigned int i=0; i<ranges.size(); i+=2) { // loop over defined ranges
      ranges[i] = min(maxIdx,ranges[i]);
      ranges[i+1] = min(maxIdx,ranges[i+1]);
      for (unsigned int idx=ranges[i]; idx<=ranges[i+1]; idx++) { // look over images in current range
        mapper.readScan(idx);
        if (!gps) mapper.registerScan(nbRegSubsamples, deskew);
        mapper.addScan(adept, deskew, remOldMap, remFarMap, gps);
      }
    }
    if (filter)
      mapper.filterMap();
    if (poVM.count("save")) {
      cout << endl << "saving to file " << savefile << "..." << flush;
      mapper.saveMap(savefile);
      cout << "done" << flush;
    }
  }

  // show gui if requested
  if (showGUI) {
    Gui3DQt::Gui myGui("VelodyneMapper", argc, argv);
    //myGui.registerVisualizer(new VisualizerCamControl(*myGui.getQGlWidget()));
    Visualizer3DMapper *mapperVis = new Visualizer3DMapper(mapper);
    mapperVis->ui.hsSubsampling->setValue(nbRegSubsamples);
    mapperVis->ui.cbAdept->setChecked(true);
    mapperVis->ui.cbUnwarp->setChecked(true);
    mapperVis->ui.cbRemOld->setChecked(remOldMap);
    mapperVis->ui.cbRemFar->setChecked(remFarMap);
    mapperVis->ui.cbGPS->setChecked(gps);
    if (!ranges.empty())
      mapperVis->ui.hsImageNb->setValue(ranges.back()+1); // both ranges/ui-indexing starts at 0, so mark next image
    sleep(1);
    myGui.registerVisualizer((Gui3DQt::Visualizer*)mapperVis, "mapper");
    cout << endl << endl << "bringing up GUI..." << flush; sleep(1);
    myGui.exec();
  }

  cout << endl;
  return 0;
}
