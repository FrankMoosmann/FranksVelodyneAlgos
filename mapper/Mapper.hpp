#ifndef VISUALIZER3D
#define VISUALIZER3D

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <MatrixDefs.hpp>
#include <CsvReader.hpp>
#include <LidarImageProjector.hpp>

#include "Frame.hpp"
#include "Map.hpp"

class Mapper
{
public:
  Mapper(std::string imgDirectory, double mapResolution, double maxDist, const LidarImageProjector &proj, std::string rem3DptsFile);
  virtual ~Mapper();

  unsigned int getImgCount() const {return imageNames.size();};
  std::string getImgName(unsigned int idx) const {return imageNames[idx].first;};
  std::string getIntensName(unsigned int idx) const {return imageNames[idx].second;};
//  void processImage(unsigned int idx, unsigned int nbSamples, bool unwarp, bool adept);
//  void processImages(unsigned int startIdx, unsigned int endIdx, unsigned int nbSamples, bool unwarp, bool adept);
  void reset();
  void saveMap(std::string filename);
  void exportMap(std::string filename);
  void exportMapPt3D(std::string filename, unsigned int resolutionFactor = 1, Map::ExportType expType = Map::Points);
  void loadMap(std::string filename);
  void importMap(std::string filename);
  void readScan(unsigned int idx);
  void registerScan(unsigned int nbSamples, bool unwarp);
  void reg1(unsigned int nbSamples, bool unwarp);
  void reg2();
  void reg3();
  void filterMap();
  void filterMap_init();
  void filterMap_step();
  void filterMap_finish();
  void addScan(bool adept, bool unwarp, bool remOldMap, bool remFarPts, bool useGPSonly);

  Map* getMap() {return map;};
  LFrameSPtr getCurrentFrame() {return currentFrame;};
  LFrameSPtr getLastFrame() {return lastFrame;};

private:

  static const size_t insFirstLineNb = 8; // the 8th line (i.e. at index 7) contains the first data row
  const double resolution; //!< resolution of the map
  const double maxDist; //!< maximum scanning distance to use points from
  std::vector< std::pair<std::string, std::string> > imageNames; // pair of distance/intensity image names
  CsvReader *insData;

  const LidarImageProjector &proj;
  std::string remPtsFilename;
  bool remPtsFileIsBinary;
  std::ofstream *remPtsFile;
  bool remPtsUseIntensity;
  Map *map;
  LFrameSPtr currentFrame, lastFrame;
  double mercatorScale;
  double mx0, my0, mz0;
  int64_t timestampFirst;
  unsigned int lastFrameImsLineNb;
  int64_t timeDiffBuffer[9]; //!< indexed by nbAddedScans%9, updated in registration-routine
  int64_t avgFrameDuration; //!< median out of timeDiffBuffer, updated in registration-routine

  void mapRemovedPoint(double x, double y, double z, unsigned char intensity, unsigned int hitCnt);
  void writeLanLonFile(std::string mapfilename);

};

#endif // VISUALIZER3D
