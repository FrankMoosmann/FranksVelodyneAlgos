#ifndef DATASETREADERVELOTRACKING_H
#define DATASETREADERVELOTRACKING_H

#include <vector>

#include "DatasetReader.hpp"

class DataSetReaderVeloTracking : public DataSetReader
{
public:
  DataSetReaderVeloTracking(std::string baseDir);
  virtual ~DataSetReaderVeloTracking();

  virtual std::string getBaseDirectory() const;
  virtual std::string getPngConfigFilename() const;
  virtual int getFrameCount() const;
  virtual ImuData getImuData(int index) const;
  virtual PngFile getDistancePng(int index) const;
  virtual PngFile getIntensityPng(int index) const;

private:
  std::string baseDirectory;
  int frameCount;

  std::vector<std::string> distanceTimestamps;
  std::vector<std::string> distancePngFiles;
  std::vector<std::string> intensityPngFiles;

  void scanForImages();
  void initTimestamps();
};

#endif // DATASETREADERVELOTRACKING_H
