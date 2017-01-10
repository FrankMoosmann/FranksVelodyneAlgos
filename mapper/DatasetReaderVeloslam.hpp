#ifndef DATASETREADERVELOSLAM_H
#define DATASETREADERVELOSLAM_H

#include <string>
#include <CsvReader.hpp>
#include <kogmo_time.h>
#include <vector>

class DataSetReader
{
public:
    struct ImuData {
        bool isTimestampValid = false;
        KogniMobil::kogmo_timestamp_t timestamp = 0;
        bool isPoseValid = false;
        double latDeg;
        double lonDeg;
        double altMeter;
        double yawRad;
        double pitchRad;
        double rollRad;
    };
    struct PngFile {
        bool isValid = false;
        std::string filename;
        bool isTimestampValid = false;
        KogniMobil::kogmo_timestamp_t timestamp;
    };

    virtual std::string getBaseDirectory() const = 0;
    virtual std::string getPngConfigFilename() const = 0;
    virtual int getFrameCount() const = 0;
    virtual ImuData getImuData(int index) const = 0;
    virtual PngFile getDistancePng(int index) const = 0;
    virtual PngFile getIntensityPng(int index) const = 0;
    // getCameraImage() - with timestamp
};

class DataSetReaderVeloSlam : public DataSetReader
{
public:
  DataSetReaderVeloSlam(std::string baseDir);
  virtual ~DataSetReaderVeloSlam();

  virtual std::string getBaseDirectory() const;
  virtual std::string getPngConfigFilename() const;
  virtual int getFrameCount() const;
  virtual ImuData getImuData(int index) const;
  virtual PngFile getDistancePng(int index) const;
  virtual PngFile getIntensityPng(int index) const;

private:
  std::string baseDirectory;
  int frameCount;

  CsvReader *insReader;
  static const size_t insFirstLineNb = 8; // the 8th line (i.e. at index 7) contains the first data row
  std::vector<std::string> distancePngFiles;
  std::vector<std::string> intensityPngFiles;

  void initImuReader();
  void scanForImages();
};

#endif // DATASETREADERVELOSLAM_H
