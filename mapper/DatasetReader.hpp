#ifndef DATASETREADER_H
#define DATASETREADER_H

#include <string>
#include <kogmo_time.h>

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

#endif // DATASETREADER_H
