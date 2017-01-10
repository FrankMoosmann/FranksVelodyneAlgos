#include "DatasetReaderVeloslam.hpp"

#include <limits>
#include <iostream>
#include <exception>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/bind.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::string;

DataSetReaderVeloSlam::DataSetReaderVeloSlam(std::string baseDir)
    : baseDirectory(baseDir)
{
    initImuReader();
    scanForImages();
    frameCount = std::numeric_limits<int>::max();
    if (insReader)
        frameCount = std::min(frameCount, (int)(insReader->lineCount()-insFirstLineNb+1));
    frameCount = std::min(frameCount, (int)(distancePngFiles.size()));
    frameCount = std::min(frameCount, (int)(intensityPngFiles.size()));
}

DataSetReaderVeloSlam::~DataSetReaderVeloSlam()
{
    if (insReader) delete insReader;
}

void DataSetReaderVeloSlam::initImuReader()
{
    try {
      cout << endl << "searching for INS data..." << flush;
      insReader = new CsvReader(baseDirectory+"/imu.cfg", ";", "", true); // file, separator, comments-indicators, extractEmptyLines
      cout << " -> found " << insReader->lineCount()-insFirstLineNb+1 << " entries" << endl;
    } catch (std::exception &e) {
      insReader = NULL;
      cerr << "Error reading imu.cfg: " << e.what() << endl;
    }
}

void DataSetReaderVeloSlam::scanForImages()
{
    using namespace boost::filesystem;

    cout << endl << "searching for distance images in " << baseDirectory << endl;
    const path directory(baseDirectory);
    const boost::regex scanfilter( "scan.*\\.png" ); // "\\" is transformed into "\" at compile-time
    boost::smatch what; // match result
    if (exists(directory)) {
      directory_iterator end;
      for (directory_iterator iter(directory); iter != end; ++iter) {
        path currFile = *iter;
        if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
        if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
        if (!boost::regex_match( currFile.leaf().string(), what, scanfilter)) continue; // Skip if no match
        distancePngFiles.push_back(currFile.string());
      }
    }
    sort(distancePngFiles.begin(), distancePngFiles.end()); // sort alphabetically

    for (std::string filename : distancePngFiles) {
        filename.replace(0, 4, "intens"); // replace "scan" by "intens"
        path pathIntens = directory / filename;
        if ((!exists(pathIntens)) || (!is_regular_file(pathIntens)) || (is_directory(pathIntens)))
          intensityPngFiles.push_back("");
        else
          intensityPngFiles.push_back(pathIntens.string());
    }

    cout << " -> found " << distancePngFiles.size() << " files";
}

std::string DataSetReaderVeloSlam::getBaseDirectory() const
{
    return baseDirectory;
}

std::string DataSetReaderVeloSlam::getPngConfigFilename() const
{
    return baseDirectory + "/img.cfg";
}

int DataSetReaderVeloSlam::getFrameCount() const
{
    return frameCount;
}

DataSetReaderVeloSlam::ImuData DataSetReaderVeloSlam::getImuData(int index) const
{
    ImuData result;
    if (index >= frameCount) {
        return result;
    }
    int insLineNb = index + insFirstLineNb;
    // 0 - time
    // 1 - #sats
    // 2-4 - lat, long, alt (DEG, meter)
    // 5-7 - yaw,pitch,roll(RAD)
    // 8-13 - velocities, angular and translational
    // 14-16 - x,y,z
    if (insReader->fieldCount(insLineNb) > 0) { // record time exists
        result.timestamp = KogniMobil::kogmo_timestamp_from_string(insReader->get(insLineNb,0).c_str());
        result.isTimestampValid = true;
    }
    if (insReader->fieldCount(insLineNb) > 7) {
      result.latDeg = insReader->get<double>(insLineNb,2);
      result.lonDeg = insReader->get<double>(insLineNb,3);
      result.altMeter = insReader->get<double>(insLineNb,4);
      result.yawRad = insReader->get<double>(insLineNb,5);
      result.pitchRad = insReader->get<double>(insLineNb,6);
      result.rollRad = insReader->get<double>(insLineNb,7);
      result.isPoseValid = true;
    }
    return result;
}

DataSetReader::PngFile DataSetReaderVeloSlam::getDistancePng(int index) const
{
    PngFile result;
    if (index < distancePngFiles.size()) {
        result.filename = distancePngFiles[index];
        result.isValid = result.filename.length() > 0;
        DataSetReader::ImuData imuData = getImuData(index);
        result.timestamp = imuData.timestamp;
        result.isTimestampValid = imuData.isTimestampValid;
    }
    return result;
}

DataSetReader::PngFile DataSetReaderVeloSlam::getIntensityPng(int index) const
{
    PngFile result;
    if (index < intensityPngFiles.size()) {
        result.filename = intensityPngFiles[index];
        result.isValid = result.filename.length() > 0;
        DataSetReader::ImuData imuData = getImuData(index);
        result.timestamp = imuData.timestamp;
        result.isTimestampValid = imuData.isTimestampValid;
    }
    return result;
}
