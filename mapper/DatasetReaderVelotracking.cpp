#include "DatasetReaderVelotracking.hpp"

#include <limits>
#include <iostream>
#include <exception>
#include <algorithm>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string/replace.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::string;
using namespace boost::filesystem;

DataSetReaderVeloTracking::DataSetReaderVeloTracking(std::string baseDir)
    : baseDirectory(baseDir)
{
    scanForImages();
    initTimestamps();
    frameCount = std::numeric_limits<int>::max();
    frameCount = std::min(frameCount, (int)(distanceTimestamps.size()));
    frameCount = std::min(frameCount, (int)(distancePngFiles.size()));
    frameCount = std::min(frameCount, (int)(intensityPngFiles.size()));
}

DataSetReaderVeloTracking::~DataSetReaderVeloTracking()
{
}

void DataSetReaderVeloTracking::scanForImages()
{
    const path distDataDirectory = path(baseDirectory) / "velodyne_dist" / "data";
    cout << endl << "searching for distance images in " << distDataDirectory.string() << endl;
    const boost::regex scanfilter( ".*\\.png" ); // "\\" is transformed into "\" at compile-time
    boost::smatch what; // match result
    if (exists(distDataDirectory)) {
      directory_iterator end;
      for (directory_iterator iter(distDataDirectory); iter != end; ++iter) {
        path currFile = *iter;
        if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
        if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
        if (!boost::regex_match( currFile.leaf().string(), what, scanfilter)) continue; // Skip if no match
        distancePngFiles.push_back(currFile.string());
      }
    }
    sort(distancePngFiles.begin(), distancePngFiles.end()); // sort alphabetically

    for (std::string filename : distancePngFiles) {
        boost::replace_all(filename, "velodyne_dist", "velodyne_intens");
        path pathIntens = path(filename); //distDataDirectory / filename;
        if ((!exists(pathIntens)) || (!is_regular_file(pathIntens)) || (is_directory(pathIntens)))
          intensityPngFiles.push_back("");
        else
          intensityPngFiles.push_back(pathIntens.string());
    }

    cout << " -> found " << distancePngFiles.size() << " files";
}

void DataSetReaderVeloTracking::initTimestamps()
{
  std::string timestampsFilename = (path(baseDirectory) / "velodyne_dist"/ "timestamps.txt").string();
  std::ifstream timestampsfile(timestampsFilename.c_str());
  std::string line;
  while (std::getline(timestampsfile, line))
  {
     distanceTimestamps.push_back(line);
  }
//  std::copy(std::istream_iterator<std::string>(timestampsfile),
//              std::istream_iterator<std::string>(),
//              std::back_inserter(distanceTimestamps));
  cout << " -> found " << distanceTimestamps.size() << " timestamps" << endl;
}

std::string DataSetReaderVeloTracking::getBaseDirectory() const
{
    return baseDirectory;
}

std::string DataSetReaderVeloTracking::getPngConfigFilename() const
{
    return (path(baseDirectory) / "velodyne_dist"/ "velodyne_config.txt").string();
}

int DataSetReaderVeloTracking::getFrameCount() const
{
    return frameCount;
}

DataSetReaderVeloTracking::ImuData DataSetReaderVeloTracking::getImuData(int index) const
{
    ImuData result;
    return result;
}

DataSetReader::PngFile DataSetReaderVeloTracking::getDistancePng(int index) const
{
    PngFile result;
    if (index < distancePngFiles.size() && index < distanceTimestamps.size()) {
        result.filename = distancePngFiles[index];
        result.isValid = result.filename.length() > 0;
        result.timestamp = KogniMobil::kogmo_timestamp_from_string(distanceTimestamps[index].c_str());
        result.isTimestampValid = result.timestamp > 0;
    }
    return result;
}

DataSetReader::PngFile DataSetReaderVeloTracking::getIntensityPng(int index) const
{
    PngFile result;
    if (index < intensityPngFiles.size() && index < distanceTimestamps.size()) {
        result.filename = intensityPngFiles[index];
        result.isValid = result.filename.length() > 0;
        result.timestamp = KogniMobil::kogmo_timestamp_from_string(distanceTimestamps[index].c_str());
        result.isTimestampValid = result.timestamp > 0;
    }
    return result;
}
