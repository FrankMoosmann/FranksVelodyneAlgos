#include "PNGReader.hpp"

#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>

#include <kogmo_time.h>
#include <FMUtils.hpp>
#include <LidarImageProjectorPNG.hpp>
#include <PngDistanceImage.hpp>
#include <CsvReader.hpp>
#include <TrackingUtils.hpp> //distanceFromPng()

#include "ParameterHeap.hpp"

using namespace std;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;
using namespace boost::filesystem;

const size_t insFirstLineNb = 8; // the 8th line (i.e. at index 7) contains the first data row

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////        PNGReader        ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

PNGReader::PNGReader(string dir)
{
  try {
    determineFilesExtractV1(dir);
  } catch (exception) {
    determineFilesExtractV2(dir);
  }
  cout << " -> found " << frameSources.size() << " files" << endl;
}

PNGReader::~PNGReader()
{
}

void PNGReader::determineFilesExtractV1(std::string dir)
{
  cout << "searching for distance images in " << dir << endl;
  const path directory(dir);
  const boost::regex scanfilter( ".*\\.png" ); // "\\" is transformed into "\" at compile-time
  const boost::regex maskfilter( "mask\\.png" ); // "\\" is transformed into "\" at compile-time
  const boost::regex intensfilter( "intens.*\\.png" ); // "\\" is transformed into "\" at compile-time
  boost::smatch what; // match result
  list<string> dstImgFiles;
  if (exists(directory)) {
    directory_iterator end;
    for (directory_iterator iter(directory); iter != end; ++iter) {
      path currFile = *iter;
      if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
      if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
      if (!boost::regex_match(currFile.filename(), what, scanfilter)) continue; // Skip if no match
      if (boost::regex_match(currFile.filename(), what, maskfilter)) continue; // Skip if no match
      if (boost::regex_match(currFile.filename(), what, intensfilter)) continue; // Skip if no match
      dstImgFiles.push_back(currFile.string()); // File matches, store it
    }
  }
  dstImgFiles.sort(); // sort alphabetically

  CsvReader insData(dir+"/imu.cfg", ";", "", true); // file, separator, extractEmptyLines
  if ((insData.lineCount() < dstImgFiles.size()+insFirstLineNb-1) || (insData.lineCount() > dstImgFiles.size()+insFirstLineNb))
    cerr << "mismatch between number of images (" << dstImgFiles.size() << ") and number of lines in IMU data file ("<< insData.lineCount()+1-insFirstLineNb <<")" << endl;

  unsigned int lineNb = insFirstLineNb;
  BOOST_FOREACH(string dstImgFile, dstImgFiles) {
    FrameSource source;
    source.dstImg = dstImgFile;
    source.intImg = "";
    // read IMU file,  0: time, 5-7: yaw,pitch,roll(RAD), 14-16: x,y,z
    if ((insData.lineCount() > lineNb) && (insData.fieldCount(lineNb) > 16)) {
      source.frameStartTime = KogniMobil::kogmo_timestamp_from_string(insData.get(lineNb,0).c_str());
      source.positionHTM2w = YawPitchRollXYZ_2_HTM(insData.get<double>(lineNb,5), insData.get<double>(lineNb,6), insData.get<double>(lineNb,7), insData.get<double>(lineNb,14), insData.get<double>(lineNb,15), insData.get<double>(lineNb,16));
    }
    frameSources.push_back(source);
    ++lineNb;
  }

  imgCfgFile = dir+"/img.cfg";
}

void getFiles(const path &dir, const boost::regex &filter, vector<string> &targetlist)
{
  directory_iterator end;
  boost::smatch what; // match result
  for (directory_iterator iter(dir); iter != end; ++iter) {
    path currFile = *iter;
    if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
    if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
    if (!boost::regex_match(currFile.filename(), what, filter)) continue; // Skip if no match
    targetlist.push_back(currFile.string()); // File matches, store it
  }
}

void PNGReader::determineFilesExtractV2(std::string dir)
{
  cout << "searching for distance images in " << dir << endl;
  const boost::regex pngfilter( ".*\\.png" ); // "\\" is transformed into "\" at compile-time
  const boost::regex txtfilter( ".*\\.txt" ); // "\\" is transformed into "\" at compile-time
  vector<string> dstImgFiles;
  vector<string> intImgFiles;
  //vector<string> oxtFiles;
  if (!exists(path(dir) / "velodyne_dist" / "data" ))
    throw invalid_argument("directory does not exist");
  path timestampsfile = path(dir) / "velodyne_dist" / "timestamps_start.txt" ;
  getFiles(path(dir) / "velodyne_dist" / "data" , pngfilter, dstImgFiles);
  try {
    getFiles(path(dir) / "velodyne_intens" / "data" , pngfilter, intImgFiles);
  } catch (exception &e) {
    // intensity-directory is optional
  }
  //getFiles(path(dir) / "oxts" / "data" , txtfilter, oxtFiles);
  sort(dstImgFiles.begin(), dstImgFiles.end()); // sort alphabetically
  sort(intImgFiles.begin(), intImgFiles.end()); // sort alphabetically
  //sort(oxtFiles.begin(), oxtFiles.end()); // sort alphabetically
  CsvReader timeData(timestampsfile.string(), ";", "", false); // file, separator, extractEmptyLines

  if (   (dstImgFiles.size() != timeData.lineCount()))
    //  || (dstImgFiles.size() != oxtFiles.size())
    //  || (dstImgFiles.size() != intImgFiles.size()))
    throw range_error("unequal number of files/entries");

  for (size_t i=0; i<dstImgFiles.size(); ++i) {
    //CsvReader insData(oxtFiles[i], " ", "", false); // file, separator, extractEmptyLines
    FrameSource source;
    source.dstImg = dstImgFiles[i];
    source.intImg = intImgFiles.size()>i ? intImgFiles[i] : "";
    source.frameStartTime = KogniMobil::kogmo_timestamp_from_string(timeData.get(i,0).c_str());
    //source.positionHTM2w = YawPitchRollXYZ_2_HTM(insData.get<double>(lineNb,5), insData.get<double>(lineNb,6), insData.get<double>(lineNb,7), insData.get<double>(lineNb,14), insData.get<double>(lineNb,15), insData.get<double>(lineNb,16));
    frameSources.push_back(source);
  }

  imgCfgFile = dir+"/velodyne_dist/velodyne_config.txt";
}

void PNGReader::readImage(size_t fileIdx, FrameSPtr cf, FrameSPtr lf)
{
  using namespace HomogeneousTransformationMatrix;

  const FrameSource &source = frameSources[fileIdx];

  PngDistanceImage dstImg(source.dstImg);
  distanceFromPng(cf->distance, dstImg);
  // TODO: if exists(frameSources[fileIdx].intImg) then load intensities
  cf->intensity.fill(0.0f);
  cf->recTime = source.frameStartTime;
  cf->positionHTM2w = source.positionHTM2w;
  cf->positionHTM2v = invert_HTM(cf->positionHTM2w);
  cf->egoEstimationHTM2w = cf->positionHTM2w;
  cf->sourcefile = source.dstImg;

  if (lf->valid) {
    double deltaTSec = (double)(cf->recTime - lf->recTime) / 1000000000.0;
    if (deltaTSec < 0.00001) {
      cout << endl << "warning: no time delay between frames!" << flush;
      deltaTSec = ParameterHeap::get()->defaultDeltaT; // set default 10Hz if invalid
    }
    cf->speed = sqrt(pow(cf->positionHTM2w(0,3)-lf->positionHTM2w(0,3),2)+pow(cf->positionHTM2w(1,3)-lf->positionHTM2w(1,3),2)+pow(cf->positionHTM2w(2,3)-lf->positionHTM2w(2,3),2)) / deltaTSec;
    cf->timeDiffSec = deltaTSec;
  } else {
    cf->speed = 0.0;
    cf->timeDiffSec = 0.0;
  }

  cf->valid = true;
}


