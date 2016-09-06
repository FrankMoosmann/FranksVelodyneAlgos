/*
 * FrameData.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: moosmann
 */

#include "FrameData.hpp"

#include <PngDistanceImage.hpp>
#include <TrackingUtils.hpp>
#include <LidarImageFeatures.hpp>

#include <fstream>
#include <iostream>

using namespace std;

bool compare(const TrackData& t1, const TrackData& t2) {
  return (t1.uid < t2.uid);
}

#define TRYCATCH(a,b) try{a} catch(exception){b}


FrameData::FrameData(std::string imgFile_, std::string trackfile_, std::string trackBBfile_, const LidarImageProjector &proj, bool loadPoints )
  : imgFile(imgFile_)
  , distanceImg(proj.getImgHorizSize(),proj.getImgVertSize())
  , pointcloud(proj.getImgHorizSize(),proj.getImgVertSize())
{
  PngDistanceImage dstImg(imgFile_);
  distanceFromPng(distanceImg, dstImg);
  LidarImageFeatures::transformTo3D(distanceImg, pointcloud, proj);

  ifstream tfile(trackfile_.c_str());
  ifstream bfile;
  TRYCATCH(bfile.open(trackBBfile_.c_str());,)
  string lineT,lineB;
  while (!tfile.eof()) {
    getline(tfile,lineT);
    TRYCATCH(getline(bfile,lineB);,lineB="";)
    if (lineT.length() > 0) {
      tracks.push_back(TrackData(lineT, lineB, loadPoints));
    }
  }
  sort(tracks.begin(), tracks.end(), compare);
  //cout << endl << "created " << imgFile << flush;
}

FrameData::~FrameData()
{
  //cout << endl << "freeing " << imgFile << flush;
}

size_t FrameData::getMemUsage() {
  size_t s = 0;
  s += distanceImg.getMemorySize();
  s += pointcloud.getMemorySize();
  BOOST_FOREACH(const TrackData &track, tracks)
    s += track.getMemUsage();
  return s;
}
