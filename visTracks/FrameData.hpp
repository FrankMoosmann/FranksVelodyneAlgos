/*
 * FrameData.hpp
 *
 *  Created on: Mar 20, 2012
 *      Author: moosmann
 */
#ifndef FRAMEDATA_H_
#define FRAMEDATA_H_

#include <vector>
#include <LidarImage.hpp>
#include <MatrixDefs.hpp>
#include <LidarImageProjector.hpp>
#include "TrackData.hpp"

class FrameData {
public:
  FrameData(std::string imgFile, std::string trackfile, std::string trackBBfile, const LidarImageProjector &proj, bool loadPoints);
  virtual ~FrameData();

  std::string imgFile;
  LidarImage<double> distanceImg;
  LidarImage<matrixTools::DVector> pointcloud;
  std::vector<TrackData> tracks;

  size_t getMemUsage();
};

#endif /* FRAMEDATA_H_ */
