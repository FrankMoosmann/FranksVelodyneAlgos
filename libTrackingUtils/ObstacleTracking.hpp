#ifndef OBSTACLETRACKING_H_
#define OBSTACLETRACKING_H_

#include <vector>
#include "MatrixDefs.hpp"
#include "LidarFrameBuffer.hpp"
#include "LidarImageProjector.hpp"
#include "AlgoRegistration.hpp"

/*!
 * \class ObstacleTracking
 * \brief This class implements object detection and tracking on range data
 */
class ObstacleTracking
{
public:
  ObstacleTracking(const unsigned int buffersize, const LidarImageProjector &imageProjector);
  virtual ~ObstacleTracking();
  
  // access functions:
  /*! returns one of the frames buffered
   * \param relPastFrame 0..buffersize-1 to get already captured data (0=current), negative values to get future buffer
   * \throws range_error if relPastFrame is >= buffersize
   */
  unsigned int bufferSize() const {return frameBuffer.size();};
  void reset(const unsigned int buffersize);
  FrameSPtr getFrame(int relPastFrame);
  void next(); //!< increments buffer pointer;

  // functions for frame-wise tracking:
  void processFrame(unsigned int relPastFrame = 0); //!< does everything a frame needs. for more control the single steps involved can also be executed separately

  // single step functions for processFrame:
  void calcFeatures(unsigned int relPastFrame = 0); //!< calculates the features needed by the next methods
  void calcFeatureMatches(unsigned int relPastFrame = 0); //!< calculates pixel-correspondences between two frames
  void calcRegistration(unsigned int relPastFrame = 0); //!< calculates and applies measurements for each track
  void calcSegmentation(unsigned int relPastFrame = 0); //!< split pointcloud into parts
  void calcTracks(unsigned int relPastFrame = 0); //!< manage tracks (create new, delete, merge, split) based upon segmentation

  // for debugging purposes made public:
  void registrationInit(unsigned int relPastFrame = 0);
  bool registrationNext(bool retFalseOnLayerChange = false); // returns false if nothing to be done
  TrackRegistration::SPtr registrator; // will point to frames that it was instantiated on (not neccessarily the current frame)
private:
  const LidarImageProjector &imageProjector; //provides coordinate mapping from 3D to 2D and back
  LidarFrameBuffer frameBuffer; //contains laser data, features, and tracks
};

#endif /*OBSTACLETRACKING_H_*/
