#ifndef LIDARFRAMEBUFFER_H_
#define LIDARFRAMEBUFFER_H_

#include <vector>
#include "LidarFrame.hpp"

/*!
 * \class LidarFrameBuffer
 * \brief This class efficiently buffers several frames, to allow for some history
 */
class LidarFrameBuffer
{
public:
  LidarFrameBuffer(const unsigned int buffersize, const unsigned int imgHorizSize, const unsigned int imgVertSize);
  virtual ~LidarFrameBuffer();
  
  /*! returns one of the frames buffered
   * \param relPastFrame 0..buffersize-1 to get already captured data (0=current), negative values to get future buffer
   * \throws range_error if relPastFrame is >= buffersize
   */
  FrameSPtr get(int relPastFrame = 0);
  FrameSPtr next(); //!< increments pointer, can be called to store new frame
  unsigned int size() const {return frames.size();}; 
  
private:
  unsigned int current;
  std::vector<FrameSPtr> frames; // store timestamp along???
  
};

#endif /*LIDARFRAMEBUFFER_H_*/
