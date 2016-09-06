#include "LidarFrameBuffer.hpp"

#include <stdexcept>

using namespace std;

LidarFrameBuffer::LidarFrameBuffer(const unsigned int buffersize, const unsigned int imgHorizSize, const unsigned int imgVertSize)
{
  cout << "constructing frame buffer of size " << buffersize << flush;
  frames.reserve(buffersize);
  for (unsigned int i=0; i<buffersize; ++i) {
    FrameSPtr f(new LidarFrame(imgHorizSize, imgVertSize));
    if (i == 0)
       cout << " using " << f->memorySize()/1024 << " kB each..." << flush;
    frames.push_back(f);
  }
  cout << "done" << endl;
}

LidarFrameBuffer::~LidarFrameBuffer()
{
}

FrameSPtr LidarFrameBuffer::get(int relPastFrame)
{
  if (relPastFrame >= (int)frames.size()) throw range_error("LidarFrameBuffer::get: required frame is out of range");
  relPastFrame = (current+frames.size()-relPastFrame) % frames.size();
  return frames[relPastFrame];
}

FrameSPtr LidarFrameBuffer::next()
{
  current = (current+1) % frames.size();
  return frames[current];
}
