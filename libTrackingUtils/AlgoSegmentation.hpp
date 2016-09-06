#ifndef ALGOSEGMENTATION_H_
#define ALGOSEGMENTATION_H_

#include "LidarFrame.hpp"


/*!
 * \brief method to segment a frame
 * 
 * The method reads several images of the frame and writes the "segmentPtr" image
 * "0" means no segment was assigned, otherwise a valid pointer to a segment is present.
 * All Segments are additionally stored in the segmentList
 * 
 * \param lf  The frame used for reading and storing
 */
void segment(FrameSPtr lf, bool verbose = true);

typedef TwoLayerGraph<const LidarSegment*, const LidarSegment*> SegmentationLinks;
SegmentationLinks linkSegments(const LidarImage<LidarSegment> &seg1, const LidarImage<LidarSegment> &seg2);
void saveLinks(const SegmentationLinks &links, std::string filename);
void linkAndSave(FrameSPtr frame1, FrameSPtr frame2);

//! \brief recolor segments
void recolorSegments(FrameSPtr lf);

////////////////////////////////////////////////////
///////////       Utiliy Function         //////////
////////////////////////////////////////////////////
inline void crossProduct(double ax, double ay, double az, double bx, double by, double bz, double &rx, double &ry, double &rz) {
  rx = ay*bz-az*by;
  ry = az*bx-ax*bz;
  rz = ax*by-ay*bx;
};
inline double dotProduct(double ax, double ay, double az, double bx, double by, double bz) {
  return (ax*bx) + (ay*by) + (az*bz);
};
inline double length(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
};
inline void normalize(double &x, double &y, double &z, double targetlength = 1.0)
{
  double f = targetlength / length(x,y,z);
  x *= f;
  y *= f;
  z *= f;
};

#endif /*ALGOSEGMENTATION_H_*/
