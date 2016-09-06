#ifndef LFRAME_H_
#define LFRAME_H_

#include <vector>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <MatrixDefs.hpp>
#include <LidarImage.hpp>
#include "Surface.hpp"

class Frame;
typedef boost::shared_ptr< Frame > LFrameSPtr;

namespace mdefs = matrixTools;
/*! 
 * \class LidarImageFrame
 * \brief This is the core class containing all data for signal processing on one Lidar Frame
 * 
 * An object of this class represents one Lidar frame and all associated data. More specifically it holds
 * - The associated position of the car in Cartesian absolute coordinates (in form of a HTM)
 * - The Lidar data itself as (1) distance image, (3) absolute Cartesian images for x,y,z respectively
 * - The first derivative in horiz and vert direction of the distance values
 * - A "pixel connectivity" image
 * - Associated normal vectors
 * - And more....
 */
class Frame
{
public:
  Frame(const unsigned int horizSize, const unsigned int vertSize);
  virtual ~Frame();

public:
  bool                                   valid; //!< indicates if at least basic data (recTime, pos, speed, distance+intensity) is valid
  int64_t                               recTime; //!< timestamp of recording time (in nanoseconds)
  int64_t                               timeDiffToLast; //!< time difference to last scan
  mdefs::DMatrix                        positionHTM2v; //!< vehicle position as HTM. point in world cs pw can be transformed into vehicle cs by pv=HTM*pw
  mdefs::DMatrix                        positionHTM2w; //!< vehicle position as HTM. point in vehicle cs pv can be transformed into world cs by pw=HTM*pv
  mdefs::DMatrix                        insPositionHTM2w; //!< vehicle position as HTM. point in vehicle cs pv can be transformed into world cs by pw=HTM*pv
  mdefs::DMatrix                        moveHTM; //!< move during current scan
  mdefs::DMatrix                        diffHTMToLast; //!< move between last scan and current scan (equal to moveHTM if no frame was missing)
  LidarImage<double>                    distance; //!< distance in meter, DBL_MAX indicating invalid measurement
  LidarImage<unsigned char>             intensity; //!< intensity ranging from 0 to 255
  LidarImage<mdefs::DVector>            point3D; //!< relative ego-based coordinates in meter, (DBL_MAX,DBL_MAX,DBL_MAX) indicating invalid measurement (in that case distance is also DBL_MAX)
  LidarImage<double>                    distDerivativeH; //!< distance differences in meter, DBL_MAX if invalid
  LidarImage<double>                    distDerivativeV;
  LidarImage<double>                    connectivityH; //!< weighted connections, values between 0.0 and 1.0
  LidarImage<double>                    connectivityV;
  LidarImage<mdefs::DVector>            normal3D; //!< relative ego-based coordinates in meter, (0,0,0) indicating invalid measurement
  LidarImage<double>                    normalStdDevRAD; //!< the standard deviation of the angle of the normal vector in radian
  LidarImage<mdefs::DMatrix>            normalVariance3D; //!< the standard deviation of the angle of the normal vector in 3D coordinates
  LidarImage<double>                    normalConfidence; //!< a value between 0 (unreliable/not valid) and 1 (reliable) indicating how precise the normal3D estimate is. Can be interpreted as c*exp(-covar)
  LidarImage<Surface::SPtr>              surface; //!< pointers to surfaces of the map. stored here to be able to modify them later
  LidarImage<mdefs::DVector>            tmpVec3D;
  LidarImage<double>                    tmpFloat1;
  
  std::string getTime();

private: // hide, so they cannot be used
  Frame()  {};
  Frame& operator=(const Frame &other) {(void)other; return *this;};
  void reset(); // set default values to all entries

};

#endif /*LFRAME_H_*/
