#ifndef LIDARIMAGEFEATURES_H_
#define LIDARIMAGEFEATURES_H_

#include <cfloat>
#include <cmath>

#include "MatrixDefs.hpp"
#include "LidarImage.hpp"
#include "LidarImageProjector.hpp"

namespace LidarImageFeatures {

enum Neighborhood { NBH4, NBH8, NBH24 };

class NormalConfidenceLookup {
public:
  static NormalConfidenceLookup* get(double resolHAngRAD, double resolVAngRAD, double stdDevDist, double stdDevHAngRAD, double stdDevVAngRAD);
  //! returns the expected standard deviation of the normal vector's angle in randian
  double getStdDevRAD(double dist, double distHoriz, double distVert) const;
  //! returns the expected covariance of the normal vector's coordinates (assuming yaw/pitch = 0/0 of the center pixel)
  matrixTools::DMatrix getCovar(double dist, double distHoriz, double distVert) const;
private:
  NormalConfidenceLookup(double resolHAngRAD, double resolVAngRAD, double stdDevDist, double stdDevHAngRAD, double stdDevVAngRAD);
  ~NormalConfidenceLookup();
  bool loadFromFile(std::string filename = "normalConfidence.lookup");
  bool saveToFile(std::string filename = "normalConfidence.lookup");
  void saveToTextFile(std::string filename); //!< saves data to a text-file for plotting (values are in degrees!)

  static NormalConfidenceLookup* uniqueInstance; // volatile -> might be changed by other thread

  double resolHAngRAD;
  double resolVAngRAD;
  double stdDevDist;
  double stdDevHAngRAD;
  double stdDevVAngRAD;
  double tableDistance; //!< reference distance at which the lookup table was generated in meter
  double lrange; //!< lower range of lookup-table distance interval in meter
  double urange; //!< upper range of lookup-table distance interval in meter
  unsigned int lookupsize; //!< size of the lookup-table = number of entries per dimension
  double stepsize; //!< angle increment between entries in the lookup-table in rad
  double *lookuptableAngle; //!< the lookup-table, values are in rad
  // As a covariance is always symmetric it would be sufficient to store 6 values only
  double *lookuptableCovar; //!< the lookup-table of the 3D covariance matrix, all 9 values flattened

  // temporary variables that are used frequently: allocated only once
  mutable double fac, dhRel, dvRel;
  mutable unsigned int ih, iv;
  mutable matrixTools::DVector center, top, left, normal;
  mutable matrixTools::DMatrix covarTop3, covarLeft3, Samples, SamplesT;
  mutable matrixTools::DVector add;
};

/*!
 * \brief method to post-process a distance image
 *
 * This function interpolates missing pixels and smoothes the values
 *
 * \param distanceImage the image to work on
 * \param maxDist remove measurements with distances above this threshold
 * \param horizInterpolPixThresh don't interpolate over more than this number of pixels, horizontally
 * \param vertInterpolPixThresh don't interpolate over more than this number of pixels, vertically
 * \param horizInterpolDistDiffThresh don't interpolate if distance difference exceeds this threshold, horizontally/vertically
 * \param smooth if set true, image will be smoothed after interpolation
 * \param smoothMaxDist distance-threshold in meter for using neighboring pixel for smoothing
 * \param tmpBuffer a temporary buffer with same size as distanceImage. If not given, it will be temporarily allocated within the method
 */
void postprocessDistanceImage(LidarImage<double> &distanceImage, double maxDist = DBL_MAX, int horizInterpolPixThresh = 2, int vertInterpolPixThresh = 2, double horizInterpolDistDiffThresh = 0.5, bool smooth = true, double smoothMaxDist = 0.2, LidarImage<double> *tmpBuffer = NULL);

/*!
 * \brief method to transform a distance image into a 3D point cloud
 *
 * \param distanceImage the image used as source. Pixels will be set invalid, if projection fails
 * \param point3D the target 3D point cloud represented as 2D array of 3D-Vectors
 * \param projector the projector used to project values
 * \param pointVariance3D if not NULL covariance matrices of the 3D points will be stored. needs the following three parameters
 * \param stdDevDist standard deviation of distance measurements, used when normalVariance3D is not NULL
 * \param stdDevHAngRAD standard deviation of the angles of measurements in horizontal direction, used when normalVariance3D is not NULL
 * \param stdDevVAngRAD standard deviation of the angles of measurements in vertical direction, used when normalVariance3D is not NULL
 */
void transformTo3D(LidarImage<double> &distanceImage, LidarImage<matrixTools::DVector> &point3D, const LidarImageProjector &projector, LidarImage<matrixTools::DMatrix> *pointVariance3D = NULL,
    double stdDevDist = 0.0, double stdDevHAngRAD = 0.0, double stdDevVAngRAD = 0.0);

/*!
 * \brief method to calculate the difference between two images
 * In case a value is invalid (DBL_MAX), the returned difference will also be DBL_MAX
 *
 * \param first one image
 * \param second another image
 * \param result result of first - second
 */
void difference(const LidarImage<double> &first, const LidarImage<double> &second, LidarImage<double> &result);

/*!
 * \brief method to calculate derivatives of the frame
 * 
 * The method reads the "distanceImage" image and writes "derivH" and "derivV"
 * The derivative values are - so far - simple differences (in meter). This might be changed
 * either to have distance values comparable to the euclidean distance of the 3D coordinates
 * or to have distance-difference / angle-difference.
 * In case distance values are invalid (DBL_MAX), the returned derivatives will also be DBL_MAX
 * 
 * \param distanceImage image holding distance values of size (h x v)
 * \param derivH horizontal derivatives will be stored herein, size (h-1 x v)
 * \param derivV vertical derivatives will be stored herein, size (h x v-1)
 */
void derivative(const LidarImage<double> &distanceImage, LidarImage<double> &derivH, LidarImage<double> &derivV);

/*!
 * \brief method to weight connections of a distance image
 * 
 * The method weights connections between pixels, by comparing the distance values to the neighbors
 * Based on the uncertainty of the scanner it tries to figure out if the sampling theorem was violated or not
 * The returned weights are between 0.0 (do not use this connection, also in case of invalid measurements) and 1.0 (valid neighbor).
 * 
 * \param distanceImage the distance readings in meter of size (h x v)
 * \param derivH the horizontal derivative of the distance readings, size (h-1 x v)
 * \param derivV the vertical derivative of the distance readings, size (h x v-1)
 * \param connWeightsH the returned weights for horizontal neighbors, size (h-1 x v)
 * \param connWeightsV the returned weights for vertical neighbors, size (h x v-1)
 * \param MAX_DST_PM don't use connections with a maximum distance difference (per distance) of this value
 * \param MAX_DST_FAC don't use connections with a distance-factor higher than this value
 */
void connectionWeigths(const LidarImage<double> &distanceImage, const LidarImage<double> &derivH, const LidarImage<double> &derivV, LidarImage<double> &connWeightsH, LidarImage<double> &connWeightsV, double MAX_DST_PM, double MAX_DST_FAC, double A, double B, double C, double D);

/*!
 * \brief method to calculate normal vectors
 * 
 * Based upon the 3D points normals vectors are calculated using a 4-neighborhood.
 * 
 * \param distanceImage distance readings
 * \param points 3d point coordinates
 * \param normals the returned normal vectors, (0,0,0) if invalid, otherwise euclidean length = 1
 * \param connWeightsH if pointer valid, horizontal connection weights will be used
 * \param connWeightsV if pointer valid, vertical connection weights will be used
 * \param normalStdDevRAD if pointer valid, an estimate of the standard deviation of the normal vectors' angle will also be returned
 * \param normalVariance3D if point valid, the 3x3 covariance matrix of the normal vectors will be calculated
 * \param normalConfidence if pointer valid, a confidence estimate will also be returned, 0=invalid, 0=valid
 * \param tmpVec3D if pointer valid, this temporary buffer will decrease computation time
 * \param tmpDbl if pointer valid, this temporary buffer will decrease computation time
 * \param projector needed to estimate normalVariance3D
 * \param DIST_DIFF if neighbor is relatively further away than this threshold, use one of the next constants to weight connections. This can be used to enforce outer edges (because in this case the connection weighting will be low)
 * \param W_FAC_CREATE connection-weight on creation (see DIST_DIFF)
 * \param resolHAngRAD angular resolution of image in horizontal direction
 * \param resolVAngRAD angular resolution of image in vertical direction
 * \param stdDevDist standard deviation of distance measurements
 * \param stdDevHAngRAD standard deviation of the angles of measurements in horizontal direction
 * \param stdDevVAngRAD standard deviation of the angles of measurements in vertical direction
 * \param nbNConfMedianPasses number of median passes across normal confidence image
 * \param nConfNeighb type of neighborhood used for median on normal confidence
 * \param nConfMin if true the median is assigned via the minimum operator with the previous normal confidence values
 */
void normals(const LidarImage<double> &distanceImage, const LidarImage<matrixTools::DVector> &points,
    LidarImage<matrixTools::DVector> &normals, const LidarImage<double> *connWeightsH = NULL,
    const LidarImage<double> *connWeightsV = NULL, LidarImage<double> *normalStdDevRAD = NULL, LidarImage<matrixTools::DMatrix> *normalVariance3D = NULL,
    LidarImage<double> *normalConfidence = NULL, LidarImage<matrixTools::DVector> *tmpVec3D = NULL, LidarImage<double> *tmpDbl = NULL,
    const LidarImageProjector *projector = NULL,
    double DIST_DIFF = 0.03, double W_FAC_CREATE = 2.0, /*double W_FAC_SMOOTH = 2.0,*/ double resolHAngRAD = 0.007,
    double resolVAngRAD = 0.007, double stdDevDist = 0.015, double stdDevHAngRAD = 0.0015, double stdDevVAngRAD = 0.00017,
    unsigned int nbNConfMedianPasses = 1, Neighborhood nConfNeighb = NBH4, bool nConfMin = false, bool verbose = true);

/*!
 * \brief method to median-filter an image
 *
 * each pixel gets assigned the median of itself and the neighboring 4 pixels
 * the first version is an optimized version for 4-neighborhood
 *
 * \param data data that shall be filtered
 * \param tmpBuffer temporary buffer
 * \param invalidVal value considered as invalid, thus it is ignored
*/
void median(const LidarImage<double> &data, LidarImage<double> &target, double invalidVal = DBL_MAX);
void median(const LidarImage<double> &data, LidarImage<double> &target, Neighborhood nbh, double invalidVal = DBL_MAX);

/*!
 * \brief method to take element-wise maximum
 *
 * it is save to use as target the same as data1/data2
*/
template <typename T>
void maximum(const LidarImage<T> &data1, const LidarImage<T> &data2, LidarImage<T> &target) {
  assert(data1.getHorizSize() == data2.getHorizSize());
  assert(data1.getVertSize()  == data2.getVertSize() );
  assert(data1.getHorizSize() == target.getHorizSize());
  assert(data1.getVertSize()  == target.getVertSize() );
  typename LidarImage<T>::const_iterator d1It = data1.begin();
  typename LidarImage<T>::const_iterator d1End = data1.end();
  typename LidarImage<T>::const_iterator d2It = data2.begin();
  typename LidarImage<T>::const_iterator d2End = data2.end();
  typename LidarImage<T>::iterator dtIt = target.begin();
  typename LidarImage<T>::iterator dtEnd = target.end();
  while ((d1It != d1End) && (d2It != d2End) && (dtIt != dtEnd)) {
    *dtIt = max(*d1It, *d2It);
    ++d1It;
    ++d2It;
    ++dtIt;
  }
}

/*!
 * \brief method to take element-wise minimum
 *
 * it is save to use as target the same as data1/data2
*/
template <typename T>
void minimum(const LidarImage<T> &data1, const LidarImage<T> &data2, LidarImage<T> &target) {
  assert(data1.getHorizSize() == data2.getHorizSize());
  assert(data1.getVertSize()  == data2.getVertSize() );
  assert(data1.getHorizSize() == target.getHorizSize());
  assert(data1.getVertSize()  == target.getVertSize() );
  typename LidarImage<T>::const_iterator d1It = data1.begin();
  typename LidarImage<T>::const_iterator d1End = data1.end();
  typename LidarImage<T>::const_iterator d2It = data2.begin();
  typename LidarImage<T>::const_iterator d2End = data2.end();
  typename LidarImage<T>::iterator dtIt = target.begin();
  typename LidarImage<T>::iterator dtEnd = target.end();
  while ((d1It != d1End) && (d2It != d2End) && (dtIt != dtEnd)) {
    *dtIt = std::min(*d1It, *d2It);
    ++d1It;
    ++d2It;
    ++dtIt;
  }
}

////////////////////////////////////////////////////
///////////       Utiliy Function         //////////
////////////////////////////////////////////////////

//! returns 1 for x << 0, and 0 for x >> 0, and 0.5 for x = 0
inline double sigmoidLikeShiftedUp(double x)
{
  return 0.5f + ((-0.5f*x)/sqrt(1+x*x));
}
//! returns 1 for x << thresh, and 0 for x >> thresh, |narrowFac|>1 makes the threshold harder, negative narrowFac inverses result (0 for x<<thresh, 1 for x>>thresh)
inline double sigmoidLikeSoftThresh(double x, double thresh, double narrowFac, bool normalizeAtZero = false)
{
  x = (x-thresh)*narrowFac;
  if (normalizeAtZero)
    return sigmoidLikeShiftedUp(x)/sigmoidLikeShiftedUp(0.0f); // normalize so value f(0)=1
  else
    return sigmoidLikeShiftedUp(x);
}
/*! tries to estimate if two pixels/measurements belong to the same surface
 * \param d1 distance of one pixel of connection
 * \param d2 distance of other pixel of connection
 * \param diff distance-difference between the pixels
 * \param diffL, diffR distance differences of the neighboring connections
*/
inline double connectivity(double d1, double d2, double diff, double diffL, double diffR, double MAX_DST_PM, double MAX_DST_FAC, double A, double B, double C, double D) {
  if ((d1 == DBL_MAX) || (d2 == DBL_MAX) || (fabs(d1) < 0.01) || (fabs(d2) < 0.01))
    return 0.0;

//  if (fabs(diff) < 0.1f) // always accept a small connection, so only check further if sufficient large connection (avoids singularity if diffL/diffR are also small)
  if (fabs(diff) < D) // always accept a small connection, so only check further if sufficient large connection (avoids singularity if diffL/diffR are also small)
    return 1.0;

  if ((fabs(diffL) < 0.01f) || (fabs(diffR) < 0.01f)) // always reject connection if it is not small but one neighboring connection is (avoids division by zero)
    return 0.0;

  double dist = std::min(d1, d2);
//  const double nFac = 0.1;
  const double nFac = A - B*exp(-dist*C); // more relative noise if closer -> adjust inclination
  double keep = sigmoidLikeSoftThresh(fabs(diff)/dist, MAX_DST_PM, 2/MAX_DST_PM); // MAX_DST_PM -> angle of surface plane
  keep = fmin(keep, sigmoidLikeSoftThresh(fabs((diff-diffL)/diffL), MAX_DST_FAC, nFac));
  keep = fmin(keep, sigmoidLikeSoftThresh(fabs((diff-diffR)/diffR), MAX_DST_FAC, nFac));
  return keep;
};

} //! end namespace

#endif /*LIDARIMAGEFEATURES_H_*/
