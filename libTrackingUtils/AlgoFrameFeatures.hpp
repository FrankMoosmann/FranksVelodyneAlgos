#ifndef ALGOFRAMEFEATURES_H_
#define ALGOFRAMEFEATURES_H_

#include <cfloat>
#include <cmath>

#include "LidarFrame.hpp"
#include "LidarImageProjector.hpp"

/*!
 * \brief method to post-process pure distance data
 *
 * Initially a new frame (with valid-flag set to true)
 * only has to hold correct distance/position data.
 * This function post-processes this data by smoothing the values
 * (if enabled via ParameterHeap) and transforms them into 3D
 *
 * \param frame  The frame used for reading and storing
 */
void postprocessDistanceData(FrameSPtr lf, const LidarImageProjector &projector);

/*!
 * \brief method to calculate derivatives of the frame
 * 
 * The method reads the "distance" image and writes "distDerivativeH" and "distDerivativeV"
 * The derivative values are - so far - simple differences (in meter). This might be changed
 * either to have distance values comparable to the euclidean distance of the 3D coordinates
 * or to have distance-difference / angle-difference.
 * Invalid values result from invalid distance values
 * 
 * \param frame  The frame used for reading and storing
 */
void calculateDerivative(FrameSPtr frame);

/*!
 * \brief method to calculate difference of distance values relative to last frame
 * 
 * \param lf  The last frame (access: read-only)
 * \param cf  The current frame
 */
void calculateDistDiff(FrameSPtr lf, FrameSPtr cf);

/*!
 * \brief method to calculate connections of the frame
 * 
 * The method reads the "distance", "distDerivativeH" and "distDerivativeV" image
 * and writes "connectivityH" and "connectivityV" images.
 * The connection values are between 0.0 (no conn., also in case of invalid measurements) and 1.0 (conn.).
 * 
 * \param frame  The frame used for reading and storing
 */
void calculateConnections(FrameSPtr frame);

/*!
 * \brief method to calculate normals of the frame
 * 
 * The method reads the "distance", "distDerivativeH", "distDerivativeV", "connectivityH" and "connectivityV" image
 * and writes "normalX", "normalY" and "normalZ" images.
 * The normals have euclidean length = 1 or (0,0,0) if invalid
 * 
 * \param frame  The frame used for reading and storing
 * \param projector Object for projecting 3D points onto the image
 * \param full   set to false if covariances should not be calculated
 */
void calculateNormals(FrameSPtr frame, const LidarImageProjector &projector, bool full = true, bool verbose = true);

/*!
 * \brief method to calculate features for matching
 * 
 * The method reads the "distance" image and writes the "matchingFeatures" image
 * The written feature vectors are 4-dimensional with real-valued scalar values
 * For invalid measurements the FV will be empty
 * 
 * \param frame  The frame used for reading and storing
 */
void calculateMatchingFeatures(FrameSPtr frame);

/*!
 * \brief calculates the "visibleInNextFrame" image, indicating if a pixel is visible in the next frame, assuming a static scene
 */
void visibilityNextFrame(FrameSPtr cf, FrameSPtr nf, LidarImageProjector &projector);
void visibilityPreviousFrame(FrameSPtr cf, FrameSPtr pf, LidarImageProjector &projector);
void visibilityOtherFrame(FrameSPtr cf, FrameSPtr of, LidarImageProjector &projector, LidarImage<bool> &visibilityResult); // used by upper two methods


#endif /*ALGOFRAMEFEATURES_H_*/
