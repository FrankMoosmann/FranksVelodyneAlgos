/*!
    \file   ParameterHeap.h
    \brief  Provides a class for storing all parameters
    \author  Frank Moosmann (<moosmann@mrt.uka.de>),
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef PARAMETERHEAP_H_
#define PARAMETERHEAP_H_

//#include <boost/thread.hpp>
#include <boost/utility.hpp>

#include "LidarImageFeatures.hpp"
#include "ParamStructCreation.hpp"
#include "MatrixDefs.hpp"

namespace mdefs = matrixTools;

#define PARAMHEAP_BASE_ARR_1 (22, ( \
       (std::string,    outDir,                         "",      "specify output directory, output disabled if empty", 0) \
      ,(double,         lidarHAngResolRAD,              0.007,   "angular resolution in horizontal direction in radian", 0) \
      ,(double,         lidarVAngResolRAD,              0.007,   "angular resolution in vertical direction in radian", 0) \
      ,(double,         framegenMaxDistance,            60,      "remove measurements with distances above this threshold", 0) \
      ,(int,            framegenInterpolHPixThresh,     5,       "Don't interpolate over more than this number of pixels, horizontally", 0) \
      ,(int,            framegenInterpolVPixThresh,     2,       "Don't interpolate over more than this number of pixels, vertically", 0) \
      ,(double,         framegenInterpolMetPixThresh,   0.5,     "Don't interpolate if distance difference exceeds this threshold, horizontally/vertically", 0) \
      ,(bool,           framegenSmoothDstImg,           true,    "Flag, whether to smooth/blur image or not (in order to filter noise)", 0) \
      ,(double,         framegenSmoothDstThresh,        0.2,     "Distance-Threshold in meter for using neighboring pixel for smoothing", 0) \
      ,(double,         connMaxDst,                     0.2,     "Maximum absolute distance difference allowed for pixel-connections (per meter displacement from the car)", 0) \
      ,(double,         connMaxDstRatio,                0.2,     "Maximum difference ratio for a connection allowed (relative to neighboring connections)", 0) \
      ,(double,         connFacA,                       1.0,     "", 0) \
      ,(double,         connFacB,                       1.0,     "", 0) \
      ,(double,         connFacC,                       0.02,    "", 0) \
      ,(double,         connFacD,                       0.1,     "", 0) \
      ,(double,         featureConnectWeightCreate,     2.0,     "default weight for creating normals of connections that are farther away", 0) \
      ,(unsigned int,   featurenNbConfMedianPasses,     1,       "number of median-passes on normal confidence", 0) \
      ,(bool,           featureConfMin,                 false,   "whether to assign median via the minimum-operator on the previous confidence values", 0) \
      ,(unsigned int,   vis3DPointSize,                 2,       "default point size in 3D visualization", 0) \
      ,(unsigned int,   visTrackHistory,                5,       "number of past track states to render at most", 0) \
      ,(double,         visMaxTrackVelVarPos,           0.0,     "maximum positional variance allowed (m/s)^2, currently unused???", 0) \
      ,(double,         visMaxTrackVelVarAng,           0.0,     "maximum angular variance allowed (rad/s)^2, if >0 renders probabilistic-trace, =0 icp-trace", 0) \
      ) )
#define PARAMHEAP_BASE_ARR_2 (23, ( \
       (bool,           segCrit1Active,                 true,    "", 0) \
      ,(bool,           segCrit2Active,                 true,    "", 0) \
      ,(bool,           segCrit3Active,                 true,    "", 0) \
      ,(unsigned int,   segVariantTilt,                 1,       "", 0) \
      ,(unsigned int,   segVariantTwist,                1,       "", 0) \
      ,(double,         segConvThresh,                  6.0,     "Convexity Threshold in degrees [0..90]", 0) \
      ,(double,         segConvFac,                     1.0,     "Convexity factor = threshold hardness", 0) \
      ,(double,         segTwistThresh,                 20.0,    "Twisting Threshold in degrees [0..90]", 0) \
      ,(double,         segTwistFac,                    5.0,     "Twisting factor = threshold hardness", 0) \
      ,(double,         segNormThresh,                  18.0,    "Threshold for normal similarity in degrees [0..90] = degrees at which score=0.5", 0) \
      ,(double,         segMinScore,                    0.4,     "Minimum score needed to connect pixels [0..1]", 0) \
      ,(unsigned int,   segMinSizePt,                   10,      "Minimum number of pixels needed for a segment to be accepted", 0) \
      ,(double,         segMinSizeM,                    0.2,     "Minimum side-length in meters for a segment to be accepted", 0) \
      ,(bool,           segUseNormConf,                 false,   "Flag, whether to use normal confidences for segmentation", 0) \
      ,(bool,           segNcdSquare,                   false,   "whether to square normal confidence (when segUseNormConf==true)", 0) \
      ,(bool,           segUseMultiscale,               false,   "Flag, whether to evaluate connections of length 2", 0) \
      ,(bool,           segUseTriangles,                true,    "Flag, whether to evaluate triangular connections for segmentation", 0) \
      ,(bool,           segUseTrackIDs,                 false,   "Flag, whether to bias segmentation by using projected tracks", 0) \
      ,(double,         segLinSatVelWeightThresh,       1.0,     "Threshold for weight2 when trackIDs are used", 0) \
      ,(double,         projectionDistTolerance,        0.2,     "Threshold for assigning projected points to pixels, in meter", 0) \
      ,(double,         matchWeightFeat,                2.0,     "Set this > 1 to weight features more than Cartesian distance", 0) \
      ,(double,         matchWeightZ,                   5.0,     "Set this > 1 to weight z-coordinate, thus preferring neighbors with similar height, also used for nearest-neighbor-search", 0) \
      ,(double,         matchMaxFSpaceDist,             10.0,    "Maximum distance in feature space to accept match (in meter in 7dim feature space)", 0) \
      ) )
#define PARAMHEAP_BASE_ARR_3 (22, ( \
       (double,         regCorrespDistWithHalfWeight,   5.0,     "weight = 0.5 at this distance", 0) \
      ,(unsigned int,   regMinSamplingCnt,              500,     "Minimum point count used for registration, logarithmic downsampling is used for tracks with more points", 0) \
      ,(bool,           regUseFeatMatch,                false,   "If true FeatureMatching is used when initial distance is too high", 0) \
      ,(double,         regFMatchAvgDistThresh,         0.05,    "// in meter", 0) \
      ,(double,         regMaxRelTransThresh,           0.25,    "Maximum distance for registration = (de)acceleration[m/s^2] * delta_t[s]^2.  E.g. 1g * (0.1s)^2 = 9.81*0.01 = 0.09", 0) \
      ,(unsigned int,   reg4DPtCountThresh,             10000,   "Point count of track from which full 6D estimate is performed, otherwise 4D only", 0) \
      ,(bool,           regIcpUseNormalWeight,          false,   "Flag, whether to use DistDiff image to weight ICP correspondences", 0) \
      ,(bool,           regIcpUseDistDiffWeight,        false,   "Flag, whether to use DistDiff image to weight ICP correspondences", 0) \
      ,(bool,           regIcpUseOcclWeight,            true,    "Flag, whether to give occluded correspondences lower weight in ICP algorithm", 0) \
      ,(double,         regIcpOcclWeight,               0.0001,  "If regIcpUseOcclWeight is set, correspondences are weighted by this parameter", 0) \
      ,(double,         regIcpVIPWeight,                2.0,     "Initial weight for VIP points", 0) \
      ,(double,         regIcpRegularizationCoeff,      0.000001,"coefficient used to regularize matrix", 0) \
      ,(unsigned int,   regIcpNbCorresp,                3,       "number of correspondences used for each model point", 0) \
      ,(double,         regMaxRelStateDiff,             2.0,     "masimum mahalanobis distance of measurement to predicted state so that measurement is accepted, 68% are within 1sigma, 90% are within 1.64sigma, 95% within 2sigma, 99.7% within 3sigma", 0) \
      ,(bool,           deepCopyTracks,                 true,    "Flag, whether Track objects are copied completely (good for debugging purposes) or whether the old track objects will be used and modified", 0) \
      ,(double,         segTrackLinkThresh,             0.1,     "Minimum relative pixel-ratio to create link between track and segment [0..1]", 0) \
      ,(double,         trkDelVarThresh,                2*2,      "Maximum allowed variance in m^2, otherwise track is deleted", 0) \
      ,(double,         trkDelUpCntThresh,              30,      "Max number without updates before track is deleted (seconds*10 for velodyne)", 0) \
      ,(double,         trkPGridWidth,                  0.1,     "cartesian grid-size for accumulated point cloud (meter)", 0) \
      ,(double,         trkNGridWidth,                  (M_PI/3),  "angular grid-size for accumulated point cloud (RAD)", 0) \
      ,(unsigned int,   trkPGridSearchRad,              3,       "search radius as number of cells", 0) \
      ,(unsigned int,   trkFeatureVersion,              2,       "1=50D/30D, 2=52D/32D", 0) \
      ) )
#define PARAMHEAP_BASE_ARR_4 (8, ( \
       (std::string,    stage1Classifier,               "stage1_classifier", "file with stage-1-classifier", 0) \
      ,(std::string,    stage1Normalizer,               "stage1_hyperplane.txt", "file with stage-1 normalization information", 0) \
      ,(std::string,    stage2Hyperplane,               "stage2_hyperplane.txt", "file with stage-2 hyperplane incl normalization information", 0) \
      ,(bool,           useDynamicMapping,              true,    "activate accumulation of track point cloud", 0) \
      ,(bool,           useSimpleMapping,               false,   "in combination with useDynamicMapping: if true point cloud is just accumulated, otherwise density checks are made to avoid infinite increase in size", 0) \
      ,(bool,           useOverwriteMapping,            false,   "in combination with useDynamicMapping: if true point cloud is replaced by new segment-point-cloud", 0) \
      ,(bool,           disableLocalization,            false,   "set true if the sensor always remains at the same position", 0) \
      ,(bool,           disableTracking,                false,   "set true to only estimate sensor motion", 0) \
      ) )
      // regIcpVIPWeight influences together with regIcpOcclWeight and regIcpNbCorresp the decision on when ICP is successful (sum_weight > 1.0), hence, 1/(regIcpNbCorresp*regIcpVIPWeight*regIcpOcclWeight)=1667 invisible VIP points are sufficient for successful registration
CREATE_PARAM_STRUCT(ParamHeapBase1, "Parameter Heap part 1", PARAMHEAP_BASE_ARR_1)
CREATE_PARAM_STRUCT(ParamHeapBase2, "Parameter Heap part 2", PARAMHEAP_BASE_ARR_2)
CREATE_PARAM_STRUCT(ParamHeapBase3, "Parameter Heap part 3", PARAMHEAP_BASE_ARR_3)
CREATE_PARAM_STRUCT(ParamHeapBase4, "Parameter Heap part 4", PARAMHEAP_BASE_ARR_4)
//CREATE_PARAM_STRUCT_STREAMOP(ParamHeapBase, PARAMHEAP_BASE_ARR, " ", 0) // use only in cpp!!




/*!
  \class ParameterHeap
  \brief Class for storing all parameters

  To avoid passing global parameter through routines, this
  class provides a central place for storing them.
  The implementation as Singleton thereby insures the
  existence of only once global instance (obtainable solely
  by the function get())
*/
class ParameterHeap : boost::noncopyable, public ParamHeapBase1, public ParamHeapBase2, public ParamHeapBase3, public ParamHeapBase4 {
public:
  static ParameterHeap* get();

  // parameters, made public for easy read/write:
  const double    defaultDeltaT; //!< default time-delta between frames
  const double    lidarDistStdDev; //!< Uncertainty of distance measurements
  const double    lidarHAngStdDevRAD; //!< Uncertainty of horizontal angle of measurements in radian
  const double    lidarVAngStdDevRAD; //!< Uncertainty of vertical angle of measurements in radian
  const double     lidarMeasStdDevConst; //!< Assumed standard deviation in x-y-z of lidar scanner on whole range, in meter
  const double     lidarMeasStdDevPerMeter; //!< Assumed additional standard deviation in x-y-z of lidar scanner measurements per distance-meter, in meter
  LidarImageFeatures::Neighborhood featureConfNeighb; //!< neighborhood-type of median
  double          regIcpMinErrOK; //!< error under which ICP result is always accepted", 0)
  double          predictNoise;
  double          trkDelMaxDist;  //!< Maximum distance from vehicle in meter", 0)

  mdefs::DSMatrix trkInitVariance; //!< The initial covariance-matrix for each new track.
  mdefs::DSMatrix trkPredVariance; //!< The noise-matrix added at each prediction step to covariance matrix

  template <class vector_t>
  void setSelection1(const vector_t &v);

  boost::program_options::options_description fullParameterDescription;

private:
  ParameterHeap();
  ~ParameterHeap();
//  static boost::mutex instanceMutex; // mutable -> changeable even if const object
  static ParameterHeap* uniqueInstance; // volatile -> might be changed by other thread
};

///////////////////////////////////////////////////////////////////////////////////
///////           template implementation           ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

template <class vector_t>
void ParameterHeap::setSelection1(const vector_t &pvec)
{
  if (pvec.size() != 25) {
    std::cerr << std::endl << "ParameterHeap::setSelection1: vector has size " << pvec.size() << std::flush;
    return;
  }
  connMaxDst                  = pvec[0]; // connMaxDst(0.2)
  connMaxDstRatio             = pvec[1]; // connMaxDstRatio(2.0)
  connFacA                    = pvec[2]; // connFacA(1.0)
  connFacB                    = pvec[3]; // connFacB(1.0)
  connFacC                    = pvec[4]; // connFacC(0.02)
  segConvThresh               = pvec[5]; // segConvThresh(6.0) DEG
  segConvFac                  = pvec[6]; // segConvFac(1.0)
  segTwistThresh              = pvec[7]; // segTwistThresh(20.0) DEG
  segTwistFac                 = pvec[8]; // segTwistFac(5.0)
  segNormThresh               = pvec[9]; // segNormThresh(10.0) DEG at which score=0.5
  segMinScore                 = pvec[10]; // segMinScore(0.45)
  segMinSizePt                = (unsigned int)(0.5+std::max(3.0,pvec[11])); // segMinSize(10)
  connFacD                    = pvec[12]; // connFacC(0.02)
  featurenNbConfMedianPasses  = (unsigned int)(0.5+std::max(0.0,pvec[13])); // 1
  segCrit1Active              = (pvec[16] > 0.5); // true
  segCrit2Active              = (pvec[17] > 0.5); // true
  segCrit3Active              = (pvec[18] > 0.5); // true
  segVariantTilt              = (unsigned int)(0.5+pvec[19]); // 1,2,3,[deact]
  segVariantTwist             = (unsigned int)(0.5+pvec[20]); // 1,2,[deact]
  switch ((int)(pvec[21]+0.5)) { // [4],8,24
    case 8: featureConfNeighb = LidarImageFeatures::NBH8; break;
    case 24: featureConfNeighb = LidarImageFeatures::NBH24; break;
    default: featureConfNeighb = LidarImageFeatures::NBH4;
  }
  featureConfMin              = (pvec[22] > 0.5); // false
  segUseMultiscale            = (pvec[14] > 0.5); // false
  segUseTriangles             = (pvec[24] > 0.5); // true
  segUseNormConf              = (pvec[15] > 0.5); // false
  segNcdSquare                = (pvec[23] > 0.5); // false
//  featureConnectWeightCreate = 2.0; // DONt optimize, should be quite ok
//  segUseTrackIDs - DONt change bool value
//  segLinSatVelWeightThresh - DONt change, only used if segUseTrackIDs == true
}

#endif /* PARAMETERHEAP_H_ */
