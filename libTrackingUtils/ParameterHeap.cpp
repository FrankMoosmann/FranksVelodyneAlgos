/*
 * ParameterHeap.cpp
 *
 *  Created on: Nov 6, 2009
 *      Author: moosmann
 */

#include "ParameterHeap.hpp"

ParameterHeap* ParameterHeap::uniqueInstance = 0;
//boost::mutex instanceMutex = boost::mutex();

using namespace matrixTools;
using namespace std;

ParameterHeap::ParameterHeap()
  : defaultDeltaT(0.1) // 0.1seconds
   ,lidarDistStdDev(0.015) // 1.5cm on whole range
// change the following in accordance to pixel quantization:
   ,lidarHAngStdDevRAD(0.0015) // 0.09° due to angle sensor
   ,lidarVAngStdDevRAD(0.00017) // 0.01° due to calibration error (angle is fixed)
   ,lidarMeasStdDevConst(lidarDistStdDev) // constant on whole range plus
   ,lidarMeasStdDevPerMeter(0.0016) // distance-dependent value (from angle uncertainty)
   ,featureConfNeighb(LidarImageFeatures::NBH4)
   ,regIcpMinErrOK(3*lidarDistStdDev)
   ,predictNoise(regMaxRelTransThresh)
   ,trkDelMaxDist(framegenMaxDistance*0.95) // in meter
   ,trkInitVariance(ublas::zero_matrix<double>(12,12)) // meter*meter
   ,trkPredVariance(ublas::zero_matrix<double>(12,12)) // meter*meter
{
  fullParameterDescription.add(ParamHeapBase1::parameterDescription);
  fullParameterDescription.add(ParamHeapBase2::parameterDescription);
  fullParameterDescription.add(ParamHeapBase3::parameterDescription);
  fullParameterDescription.add(ParamHeapBase4::parameterDescription);
  // trkInitVariance(0,0)..(5,5): no position uncertainty, as coordinate system is chosen arbitrarily
  trkInitVariance(6,6) = M_PI*M_PI*0.25; // expected deviation of 90°/s = 0.5PI/s => variance = 0.25*PI^2
  trkInitVariance(7,7) = M_PI*M_PI*0.25;
  trkInitVariance(8,8) = M_PI*M_PI*0.25;
  trkInitVariance(9,9) = 100.0; // expected deviation of 36 km/h = 10m/s => variance = 100
  trkInitVariance(10,10) = 100.0;
  trkInitVariance(11,11) = 20.0;
  double val = 0.01; // angular error = 1/2*a*t^2?
  trkPredVariance(0,0) = val;
  trkPredVariance(1,1) = val;
  trkPredVariance(2,2) = val;
  double aPosMax = 3.0; // expected maximum positional acceleration [m/s^2]
  val = 0.5 * aPosMax * defaultDeltaT * defaultDeltaT; // position-error due to acceleration = delta_p = 0.5 * a_max[m/s^2] * delta_t[s]^2
  trkPredVariance(3,3) = val;
  trkPredVariance(4,4) = val;
  trkPredVariance(5,5) = val;
  val = aPosMax / 20.0; // angular acceleration (circle) =  tangential acceleration[m/s^2] / radius[m]
  trkPredVariance(6,6) = val;
  trkPredVariance(7,7) = val;
  trkPredVariance(8,8) = val;
  val = aPosMax * defaultDeltaT; // velocity-error due to acceleration = delta_v = a_max[m/s^2] * delta_t[s]
  trkPredVariance(9,9) = val;
  trkPredVariance(10,10) = val;
  trkPredVariance(11,11) = val;
}

ParameterHeap::~ParameterHeap() {
}

ParameterHeap* ParameterHeap::get() {
//  if (uniqueInstance == 0) {
//    boost::mutex::scoped_lock lock(instanceMutex);
    if (uniqueInstance == 0) { // double check in case multithreaded apps call it in parallel
      uniqueInstance = new ParameterHeap();
    }
//  }
  return uniqueInstance;
}
