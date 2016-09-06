/*
 * LidarImageProjector.cpp
 *
 *  Created on: Jun 15, 2010
 *      Author: moosmann
 */

#include "LidarImageProjector.hpp"

//bool LidarImageProjector::getImageIndex(double xa, double ya, double za, int &hi, int &vi, double &dist) const
//{
//  boost::numeric::ublas::vector<double> pa(4); pa(0) = xa; pa(1) = ya; pa(2) = za; pa(3) = 1.0;
//  boost::numeric::ublas::vector<double> pr = prod(positionHTM2v, pa);
//  return getImageIndexRel(pr(0), pr(1), pr(2), hi, vi, dist);
//}

bool LidarImageProjector::getImageIndexRel(double xr, double yr, double zr, int &hi, int &vi, double &dist) const
{
  double xrs = xr*xr;
  double yrs = yr*yr;
  double zrs = zr*zr;
  //cout << boost::format("x/y/z = %f / %f / %f --> ") % xrs % yrs % zrs;
  double yawRAD = atan2(yr,xr); // atan2 returns in [-pi..+pi]
  double pitchRAD = atan2(-zr,sqrt(xrs+yrs)); // returns in [-pi..+pi]

  dist = sqrt(xrs+yrs+zrs);
  return getImageIndexYP(yawRAD, pitchRAD, hi, vi);
}

//void LidarImageProjector::get3DCoord(int hi, int vi, double dist, double &xa, double &ya, double &za) const
//{
//  ublas::vector<double> pr(4);
//  get3DCoordRel(hi, vi, dist, pr(1), pr(2), pr(3));
//  pr(3) = 1;
//  ublas::vector<double> pa = prod(positionHTM2w, pr);
//  xa = pa(0);
//  ya = pa(1);
//  za = pa(2);
//}

void LidarImageProjector::get3DCoordRel(int hi, int vi, double dist, matrixTools::DVector &pt) const
{
  double x, y, z;
  get3DCoordRel(hi, vi, dist, x, y, z);
  pt.resize(3);
  pt(0) = x;
  pt(1) = y;
  pt(2) = z;
}

