// Author:         Andreas Geiger <geiger@kit.edu>

#if !defined(CONVERT_COORDINATES_HPP)
#define CONVERT_COORDINATES_HPP

#include <cmath>

/*!
  * \file convert_coordinates.hpp
  *
  * \brief provides functions to convert global lat/long into local cartesian x/y coordinates
  *
  * the following functions map lat/long coordinates to the euclidean mercator coordinate system mx/my
  * lat/long are always specified in degrees (DEG)
  * mx is the distance from the central meridian towards EAST
  * my is the distance from the equator towards NORTH
  * mercator cordinates mx/my correspond to metric coordinates x/y ONLY at the equator (scale=1)
  * at other latitudes a scale factor can be used to get metric relations
  * for more details see GCDC/docs/coordinate_systems.pdf
  * \note that in GCDC x is towards NORTH, y towards WEST!!!
  */

namespace convert_coordinates {

  const double EARTH_RADIUS_EQUA = 6378137.0; // earth radius at equator [m]
  const int _tileSize = 256;
  const double _initialResolution = 2.0 * M_PI * EARTH_RADIUS_EQUA / (double)_tileSize; // pixel-size in meter at level 0
  const int _originShift = M_PI * EARTH_RADIUS_EQUA; // tiles-pixels start at lower left corner but level 0 tile is centered at lat/lon=0/0

  inline double levelResolution(const int &levelOfDetail) { // pixel-size in meter for given level
      return _initialResolution / (double)(pow((int)2,levelOfDetail));
  }

  // inlined functions to avoid multiple definitions

  /*! \brief convert latitude to scale, which is needed by mercator transformations
   *  \param lat latitude in degrees (DEG)
   *  \return scale factor
   *  \note when converting from lat/lon -> mercator and back again,
   *        or vice versa, use the same scale in both transformations!
   */
  inline double lat_to_scale (double lat) {
    return cos(lat * M_PI / 180.0);
  }

  /*! \brief converts lat/lon/scale to mx/my (mx/my in meters if correct scale is given)
   */
  inline void latlon_to_mercator (double lat, double lon, double scale, double &mx, double &my) {
    mx = scale * lon * M_PI * EARTH_RADIUS_EQUA / 180.0;
    my = scale * EARTH_RADIUS_EQUA * log( tan((90.0+lat) * M_PI / 360.0) );
  }

  /*! \brief convenience function, uses lat0 to calculate appropriate scale
   */
  inline void latlon_to_scaled_mercator (double lat, double lon, double lat0, double &mx, double &my) {
    double scale = lat_to_scale( lat0 );
    mx = scale * lon * M_PI * EARTH_RADIUS_EQUA / 180.0;
    my = scale * EARTH_RADIUS_EQUA * log( tan((90.0+lat) * M_PI / 360.0) );
  }

  /*! \brief converts mx/my/scale to lat/lon (mx/my in meters if correct scale is given)
   */
  inline void mercator_to_latlon (double mx, double my, double scale, double &lat, double &lon) {
    lon = mx * 180.0 / (M_PI * EARTH_RADIUS_EQUA * scale);
    lat = 360.0 * atan( exp(my/(EARTH_RADIUS_EQUA * scale)) ) / M_PI - 90.0;
  }
  
  /*! \brief convenience function, uses lat0 to calculate appropriate scale
   */
  inline void scaled_mercator_to_latlon (double mx, double my, double lat0, double &lat, double &lon) {
    double scale = lat_to_scale( lat0 );
    lon = mx * 180.0 / (M_PI * EARTH_RADIUS_EQUA * scale);
    lat = 360.0 * atan( exp(my/(EARTH_RADIUS_EQUA * scale)) ) / M_PI - 90.0;
  }

  /*! \brief adds meters dx/dy to given lat/lon and returns new lat/lon
   */
  inline void latlon_add_meters (double lat_start, double lon_start, double dx, double dy, double &lat_end, double &lon_end) {
    double scale = lat_to_scale (lat_start);
    double mx,my;
    latlon_to_mercator (lat_start, lon_start, scale, mx, my);
    mx += dx;
    my += dy;
    mercator_to_latlon (mx, my, scale, lat_end, lon_end);
  }

  /*! \brief given two lat/lon coordinates, returns their difference in meters dx/dy
   */
  inline void latlon_diff_to_meters (double lat_start, double lon_start, double lat_end, double lon_end, double &dx, double &dy) {
    double scale = lat_to_scale (lat_start);
    double mx1,my1, mx2, my2;
    latlon_to_mercator (lat_start, lon_start, scale, mx1, my1);
    latlon_to_mercator (lat_end, lon_end, scale, mx2, my2);
    dx = mx2-mx1;
    dy = my2-my1;
  }

  inline void LatLonToMeters (const double &lat, const double &lon, double &mx, double& my) {
      mx = lon * M_PI * EARTH_RADIUS_EQUA / 180.0;
      my = log( tan((90.0 + lat) * M_PI / 360.0 )) / (M_PI / 180.0);
      my = my * M_PI * EARTH_RADIUS_EQUA / 180.0;
  }

  inline void MetersToLatLon (const double &mx, const double &my, double &lat, double& lon) {
      lon = (mx / M_PI / EARTH_RADIUS_EQUA) * 180.0;
      lat = (my / M_PI / EARTH_RADIUS_EQUA) * 180.0;
      lat = 180.0 / M_PI * (2.0 * atan( exp( lat * M_PI / 180.0)) - M_PI / 2.0);
  }

  inline void MetersToPixels (const double &mx, const double &my, const int &levelOfDetail, double &px, double& py) {
      double res = levelResolution(levelOfDetail);
      px = (mx + _originShift) / res;
      py = (my + _originShift) / res;
  }

  inline void PixelsToMeters (const double &px, const double &py, const int &levelOfDetail, double &mx, double& my) {
      double res = levelResolution(levelOfDetail);
      mx = px * res - _originShift;
      my = py * res - _originShift;
  }

  inline void PixelsToMetersKarlsruhe (const double &px, const double &py, const int &levelOfDetail, double &mx, double& my) {
      double res = levelResolution(levelOfDetail);
      mx = (px * res - _originShift)*0.65;
      my = (py * res - _originShift)*0.65;
  }

  inline void LatLonToPixels (const double &lat, const double &lon, const int &levelOfDetail, double &px, double& py) {
      double mx,my;
      LatLonToMeters(lat,lon,mx,my);
      MetersToPixels(mx,my,levelOfDetail,px,py);
  }

  inline void PixelsToLatLon (const double &px, const double &py, const int &levelOfDetail, double &lat, double& lon) {
      double mx,my;
      PixelsToMeters(px,py,levelOfDetail,mx,my);
      MetersToLatLon(mx,my,lat,lon);
  }

  inline void PixelsToTile (const double &px, const double &py, int &tx, int& ty) {
      tx = (int)(ceil(px/(double)(_tileSize)) - 1.0);
      ty = (int)(ceil(py/(double)(_tileSize)) - 1.0);
  }

  inline void MetersToTile (const double &mx, const double &my, const int &levelOfDetail, int &tx, int& ty) {
      double px,py;
      MetersToPixels(mx,my,levelOfDetail,px,py);
      PixelsToTile(px,py,tx,ty);
  }

  inline void TileToMeters (const int &tx, const int& ty, const int &levelOfDetail, double &px_min, double &px_max, double &py_min, double &py_max) {
      PixelsToMeters((double)((tx+0)*_tileSize),(double)((ty+0)*_tileSize),levelOfDetail,px_min,py_min);
      PixelsToMeters((double)((tx+1)*_tileSize),(double)((ty+1)*_tileSize),levelOfDetail,px_max,py_max);
  }

  inline void LatLonToTile (const double &lat, const double &lon, const int &levelOfDetail, int &tx, int& ty) {
      double mx,my;
      LatLonToMeters(lat,lon,mx,my);
      MetersToTile(mx,my,levelOfDetail,tx,ty);
  }

  inline void TileToLatLon (const int &tx, const int& ty, const int &levelOfDetail, double &lat_min, double &lat_max, double &lon_min, double &lon_max) {
      PixelsToLatLon((double)((tx+0)*_tileSize),(double)((ty+0)*_tileSize),levelOfDetail,lat_min,lon_min);
      PixelsToLatLon((double)((tx+1)*_tileSize),(double)((ty+1)*_tileSize),levelOfDetail,lat_max,lon_max);
  }
};

#endif

