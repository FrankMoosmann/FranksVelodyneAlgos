#ifndef PNGDISTANCEIMAGE_H_
#define PNGDISTANCEIMAGE_H_

#include <png++/png.hpp>
#include <limits>

/*!
 * \class PngDistanceImage
 * \brief A class for storing distance images as PNG files
 *
 * This class provides an easy way for storing/loading a 2D-array with distance values (in meter) as PNG image
 *
 * \author Frank Moosmann <moosmann@mrt.uka.de>
 **/
class PngDistanceImage
{
private:
  png::image<png::gray_pixel_16> image; //!< image that holds all the data and that implements save/load as PNG file
  static const png::gray_pixel_16 maxPix = 65535; // maximum coded distance value, should always be max of 16bit, thus equal to std::numeric_limits<png::gray_pixel_16>::max();
  static const double default_multiplier = 500.0; //!< multiplier to turn floating point distance values into discrete pixel values. chosen in a way that all occuring distances can be coded but resolution is still as high as possible, 1/500 = 2mm resolution
  double multiplier; //!< multiplier to turn floating point distance values into discrete pixel values. chosen in a way that all occuring distances can be coded but resolution is still as high as possible, 1/500 = 2mm resolution
  double maxDst; //!< maximum distance possible to store, equal to ((double)maxPix)/multiplier
  png::gray_pixel_16 dist2PixVal(double distance) { if ((distance > maxDst) || (distance < 0.0)) distance = 0.0; return png::gray_pixel_16(distance*multiplier);};

public:
  PngDistanceImage(std::string const &pngfile, double multiplier = default_multiplier) { setMultiplier(multiplier); load(pngfile);};
  PngDistanceImage(unsigned int width, unsigned int height, double multiplier = default_multiplier) : image(width, height) {setMultiplier(multiplier);};
  ~PngDistanceImage() {};

  static double getDefaultMultiplier() {return default_multiplier;}; //! returns the multiplier to turn PNG-file values into meters
  double getMultiplier() {return multiplier;}; //! returns the multiplier to turn PNG-file values into meters
  void setMultiplier(double mult) {multiplier=mult; maxDst=maxPix/multiplier;}; //! set the multiplier to turn PNG-file values into meters
  double maxDist() const {return maxDst;}; //!< returns the maximum distance in meters that can be encoded
  unsigned int width() const {return image.get_width();}; //!< returns the width of the image
  unsigned int height() const {return image.get_height();}; //!< returns the height of the image
  void load(std::string const &pngfile) //!< loads image from file -> width&height might change
  {
    image.read(pngfile);
    for (unsigned int row=0; row<image.get_height(); ++row) {
      for (unsigned int col=0; col<image.get_width(); ++col) {
        if (image[row][col] == 0)
          image[row][col] = maxPix;
      }
    }
  };
  void save(std::string const &pngfile) {image.write(pngfile);}; // can't declare const, because image.write is not a const member
  void setDistance(unsigned int row, unsigned int col, double distance) //!< set distance in meter
  {
    image[row][col] = dist2PixVal(distance);
  };
  void fill(double distance) //!< fills the whole image with the given distance value
  {
    png::gray_pixel_16 pixVal = dist2PixVal(distance);
    for (unsigned int row=0; row<image.get_height(); ++row) {
      for (unsigned int col=0; col<image.get_width(); ++col) {
        image[row][col] = pixVal;
      }
    }
  };
  template <typename dist_val_t> // made template to enabled both float and doulbe
  void setDistances(dist_val_t const *data_p, unsigned int alignment = 0) //!< set all image values from a data block of (minimum) size (width+alignment)*height
  {
    for (unsigned int row=0; row<image.get_height(); ++row) {
      for (unsigned int col=0; col<image.get_width(); ++col) {
        setDistance(row,col,*data_p);
        ++data_p;
      }
      data_p += alignment;
    }
  };
  template <typename dist_val_t> // made template to enabled both float and doulbe
  dist_val_t getDistance(unsigned int row, unsigned int col) const //!< get distance in meter
  {
    if (image[row][col] == maxPix)
      return std::numeric_limits<dist_val_t>::max(); // return maximum value of dist_val_t
    return (dist_val_t)(image[row][col])/(dist_val_t)multiplier;
  };
  template <typename dist_val_t> // made template to enabled both float and doulbe
  void getDistances(dist_val_t *data_p, unsigned int alignment = 0) const //!< read all image values into a data block of (minimum) size (width+alignment)*height
  {
    for (unsigned int row=0; row<image.get_height(); ++row) {
      for (unsigned int col=0; col<image.get_width(); ++col) {
        *data_p = getDistance<dist_val_t>(row,col);
        ++data_p;
      }
      data_p += alignment;
    }
  };
};


#endif /*PNGDISTANCEIMAGE_H_*/
