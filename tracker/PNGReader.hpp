#ifndef PNGREADER_H_
#define PNGREADER_H_

#include <LidarFrame.hpp>
#include <MatrixDefs.hpp>

class PNGReader
{

public:
  PNGReader(std::string directory);
  virtual ~PNGReader();

  /*!
   * \brief method reads specified files into obstacle tracker
   *
   * distance values are in meter, intensity from 0.0..1.0,
   * invalid distance values are set to DBL_MAX, invalid intensity to 0.0
   *
   * \param idx          specifies file to read, must be in range [0, getImageCount()-1]
   * \param incBufferIdx if true, data is read into the next frame and index is increased, otherwise the current frame is overwritten
   *
   * \throws range_error if vertical or horizontal size of LidarFrame images != size of LidarPreciseImager polar map
   */
  void        readImage(size_t idx, FrameSPtr currentF, FrameSPtr lastF);
  std::string getDstImageName(size_t idx) const {return frameSources[idx].dstImg;};
  std::string getImageConfigName() const {return imgCfgFile;};
  size_t      getImageCount() const {return frameSources.size();};

private:
  struct FrameSource {
    std::string dstImg; // distance-png-filenames with full path
    std::string intImg; // intensity-png-filenames with full path
    int64_t frameStartTime;
    matrixTools::DMatrix positionHTM2w;
    FrameSource() {
      dstImg = "";
      intImg = "";
      frameStartTime = 0;
      positionHTM2w = matrixTools::DIdMatrix(4);
    };
  };
  std::vector<FrameSource> frameSources;
  std::string imgCfgFile;

  void determineFilesExtractV1(std::string directory); // png-files + imu.cfg
  void determineFilesExtractV2(std::string directory); // velodyne_dist + velodyne_intens + oxts subdirectories
};

#endif // PNGREADER_H_

