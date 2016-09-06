#ifndef ALGOVISUALIZATION_H_
#define ALGOVISUALIZATION_H_

#include <map>
#include <boost/tuple/tuple.hpp>
#include <QImage>
#include <GL/glut.h>

#include "LidarImage.hpp"
#include "LidarFrame.hpp"
#include "LidarSegment.hpp"
#include "PrincipalComponentAnalysis.hpp"
#include "MatrixDefs.hpp"
#include "LidarImageVisualization.hpp"
#include "AlgoRegistration.hpp"


class AlgoVisualization
{
public:
  AlgoVisualization();
  virtual ~AlgoVisualization();
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////    2D Rendering Methods   /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
        void renderMono(LidarImage<bool> &li, QImage &qi, bool invert = false);
  void renderColored(LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs = false, double valInvalid = DBL_MAX, LidarImageVisualization::ColorMode2D mode = LidarImageVisualization::HUEP3);
  void renderGrey(LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs = false, bool invert = false);
  void renderColoredNormals(LidarImage<matrixTools::DVector> &lin, LidarImage<double> &conf, QImage &qi);
  void renderSegments(LidarImage<LidarSegment> &li, QImage &qi);
  void renderTrackProjection(const LidarImage<LidarFrame::TrackIndex> &trackProj, const LidarFrame::TrackSPtrVec &tracks, QImage &qi);
  void renderTrackProjectionCnt(LidarImage<unsigned int> &tids, QImage &qi);
  void renderTranslationNorm(FrameSPtr cf, QImage &qi);

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////    3D Rendering Methods   /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  void render3DCarPos(bool renderCS = false);
  void render3DGrid();
  void render3DPoints(FrameSPtr frame, bool renderUncertaint, bool invertColors);
  void render3DNormals(FrameSPtr frame, bool renderUncertainty, bool invertColors);
  void render3DConnections(LidarImage<matrixTools::DVector> &points, LidarImage<double> &connR, LidarImage<double> &connD, LidarImage<double> *connRU = NULL, LidarImage<double> *connRD = NULL);
  void render3DSegmentationLinks(FrameSPtr frame);
  void render3DSegments(FrameSPtr frame, bool invertColors);
  void render3DCorrespondences(FrameSPtr frame, FrameSPtr nframe, bool invertColors);
  void render3DRegistrationResult(TrackRegistration::SubTrackLst::iterator track); // track must be a valid iterator, not pointing to end()
  void render3DTrackProps(PointCloudTrack *track); // draws all tracks of a given trackList, selecting color by track properties
  void render3DTrack(PointCloudTrack *track, bool full = true, double R = -1, double G = -1, double B = -1, matrixTools::DMatrix *htm2Ego = NULL); // draw track. if R/G/B < 0 use tracks' color, if htm2Ego == NULL use track's state
  void render3DTrackVariance(PointCloudTrack *track); // track must be a valid track object of the given frame
  void render3DTrackHistory(PointCloudTrack *track); // track must be a valid track object of the given frame

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////        Helper Methods       ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  void render3DNormalPlane(const matrixTools::DVector &p, const matrixTools::DVector &n, double radius = 0.2, GLUquadricObj *quadric = NULL); // render normal vector as disk around point with current color. existing quadric is used if specified
};

#endif /*ALGOVISUALIZATION_H_*/
