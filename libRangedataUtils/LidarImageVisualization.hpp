#ifndef LIDARIMAGEVISUALIZATION_H_
#define LIDARIMAGEVISUALIZATION_H_

#include <QImage>
#include <QColor>
#include <GL/glut.h>

#include "LidarImage.hpp"
#include <MatrixDefs.hpp>

namespace LidarImageVisualization
{
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////       Helper Methods      /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  inline void HSV2RGB(int H, int S, int V, int &r, int &g, int &b) { // H = 0..180; S,V=0..255
    QColor c = QColor::fromHsv(H, S, V);
    c.getRgb(&r, &g, &b);
  }

  inline void glColorHSV(int H, int S, int V) { // H = 0..360; S,V=0..255
    int r,g,b; HSV2RGB(H, S, V, r,g,b);
    glColor3f(r/255.0f,g/255.0f,b/255.0f); // color according to "speed"
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////    2D Rendering Methods   /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  void renderMono(const LidarImage<bool> &li, QImage &qi, bool invert = false);
  enum ColorMode2D {HUE, HUEP3, REDGREEN, GREENRED};
  void renderColored(const LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs = false, double valInvalid = DBL_MAX, ColorMode2D mode = HUEP3);
  void renderGrey(const LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs = false, bool invert = false);
  void renderColoredNormals(const LidarImage<matrixTools::DVector> &lin, QImage &qi, const LidarImage<double> *conf = NULL);

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////    3D Rendering Methods   /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  enum ColorMode3D {ColorFixed, ColorHeight, ColorDist};
  void render3DPoints(const LidarImage<matrixTools::DVector> &points, ColorMode3D mode=ColorFixed, unsigned int pointSize=2, double r=0.9, double g=0.9, double b=0.9, const LidarImage<double> *distance = NULL, double minDist = 2.0, double maxDist = 80.0);
  void render3DPointCovar(const LidarImage<matrixTools::DVector> &points, const LidarImage<matrixTools::DMatrix> &pointCovar);
  void render3DNormals(const LidarImage<matrixTools::DVector> &points, const LidarImage<matrixTools::DVector> &normals, const LidarImage<double> *normalConfidence = NULL, const LidarImage<double> *normalStdDevRAD = NULL, const LidarImage<matrixTools::DMatrix> *normalCovar = NULL, unsigned int lineWidth = 2, double lengthFactor = 0.1);
  void render3DConnections(const LidarImage<double> &val, unsigned int pixOff1stHoriz, unsigned int pixOff1stVert, unsigned int pixOff2ndHoriz, unsigned int pixOff2ndVert, const LidarImage<matrixTools::DVector> &points);

}

#endif /*LIDARIMAGEVISUALIZATION_H_*/
