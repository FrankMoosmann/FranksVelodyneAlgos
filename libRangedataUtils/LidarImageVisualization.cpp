#include "LidarImageVisualization.hpp"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <GL/glut.h>
#include <CovarEllipsoidRendering.hpp>

using namespace std;
using namespace matrixTools;
using namespace matrixTools::ublas;


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    2D Rendering Methods   /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void LidarImageVisualization::renderMono(const LidarImage<bool> &li, QImage &qi, bool invert)
{
  int hsize, vsize;
  li.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      unsigned int val = li.get(col,row) ? 255 : 0;
      if (invert) val = 255 - val;
      qi.setPixel( col, row, QColor(val, val, val).rgb() );
    }
  }
}

void LidarImageVisualization::renderColored(const LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs, double valInvalid, ColorMode2D mode)
{
  int hsize, vsize;
  li.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      double dat = li.get(col,row);
      double val = 0.0;
      double hue = 0.0;
      if (dat != valInvalid) { // only if pixel is valid
        if (useAbs) dat = fabs(dat);
        dat = max(minVal, min(maxVal,dat)); // cut off value, so it is within [minVal...maxVal]
        dat = 1.0 - (dat - minVal)/(maxVal-minVal); // scale to [0...1]
        switch (mode) {
          case HUE : {
            val = dat*255.0; // scale to 0...255 
            hue = 359.0 - dat*359.0; //scale to 0..359
          }; break;
          case HUEP3 : {
            val = dat*255.0; // scale to 0...255 
            hue = 359.0 - dat*dat*dat*359.0; //scale to 0..359
          }; break;
          case REDGREEN : {
            val = 255.0; 
            hue = dat*120.0; //scale to 0(red)..120(green)
          }; break;
          case GREENRED : {
            val = 255.0; 
            hue = 120.0 - dat*120.0; //scale to 0(red)..120(green)
          }; break;
        }
      }
      qi.setPixel( col, row, QColor::fromHsv(hue, 255, val).rgb() );
    }
  }
}

void LidarImageVisualization::renderGrey(const LidarImage<double> &li, QImage &qi, double minVal, double maxVal, bool useAbs, bool invert)
{
  int hsize, vsize;
  li.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16); //Format_Mono);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      double dat = li.get(col,row);
      int val = 0;
      if (dat < DBL_MAX) { // only if pixel is valid
        if (useAbs) dat = fabs(dat);
        dat = max(minVal, min(maxVal,dat)); // now value is within [minVal...maxVal]
        dat = (dat - minVal)/(maxVal-minVal); // now value is within [0...1]
        val = dat*dat*255.0; // scale to 255...0 
      }
      if (invert)
        val = 255.0 - val;
      qi.setPixel( col, row, QColor(val, val, val).rgb() );
    }
  }
}

void LidarImageVisualization::renderColoredNormals(const LidarImage<DVector> &normals, QImage &qi, const LidarImage<double> *conf)
{
  int hsize, vsize;
  normals.getSize(hsize,vsize);
  qi = QImage(hsize,vsize,QImage::Format_RGB16);
  for (int row=0; row<vsize; ++row) {
    for (int col=0; col<hsize; ++col) {
      DVector n = normals.get(col, row);
      double length = norm_2(n);
      if (length < 0.01) {
        qi.setPixel( col, row, QColor(0,0,0).rgb() );
      } else {
        double r = ((n(0)/length)+1.0)/2.0; // x-value in range -1..1 normalized to 0..1
        double g = ((n(1)/length)+1.0)/2.0; // y-value in range -1..1 normalized to 0..1
        double b = ((n(2)/length)+1.0)/2.0; // z-value in range -1..1 normalized to 0..1
        if (conf) {
          double c = conf->get(col, row);
          r*=c; g*=c; b*=c;
        }
//        cout << "x: " << nx << ",   y: " << ny << ",   z: " << nz;
//        cout << "r: " << r << ",   g: " << g << ",   b: " << b << ",  c: " << c << endl;
        qi.setPixel( col, row, QColor(r*254,g*254,b*254).rgb() );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    3D Rendering Methods   /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void LidarImageVisualization::render3DPoints(const LidarImage<DVector> &points, LidarImageVisualization::ColorMode3D mode, unsigned int pointSize, double r, double g, double b, const LidarImage<double> *distance, double minD, double maxD)
{
  if ((mode == ColorDist) && (!distance)) {
    cerr << "warning: LidarImageVisualization::render3DPoints called with mode==ColorDist but not image specified. switching to constant color." << endl;
    mode = ColorFixed;
  }
  bool invert = false;
  double DRange = maxD - minD;
  glPointSize(pointSize);
  glBegin(GL_POINTS);
  int hsize, vsize;
  points.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      const DVector &p = points.get(col, row);
      switch (mode) {
        case ColorFixed: {
          glColor3f(r, g, b);
        } break;
        case ColorHeight: {
          float height = p(2);
          if (height < 0) {
            float ratio = min(1.0,height/(-2.5));
            r = invert  ?  0.2        :  1.0-ratio;
            g = invert  ?  (1-ratio)  :  1.0-ratio;
            b = invert  ?  ratio      :  1.0;
          } else {
            float ratio = min(1.0,height/0.5);
            r = invert  ?  ratio      :  1.0;
            g = invert  ?  (1-ratio)  :  1.0-ratio; // normalize to range 0..1
            b = invert  ?  0.2        :  1.0-ratio;
          }
          glColor3f(r, g, b);
        } break;
        case ColorDist: {
          double drel = (min(max(distance->get(col, row),minD),maxD)-minD)/DRange;
          glColorHSV(fmod(sqrt(drel),0.5)*2*359, 255, 255); // H = 0..360; S,V=0..255
        } break;
      }

      if (p(0) != DBL_MAX)
        glVertex3f(p(0), p(1), p(2));
    }
  }
  glEnd();
}

void LidarImageVisualization::render3DPointCovar(const LidarImage<DVector> &points, const LidarImage<matrixTools::DMatrix> &pointCovar)
{
  int hsize, vsize;
  points.getSize(hsize,vsize);
  int divisor = max(1,hsize*vsize/3000); // maximum number of ellipsoids
  glColor3f(0.5, 0.5, 0.5);
  GLUquadricObj *quadric = gluNewQuadric();
  gluQuadricDrawStyle(quadric, GLU_LINE); // rendering style (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      if (col%divisor == 0) {
        const DVector &p = points.get(col, row);
        glPushMatrix();
        glTranslatef(p(0),p(1),p(2)); // move coordinate system to point
        render3DCovar33Ellipsoid(pointCovar.get(col, row), 1.0);
        glPopMatrix();
      }
    }
  }
  gluDeleteQuadric(quadric);
}

void LidarImageVisualization::render3DNormals(const LidarImage<matrixTools::DVector> &points, const LidarImage<matrixTools::DVector> &normals,
    const LidarImage<double> *normalConfidence, const LidarImage<double> *normalStdDevRAD, const LidarImage<matrixTools::DMatrix> *normalCovar,
    unsigned int lineWidth, double lengthFactor)
{
  const double alpha = 0.7;
  double r,g,b;
  int hsize, vsize;
  points.getSize(hsize,vsize);
  int divisor = max(1,hsize*vsize/3000); // maximum number of ellipsoids

  // render normals as lines
  glLineWidth(lineWidth);
  glBegin(GL_LINES);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      const DVector &p = points.get(col, row);
      const DVector &n = normals.get(col, row);
      double length = norm_2(n);
      if ((p(0) == DBL_MAX) || (length < 0.1))
        continue;

      if (normalConfidence) {
        double nconf = normalConfidence->get(col, row);
        r = 1.0-nconf;
        g = nconf;
        b = 0.0;
      } else {
        r = ((n(0)/length)+1.0)/2; // x-value normalized to -1..1
        g = ((n(1)/length)+1.0)/2; // y-value in range -1..1
        b = ((n(2)/length)+1.0)/2; // z-value in range -1..1
      }
      glColor4f(r, g, b, alpha);
      glVertex3f(p(0), p(1), p(2));
      glVertex3f(p(0)+n(0)*lengthFactor, p(1)+n(1)*lengthFactor, p(2)+n(2)*lengthFactor);
    }
  }
  glEnd();
  // render uncertainty
  GLUquadricObj *quadric = gluNewQuadric();
  gluQuadricDrawStyle(quadric, GLU_LINE); // rendering style (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      const DVector &p = points.get(col, row);
      const DVector &n = normals.get(col, row);
      double length = norm_2(n);
      if ((p(0) == DBL_MAX) || (length < 0.1))
        continue;

      if (normalConfidence) {
        double nconf = normalConfidence->get(col, row);
        r = 1.0-nconf;
        g = nconf;
        b = 0.0;
      } else {
        r = ((n(0)/length)+1.0)/2; // x-value normalized to -1..1
        g = ((n(1)/length)+1.0)/2; // y-value in range -1..1
        b = ((n(2)/length)+1.0)/2; // z-value in range -1..1
      }
      glColor4f(r, g, b, alpha);
      if (normalStdDevRAD) { // render cone
        if (col%divisor == 0) {
          double nSD = min(M_PI/4,normalStdDevRAD->get(col, row));
          glPushMatrix();
          glTranslatef(p(0),p(1),p(2)); // move coordinate system to point
          glRotatef(acos(n(2))*180/M_PI, -n(1), n(0), 0.0);// rotate so that normal direction is (0,0,1).
          gluCylinder(quadric, 0.0, tan(nSD)*lengthFactor, lengthFactor, 6, 1); // baseRadius, topRadius, height, slices, stacks
          glPopMatrix();
        }
      }
      if (normalCovar) { // render covar ellipsoids
        if (col%divisor == 0) {
          glPushMatrix();
          glTranslatef(p(0)+n(0)*lengthFactor,p(1)+n(1)*lengthFactor,p(2)+n(2)*lengthFactor); // move coordinate system to point
          render3DCovar33Ellipsoid(normalCovar->get(col, row), lengthFactor);
//          cout << normalCovar->get(col, row) << flush;
          glPopMatrix();
        }
      }
    }
  }
  gluDeleteQuadric(quadric);
}

void LidarImageVisualization::render3DConnections(const LidarImage<double> &val, unsigned int hOff1, unsigned int vOff1, unsigned int hOff2, unsigned int vOff2, const LidarImage<matrixTools::DVector> &points)
{
  glBegin(GL_LINES);
//  const double alpha = 0.5;
  const char alphaC = 128;
  int hsize, vsize;
  val.getSize(hsize,vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double c = val.get(col,row);
//      if (c > 0.0) {
        const DVector &p1 = points.get(col+hOff1, row+vOff1);
        const DVector &p2 = points.get(col+hOff2, row+vOff2);
        if ((p1(0) != DBL_MAX) && (p2(0) != DBL_MAX)) {
          c = fmax(0.0, fmin(1.0,c)); // cut off value, so it is within [minVal...maxVal]
          // same as ColorMode2D:REDGREEN
          double hue = c*120.0; //scale to 0(red)..120(green)
          double val = 255.0;
          // original method
//          double hue = c*120.0; //scale to 0..120
//          double val = 55.0 + c*200.0; // scale to 127...255
          int r,g,b,a;
          QColor hsv = QColor::fromHsv(hue, 255, val);
          hsv.getRgb(&r,&g,&b,&a);
          glColor4ub((char)r, (char)g, (char)b, alphaC);
          glVertex3f(p1(0), p1(1), p1(2));
          glVertex3f(p2(0), p2(1), p2(2));
        }
//      }
    }
  }
  glEnd();
}
