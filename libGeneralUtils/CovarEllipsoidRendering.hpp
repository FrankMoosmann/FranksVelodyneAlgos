#ifndef COVARELLIPSOIDRENDERING_H_
#define COVARELLIPSOIDRENDERING_H_

#include <GL/glut.h>
#include "MatrixDefs.hpp"
#include "PrincipalComponentAnalysis.hpp"

//! render 3x3 position covariance-ellipsoid around current center (0,0,0) [x,y,z]
template <class MatrixExpr>
void render3DCovar33Ellipsoid(MatrixExpr &covar, double stdDevFactor = 1.0, GLUquadricObj *quadric = NULL) {
  using namespace std;
  using namespace matrixTools;
  using namespace matrixTools::ublas;
  // assumption: covariance is a 3x3 matrix with x,y,z variances for translation
  BOOST_ASSERT(covar.size1() == 3 && "AlgoVisualization::render3DCovar33Ellipsoid: covar.size1 != 3");
  BOOST_ASSERT(covar.size2() == 3 && "AlgoVisualization::render3DCovar33Ellipsoid: covar.size2 != 3");

  // 1) determine principal axis of ellipsoid coordinate system
  DSMatrix sc = symmetric_adaptor<MatrixExpr, upper>(covar);
  PCA pca(sc);
  matrixTools::DVector eigenval(pca.getSortedEigenValues());
  matrixTools::DMatrix eigenvec(pca.getSortedEigenVectors());
  for (unsigned int i=0; i<3; ++i) {
    eigenval(i) = min(50.0,max(0.0001,stdDevFactor*sqrt(eigenval(i)))); // turn variance in std-dev and force length into interval 0.1mm..50m
  }
//  cout << "   EVec:" << eigenvec;
//  cout << "   EVal:" << eigenval;

  // 2) store axis in 4x4 matrix, column-major for OpenGL
  double axisarray[16];
  axisarray[0] = eigenvec(0,0); axisarray[4] = eigenvec(0,1); axisarray[8] = eigenvec(0,2);  axisarray[12] = 0.0;
  axisarray[1] = eigenvec(1,0); axisarray[5] = eigenvec(1,1); axisarray[9] = eigenvec(1,2);  axisarray[13] = 0.0;
  axisarray[2] = eigenvec(2,0); axisarray[6] = eigenvec(2,1); axisarray[10] = eigenvec(2,2); axisarray[14] = 0.0;
  axisarray[3] = 0.0;           axisarray[7] = 0.0;           axisarray[11] = 0.0;           axisarray[15] = 1.0;

  // 3) apply rotation matrix
  glPushMatrix();
  glMultMatrixd(axisarray);

  // 4) paint ellipsoid
  glLineWidth(0.5);
  glScalef(eigenval(0), eigenval(1), eigenval(2)); // scale ellipsoid according to std-deviation
  GLUquadricObj *myquadric = quadric;
  if (quadric == NULL) myquadric = gluNewQuadric();
  gluQuadricDrawStyle(myquadric, GLU_LINE); // rendering style of the quadric (GLU_POINT, GLU_LINE, GLU_FILL or GLU_SILHOUETTE), default: GLU_FILL
  gluSphere(myquadric, 1.0, 20, 20); // ..., slices, stacks
  if (quadric == NULL) gluDeleteQuadric(myquadric);

  // 4) change coordinate system back to Track-CS
  glPopMatrix();
};

//! render 6x6 covariance-ellipsoid around current center (0,0,0) [rx,ry,rz,tx,ty,tz]
template <class MatrixExpr>
void render3DCovar66Ellipsoid(MatrixExpr &covar, double stdDevFactor = 1.0, double yprIndicatorLength = 0.2, GLUquadricObj *quadric = NULL) {
  using namespace std;
  using namespace matrixTools;
  using namespace matrixTools::ublas;
  // assumption: covariance is a 6x6 matrix with rot-x,y,z and trans-x,y,z variances
  BOOST_ASSERT(covar.size1() == 6 && "AlgoVisualization::render3DCovarEllipsoid: covar.size1 != 6");
  BOOST_ASSERT(covar.size2() == 6 && "AlgoVisualization::render3DCovarEllipsoid: covar.size2 != 6");
//  cout << "   COVAR:" << covar;

  matrix_range< MatrixExpr > translCovar(covar, ublas::range(3, 6), ublas::range(3, 6));

  // 1) paint translation-variance indicators
  render3DCovar33Ellipsoid(translCovar, stdDevFactor, quadric);

  // 2) paint rotation-variance indicators
  double maxDev = M_PI/2.0; // 90degrees
  double minDev = 0.03; //1.7degree
  const double rRoll = 0.5*yprIndicatorLength; // radius for roll-indicator
  unsigned int drawSteps = 10;
  double alpha;
  double alphaInc;
  glLineWidth(3);
  glEnable(GL_BLEND);
  // rot-X / roll
  double stdDev = min(maxDev, max(minDev, stdDevFactor*sqrt(covar(0,0))));
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(yprIndicatorLength,0,0);
  alpha = -stdDev;
  alphaInc = 2.0*stdDev/((double)drawSteps);
  for (unsigned int i=0; i<=drawSteps; ++i) { // draws drawSteps+1 points to obtain drawSteps line segments
    glVertex3f(yprIndicatorLength,rRoll*sin(alpha),rRoll*cos(alpha));
    alpha += alphaInc;
  }
  glEnd();
//  glBegin(GL_LINE_STRIP);
//  alpha = -stdDev;
//  for (unsigned int i=0; i<=10; ++i) {
//    alpha += ((double)i)/10.0*2*stdDev;
//    glVertex3f(yprIndicatorLength,rRoll*sin(alpha),rRoll*cos(alpha));
//  }
//  glEnd();
  // rot-Y / pitch
  stdDev = min(maxDev, max(minDev, stdDevFactor*sqrt(covar(1,1))));
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(0,0,0);
  alpha = -stdDev;
  alphaInc = 2.0*stdDev/((double)drawSteps);
  for (unsigned int i=0; i<=drawSteps; ++i) { // draws drawSteps+1 points to obtain drawSteps line segments
    glVertex3f(yprIndicatorLength*cos(alpha),0,yprIndicatorLength*sin(alpha));
    alpha += alphaInc;
  }
  glEnd();
  // rot-Z / yaw
  stdDev = min(maxDev, max(minDev, stdDevFactor*sqrt(covar(2,2))));
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(0,0,0);
  alpha = -stdDev;
  alphaInc = 2.0*stdDev/((double)drawSteps);
  for (unsigned int i=0; i<=drawSteps; ++i) { // draws drawSteps+1 points to obtain drawSteps line segments
    glVertex3f(yprIndicatorLength*cos(alpha),yprIndicatorLength*sin(alpha),0);
    alpha += alphaInc;
  }
  glDisable(GL_BLEND);
  glEnd();
  glLineWidth(1);
};


#endif /*LIDARIMAGEVISUALIZATION_H_*/
