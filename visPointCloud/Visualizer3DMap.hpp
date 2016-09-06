#ifndef VISUALIZER3DMAP_H_
#define VISUALIZER3DMAP_H_

#include <cmath>
#include <deque>
#include <vector>
#include <GL/glut.h>

#include <Gui3DQt/Visualizer.hpp>
#include <Gui3DQt/PointCloudRenderer.hpp>

#include "ui_Visualizer3DMap.h"

class Visualizer3DMap : public Gui3DQt::Visualizer
{
  Q_OBJECT

public:
  Visualizer3DMap(std::string mapfile, int dummy = 0, int rangeFrom_ = -1, int rangeTo_ = -1, bool loadOnlyTraj = false, QWidget *parent = 0);
  virtual ~Visualizer3DMap();

  virtual void paintGLOpaque();
  virtual void paintGLTranslucent();

private:
  typedef Gui3DQt::PointCloudRenderer::GlVec3 GlVec3;
  typedef Gui3DQt::PointCloudRenderer::GlCol3 GlCol3;
  struct Attrib {
    unsigned char normalDirectionColorIndex; // max 256 for Hue channel (S/V=255)
    unsigned char normalConfidenceColorIndex; // max 256 for Value channel (H=S=0)
    unsigned char scannerDistColorIndex; // max 256 for Value channel (S/V=255)
  };

  Ui::Visualizer3DMapClass ui;

  std::string filename;
  int rangeFrom;
  int rangeTo;
  bool mapIsP3d;
  bool mapIsBinary;
  bool mapHasIntensity;

  Gui3DQt::PointCloudRenderer renderer;
  std::deque<GlVec3> traj;
  std::deque<GlVec3> trajINS;
  std::vector<Attrib> attribs; // buffer with all attributes used to regenerate the color3f values
  std::vector<GLuint> indices; // list-like object containing which entries/indices of the OpenGL buffers to draw
  GLuint glListIndex; // paint lists: [0]pointcloud(unused, replaced by vertex-array),
                      // [1]limit-cube, [2]coordinats, [3]traj, [4]trajIns
  unsigned int ptSize;
  float tminZ,tmaxZ; // trajectories min/max Z
  float pminX,pmaxX,pminY,pmaxY,pminZ,pmaxZ;

  void checkFileType();
  void importTrajectory(std::deque<GlVec3> &traj, std::string filename);
  void renderTraj(const std::deque<GlVec3> &traj, bool colorByHeight = false, bool renderSpheres = false);
  inline float rotateZ(const float &x, const float &y, const float &z, const float &roll, const float &pitch) { // method for approximately calculating updated z coordinate after rotation
    return z - x*tan(pitch) + y*tan(roll);
  }
  void recalcWorldMinMaxZ();


private slots:
  void sliceHeightChanged() {if (ui.rbMapSlices->isChecked()) updateMap3D();};
  void heightAngleChanged() {if (ui.rbMapSlices->isChecked() || ui.rbMapHeight->isChecked()) updateMap3D();};
  void update3D(bool doUpdate) {if (doUpdate) updateMap3D();};
  void updateMap3D();
  void updateLimit3D();
  void updateTraj3D();
  void emitStateChanged() {emit stateChanged();};
  void incPtSize() {++ptSize; updateMap3D();updateTraj3D();};
  void decPtSize() {if (ptSize > 1) --ptSize; updateMap3D();updateTraj3D();};
  void loadMap();
  void reloadTraj();
  void reloadMap();
};

#endif // VISUALIZER3DMAP_H_
