#!/bin/bash
for dir in libGui3DQt/Gui3DQt/ libGeneralUtils libRangedataUtils libTrackingUtils mapper tracker visPointCloud visTracks; do
  pushd $dir
  qmake-qt4
  make $@ # e.g. for cleaning all: ./makeAll.sh clean
  popd
done
