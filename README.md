# Franks Velodyne Algos
This software package contains the implemented algorithms described in [my thesis](http://digbib.ubka.uni-karlsruhe.de/volltexte/1000032359).
These are a "mapper" and a "tracker" (which also performs mapping) and two visualization tools for viewing the results.

## TODO
This software package is not perfect and a lot of things can be improved. The main ones are
- Bugfix: The tracker sometimes fails to estimate the egomotion (and hence tracking fails completely). This is a bug and can be solved by the next TODO, since the mapper does not fail at these points. 
- Get rid of duplicated code in the mapper+tracker. The mapper is supposed to be a sub-part of the tracker, so code could be reused. This would also add "de-skewing", i.e. scan-correction to the tracker.
- Write a common (for the mapper/tracker) abstract class for parsing data directories and subclasses that implement different versions of data directories
- many more, including optimizations for running time (different matrix library, use heavy parallelization e.g. using GPGPU), adding class-specific detectors/motion models, ...

## COMPILATION
This code was used under a Kubuntu 14.04 system. But you should be able to compile it on more recent versions as well.
It might even be possible to get it working on Windows, since I did not use any platform specific library.
You need have installed (might not be complete):
- QT5
- various boost libraries
- blas & lapack
- libpng++

Compile it via cmake:
```
mkdir build
cd build
cmake ..
make
```



## KNOWN ISSUES
- Runtime-Exceptions
  => link against a different lapack/blas implementation (e.g. Intel Math Kernel Library)
- Tracker: Failed self-localization (leading to thousands of wrong tracks) during sharp turns
  => see TODO (1 & 2)
- When extracting the archive on file systems that do not support symlinks, some files are probably lost.
  => Copy the following files within "tracker":
     71_42_stage1_classifier -> stage1_classifier
	 71_42_stage1_classifier.svmmodel -> stage1_classifier.svmmodel
	 71_42_stage1_hyperplane.txt -> stage1_hyperplane.txt
	 71_42_stage2_hyperplane.txt -> stage2_hyperplane.txt


  
## INPUT DATA
Range data in the form of range images. Directory structures see "exampleData".
For a description of the data format see libRangedataUtils/PngDistanceImage.hpp
or the provided full datasets on [KIT/MRT website](http://www.mrt.kit.edu/z/publ/download/velodynetracking/dataset.html)

## USAGE
`program --help` will display possible options. typical usage:

1) tracker from within sub-directory "tracker":
```
./tracker ../exampleData/DurlacherTor/
```
will process the specified directory within the GUI
```
./tracker ../exampleData/DurlacherTor/ 1 30 --outDir ../test_out/
./tracker ../exampleData/MuehlburgerTor/ 1 100 --outDir ../test_out/
```
will process the specified range in console-only mode and write the output the specified directory

2) mapper
```
mapper/mapper exampleData/DurlacherTor/ 
```
will process the specified directory within the GUI
```
mapper/mapper exampleData/DurlacherTor/ --range 1 20 40 50 --local --save test_map.p3d
```
will process frame 1-20 and 40-50 in console-only mode while keeping only a local map and write the map to the specified file


## OUTPUT DATA
1) tracker
for each input file/frame one output text-file is generated, containing for each track one line. each line holds space-separated:
  [UniqueID] : 0 corresponds to the static scene track
  age
  lastUpdate
  colorR
  colorG
  colorB
  12D-statevector : rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot specified relative to sensor coordinate system in rad/m/s
  12x12-stateCovariance
  number-of-initial-points : the first n of the points (below) are the initial points from the first frame of the track
  total-number-of-points/appearance-size
  appearance : point(x,y,z),normal(nx,ny,nz),normalConfidence
additionally the text file "world.p3d" is created, which holds the mapped static scene, format as specified below

2) mapper
the ".p3d" map-file is a text file holding space-separated a collection of (x y z) values

## VISUALIZATION
1) tracking-results
from within "visTracks":
```
./genBoundingBoxes.sh ../test_out/
./visTracks ../test_out/
```

2) mapping-results
from within "visPointCloud":
```
./visPointCloud ../test_out/test_map.p3d
```

## TRACKER-PARAMETERS
All parameters can be adjusted on the console or within libTrackingUtils/ParameterHeap.?pp

Re-Training of the SVM classifier and the linear model for track association is more difficult (but cannot be avoided when using different input data).
To do so:
1) Delete the tdata?.txt (tgz in the repo) and start the tracker in GUI-mode.
2) Press "START" in the "Image Selection" area in order to read+process the first N frames. It will stop as soon as the first track-management step is about to carry out.
3) Press "Init" in the "Track Creation" area to start looping through all the tracklets that are about to be merged/ignored
4) For each tracklet, either press "keep" if the tracklet corresponds to a new moving object or "merge" to merge it with the existing track selected in the combo box. You can use "-5" and "+" to navigate through the tracklets.
5) Once you're done, press "LOOP" in the "Image Selection" to loop to step 4) for the next frame. Process as many frames as possible.
The files tdata1.txt and tdata2.txt now contain the ground truth and can be used to retrain the needed classifier/hyperplanes (sorry, I lost my scripts, you have to re-write them by yourself).
