#include "Mapper.hpp"

#include <sys/time.h>
#include <cfloat>
#include <iostream>
#include <iomanip> // for setprecision
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/bind.hpp>

#include <HomogeneousTransformationMatrix.hpp>
#include <kogmo_time.h>
#include <FMUtils.hpp>
#include <TimeVal.hpp>
#include <PngDistanceImage.hpp>
#include <LidarImageProjector.hpp>
#include <LidarImageProjectorPNG.hpp>
#include <LidarImageFeatures.hpp>
#include "convert_coordinates.hpp"

#include <ParameterHeap.hpp>

using namespace std;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////              Mapper               ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

Mapper::Mapper(string dir, double mapResolution_, double maxDist_, const LidarImageProjector &proj_, std::string rem3DptsFile)
  : resolution(mapResolution_)
    ,maxDist(maxDist_)
    ,proj(proj_)
    ,remPtsFilename(rem3DptsFile)
    ,remPtsFile(NULL)
    ,map(NULL)
{
  reset(); // initialize remaining members
  if (rem3DptsFile.size() > 2) {
    remPtsFileIsBinary = remPtsFilename[remPtsFilename.size()-1] == 'b';
    remPtsUseIntensity = remPtsFilename[remPtsFilename.size()-1] == 'i' || remPtsFilename[remPtsFilename.size()-2] == 'i';
    cout << endl << "preparing 3d point cloud export";
    if (remPtsFileIsBinary) cout << " as binary file";
    if (remPtsUseIntensity) cout << " using intensity";
    cout << "...";
    if (remPtsFileIsBinary) 
      remPtsFile = new ofstream(rem3DptsFile.c_str(), ios::out | ios::binary);
    else
      remPtsFile = new ofstream(rem3DptsFile.c_str());
    cout << "done" << flush;
  }
  cout << endl << "searching for distance images in " << dir << endl;
  const path directory(dir);
  const boost::regex scanfilter( "scan.*\\.png" ); // "\\" is transformed into "\" at compile-time
  boost::smatch what; // match result
  if (exists(directory)) {
    directory_iterator end;
    for (directory_iterator iter(directory); iter != end; ++iter) {
      path currFile = *iter;
      if (!is_regular_file(currFile)) continue; // Skip if not a file //i->status()
      if (is_directory(currFile)) continue; // Skip if it is a directory // i->status()
      if (!boost::regex_match( currFile.leaf().string(), what, scanfilter)) continue; // Skip if no match
      string leafIntens = currFile.leaf().string();
      leafIntens.replace(0, 4, "intens"); // replace "scan" by "intens"
      path pathIntens = directory / leafIntens;
      if ((!exists(pathIntens)) || (!is_regular_file(pathIntens)) || (is_directory(pathIntens)))
        imageNames.push_back(pair<string, string>(currFile.string(), "")); // File matches, store it
      else
        imageNames.push_back(pair<string, string>(currFile.string(), pathIntens.string())); // File matches, store it
    }
  }
  sort(imageNames.begin(), imageNames.end()); // sort alphabetically
  cout << " -> found " << imageNames.size() << " files";
  if (remPtsUseIntensity) cout << " using intensity values";
  cout << endl << "searching for INS data..." << flush;
  try {
    insData = new CsvReader(dir+"/imu.cfg", ";", "", true); // file, separator, comments-indicators, extractEmptyLines
    cout << " -> found " << insData->lineCount()-insFirstLineNb+1 << " entries" << endl;
    if ((insData->lineCount() < imageNames.size()+insFirstLineNb-1) || (insData->lineCount() > imageNames.size()+insFirstLineNb))
      cerr << "mismatch between number of images (" << imageNames.size() << ") and number of lines in IMU data file ("<< insData->lineCount()+1-insFirstLineNb <<")" << endl;
  } catch (exception &e) {
    insData = NULL;
    cout << " -> not found" << endl;
  }
}


Mapper::~Mapper()
{
  delete map; // will write out the rest of the points into "remPtsFile"
  delete remPtsFile;
  delete insData;
}

void Mapper::saveMap(string filename)
{
  if (map) map->saveToFile(filename);
}

void Mapper::loadMap(string filename)
{
  if (map) map->loadFromFile(filename);
}

void Mapper::exportMap(string filename)
{
  if (map) map->exportMap(filename);
}

void Mapper::writeLanLonFile(std::string mapfilename)
{
  double lat0,lon0;
  convert_coordinates::mercator_to_latlon(mx0, my0, mercatorScale, lat0, lon0);
  string filename = mapfilename + ".posLatLon";
  ofstream out(filename.c_str());
  out << std::setprecision(12);
  out << lat0 << "\t" << lon0 << endl;
}

void Mapper::exportMapPt3D(string filename, unsigned int resolutionFactor, Map::ExportType expType)
{
  if (map) {
    map->exportMap(filename, resolutionFactor, expType);
    writeLanLonFile(filename);
  }
}

void Mapper::importMap(string filename)
{
  if (map) map->importMap(filename);
}

void Mapper::reset()
{
  delete map;
  map = new Map(resolution, maxDist);
  map->registerNotifier(boost::bind(&Mapper::mapRemovedPoint, this, _1, _2, _3, _4, _5));
  currentFrame.reset(); // creates empty,invalid pointer
  lastFrame.reset(); // creates empty,invalid pointer
  mercatorScale = 1.0;
  mx0 = 0.0;
  my0 = 0.0;
  mz0 = 0.0;
  timestampFirst = 0;
  lastFrameImsLineNb = 0;
  avgFrameDuration = 0;
}

void Mapper::readScan(unsigned int fileIdx)
{
  using namespace HomogeneousTransformationMatrix;
  ParameterHeap* params = ParameterHeap::get();

  cout << endl;
  cout << endl << "reading scan " << getImgName(fileIdx) << "..." << flush;
  ASSA::TimeVal startT = ASSA::TimeVal::gettimeofday();

  lastFrame.swap(currentFrame);

  // 1) load png-file and transform to frame
  PngDistanceImage dstImg(getImgName(fileIdx));
  int ihsize = dstImg.width();
  int ivsize = dstImg.height();
  bool isFirstFrameRead = (currentFrame.use_count() == 0);
  if (isFirstFrameRead) {
    currentFrame = LFrameSPtr(new Frame(ihsize, ivsize));
    lastFrame = LFrameSPtr(new Frame(ihsize, ivsize));
  }
  int hsize, vsize;
  currentFrame->distance.getSize(hsize,vsize); // should be representative to take resolution information out of distance image
  if (hsize != ihsize)
    throw range_error("Mapper::transformTo: objects have different horizontal resolution: "+to_string(hsize)+"/"+to_string(ihsize));
  if (vsize != ivsize)
    throw range_error("Mapper::transformTo: objects have different vertical resolution"+to_string(vsize)+"/"+to_string(ivsize));

  // Distance image generation
  currentFrame->distance.fill(DBL_MAX);
  double *distMemBlock = currentFrame->distance.getAdr(0, 0);
  dstImg.getDistances(distMemBlock, 2); // two bytes pixel boundary (one left, one right)
  double MAX_VAL = std::numeric_limits<double>::max();
  double *cDst = distMemBlock;
  for (int i=0; i<hsize*vsize; ++i, ++cDst)
    if (*cDst == 0.0) *cDst = MAX_VAL;
  // Now, invalid distance-values are MAX and invalid intensity values are 0.0

  // HACK to remove measurements at the back (invalid when car is rotating)
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < 15; ++col) {
      currentFrame->distance.set(col,row,DBL_MAX);
      currentFrame->distance.set(ihsize-1-col,row,DBL_MAX);
    }
  }

  // Intensity image generation
  if (getIntensName(fileIdx) != "") {
    currentFrame->intensity.fill(0);
    cout << "with intensity..." << flush;
    png::image<png::gray_pixel> image(getIntensName(fileIdx));
    for (int row = 0; row < vsize; ++row) {
      for (int col = 0; col < hsize; ++col) {
        if (currentFrame->distance.get(col,row) != DBL_MAX)
          currentFrame->intensity.set(col,row,image[row][col]);
      }
    }
  } else {
    currentFrame->intensity.fill(10); // make points visible at all
  }

  currentFrame->recTime = 0;
  if (lastFrame->valid) {
    currentFrame->positionHTM2w = lastFrame->positionHTM2w;
  } else {
    currentFrame->positionHTM2w = DIdMatrix(4);
  }

  // Read IMU to initialize transformations
  currentFrame->insPositionHTM2w = DIdMatrix(4);
  size_t imuLineNb = fileIdx+insFirstLineNb;
  if ((insData) && (insData->lineCount() > imuLineNb)) {
    // 0 - time
    // 1 - #sats
    // 2-4 - lat, long, alt (DEG, meter)
    // 5-7 - yaw,pitch,roll(RAD)
    // 8-13 - velocities, angular and translational
    // 14-16 - x,y,z
    if (insData->fieldCount(imuLineNb) > 0) { // record time exists
      currentFrame->recTime = KogniMobil::kogmo_timestamp_from_string(insData->get(imuLineNb,0).c_str());
      currentFrame->timeDiffToLast = currentFrame->recTime - lastFrame->recTime;
    }
    if (insData->fieldCount(imuLineNb) > 7) {
      double lat = insData->get<double>(imuLineNb,2);
      double lon = insData->get<double>(imuLineNb,3);
      double alt = insData->get<double>(imuLineNb,4);
      double yaw = insData->get<double>(imuLineNb,5);
      double pit = insData->get<double>(imuLineNb,6);
      double rol = insData->get<double>(imuLineNb,7);
      if (isFirstFrameRead) {
        mercatorScale = convert_coordinates::lat_to_scale(lat);
        convert_coordinates::latlon_to_mercator(lat, lon, mercatorScale, mx0, my0); // set origin to first position
        mz0 = alt;
        if (remPtsFile)
          writeLanLonFile(remPtsFilename);
      }
      double mx, my, mz;
      convert_coordinates::latlon_to_mercator(lat, lon, mercatorScale, mx, my);
      mx -= mx0;
      my -= my0;
      mz = alt - mz0;
      currentFrame->insPositionHTM2w = YawPitchRollXYZ_2_HTM(yaw, pit, rol, mx, my, mz);
      if (isFirstFrameRead) {
        currentFrame->positionHTM2w = currentFrame->insPositionHTM2w;
        lastFrame->positionHTM2w = currentFrame->positionHTM2w;
        lastFrame->positionHTM2v = invert_HTM(currentFrame->positionHTM2w);
      }
    }
  }

  currentFrame->positionHTM2v = invert_HTM(currentFrame->positionHTM2w);
  lastFrameImsLineNb = imuLineNb;
  currentFrame->valid = true;

  if (timestampFirst == 0)
    timestampFirst = currentFrame->recTime;
  if (currentFrame->recTime == 0)
    cout << "at " << "N/A " << "ms..." << flush;
  else
    cout << "at " << (currentFrame->recTime - timestampFirst)/1000000 << "ms..." << flush;
   cout << "and pos " << HTM_2_YawPitchRollXYZ(currentFrame->positionHTM2w) << "..." << flush;

   // 2) check for lost frames etc
  currentFrame->moveHTM = lastFrame->moveHTM;
  if ((currentFrame->recTime == 0) || (lastFrame->recTime == 0)) {
    currentFrame->diffHTMToLast = lastFrame->moveHTM; // no timings -> use last value
  } else {
    // in case timings are available, try to detect missing frames
    timeDiffBuffer[map->nbAddedScans%9] = currentFrame->timeDiffToLast;
    if (map->nbAddedScans == 1) {
      cout << "init avgTD=" << currentFrame->timeDiffToLast/1000000 << "ms..." << flush;
      for (unsigned int i=0; i<9; ++i)
        timeDiffBuffer[i] = currentFrame->timeDiffToLast;
    }
    int64_t sortedTimeDiffBuffer[9];
    memcpy(&sortedTimeDiffBuffer[0], &timeDiffBuffer[0], sizeof(int64_t)*9);
    sort(&sortedTimeDiffBuffer[0],&sortedTimeDiffBuffer[8]);
    avgFrameDuration = sortedTimeDiffBuffer[5];
    cout << "avgTD=" << avgFrameDuration/1000000 << "ms..." << flush;
    if ((double)currentFrame->timeDiffToLast < (double)avgFrameDuration/1.1) {
      // data loss at end -> proceed as usual, overwrite timestamp
      cout << endl << "DATA LOSS AT END (curr/avg=" << currentFrame->timeDiffToLast/1000000 << "ms/" << avgFrameDuration/1000000 << "ms)" << flush;
      currentFrame->timeDiffToLast = avgFrameDuration;
    }
    if (((double)currentFrame->timeDiffToLast < (double)avgFrameDuration*1.1)) {
      currentFrame->diffHTMToLast = lastFrame->moveHTM; // no frame loss -> use last value
    } else {
      double fac = (double)currentFrame->timeDiffToLast / double(avgFrameDuration); // should be > 1 as only frame losses should occur
      cout << endl << "MISSING FRAMES (fac=curr/avg=" << currentFrame->timeDiffToLast/1000000 << "ms/" << avgFrameDuration/1000000 << "ms=" << flush;
      cout << fac << ")..." << flush;
      DVector diffParams = HTM_2_YawPitchRollXYZ(lastFrame->moveHTM);
      currentFrame->diffHTMToLast = YawPitchRollXYZ_2_HTM(fac*diffParams); // extrapolate
    }
  }

  // 3) calculate features
  LidarImageFeatures::postprocessDistanceImage(currentFrame->distance,
      params->framegenMaxDistance,
      params->framegenInterpolHPixThresh,
      params->framegenInterpolVPixThresh,
      params->framegenInterpolMetPixThresh,
      params->framegenSmoothDstImg,
      params->framegenSmoothDstThresh,
      &currentFrame->tmpFloat1);
  LidarImageFeatures::transformTo3D(currentFrame->distance, currentFrame->point3D, proj);
  LidarImageFeatures::derivative(currentFrame->distance, currentFrame->distDerivativeH, currentFrame->distDerivativeV);
  LidarImageFeatures::connectionWeigths(currentFrame->distance, currentFrame->distDerivativeH, currentFrame->distDerivativeV,
      currentFrame->connectivityH, currentFrame->connectivityV,
      params->connMaxDst,
      params->connMaxDstRatio,
      params->connFacA,
      params->connFacB,
      params->connFacC,
      params->connFacD);
  LidarImageFeatures::normals(currentFrame->distance, currentFrame->point3D, currentFrame->normal3D,
      &currentFrame->connectivityH, &currentFrame->connectivityV, &currentFrame->normalStdDevRAD, NULL,
      &currentFrame->normalConfidence, &currentFrame->tmpVec3D, &currentFrame->tmpFloat1, &proj,
      2*params->lidarMeasStdDevConst,
      params->featureConnectWeightCreate,
      params->lidarHAngResolRAD,
      params->lidarVAngResolRAD,
      params->lidarDistStdDev,
      params->lidarHAngStdDevRAD,
      params->lidarVAngStdDevRAD);
  LidarImageFeatures::median(currentFrame->normalConfidence, currentFrame->tmpFloat1);
  LidarImageFeatures::median(currentFrame->normalConfidence, currentFrame->tmpFloat1); // 2nd pass
  currentFrame->tmpVec3D = currentFrame->point3D; // copy 3D position data

  ASSA::TimeVal neededT = ASSA::TimeVal::gettimeofday() - startT;
  cout << "finished" << flush;
  cout << " [" << neededT.millisec() << "ms]" << flush;
}

void Mapper::registerScan(unsigned int nbSamples, bool unwarp)
{
  ASSA::TimeVal startT = ASSA::TimeVal::gettimeofday();
  map->registerScan(currentFrame, lastFrame, nbSamples,unwarp);
  ASSA::TimeVal neededT = ASSA::TimeVal::gettimeofday() - startT;
  cout << " [" << neededT.millisec() << "ms]" << flush;
//  if (ui.cbUnwarp->isChecked())
//    map->registerScan(currentFrame, lastFrame, ui.hsSubsampling->value(),ui.cbUnwarp->isChecked());
}
void Mapper::reg1(unsigned int nbSamples, bool unwarp)
{
  map->reg1(currentFrame, lastFrame, nbSamples, unwarp);
}
void Mapper::reg2()
{
  map->reg2();
}
void Mapper::reg3()
{
  map->reg3(currentFrame);
}

void Mapper::filterMap()
{
//  map->filter();
  map->filter_init();
  map->filter_finish();
}
void Mapper::filterMap_init()
{
  map->filter_init();
}
void Mapper::filterMap_step()
{
  map->filter_next();
}
void Mapper::filterMap_finish()
{
  map->filter_finish();
}

void Mapper::addScan(bool adept, bool unwarp, bool remOldMap, bool remFarPts, bool useGPSonly)
{
  ASSA::TimeVal startT = ASSA::TimeVal::gettimeofday();

  if (useGPSonly) { // overwrite by IMU values
    currentFrame->positionHTM2w = currentFrame->insPositionHTM2w;
    currentFrame->positionHTM2v = invert_HTM(currentFrame->positionHTM2w);
  }

  // after registration or GPS-overwrite, only positionHTM2w/positionHTM2v can be considered as valid
  // --> re-estimate moveHTM and diffHTMToLast (used for unwarping and prediction of next frame)
  currentFrame->diffHTMToLast = ublas::prod(lastFrame->positionHTM2v, currentFrame->positionHTM2w);
  if ((currentFrame->timeDiffToLast == 0) || ((double)currentFrame->timeDiffToLast < (double)avgFrameDuration*1.1)) {
    currentFrame->moveHTM = currentFrame->diffHTMToLast;
  } else {
    double fac = (double)avgFrameDuration / double(currentFrame->timeDiffToLast); // should be < 1 as only frame losses should occur
    DVector diffParams = HTM_2_YawPitchRollXYZ(currentFrame->diffHTMToLast);
    currentFrame->moveHTM = YawPitchRollXYZ_2_HTM(fac*diffParams);
  }

  // debug:
  cout << "result: moved from " << HTM_2_YawPitchRollXYZ(lastFrame->positionHTM2w)
       << " to " << HTM_2_YawPitchRollXYZ(currentFrame->positionHTM2w)
       << " => diff = " << HTM_2_YawPitchRollXYZ(currentFrame->diffHTMToLast)
       << " => frameDiff = " << HTM_2_YawPitchRollXYZ(currentFrame->moveHTM) << flush;
  //    cout << " --> " << frame->positionHTM2w << "..." << flush;

  if ((map->nbAddedScans > 2) && ((double)currentFrame->timeDiffToLast >= (double)avgFrameDuration*1.1)) {
    cout << endl << "SKIPPING scan-add due to missing frames" << flush;
    return;
  }
  if (map->nbAddedScans == 1 && unwarp) { // unwarp last scan (i.e. 1st frame) and thus map
    DVector poseDiff = HTM_2_YawPitchRollXYZ(invert_HTM(currentFrame->moveHTM));
    cout << "unwarping 1st frame with " << poseDiff << "..." << flush;
    int hSize, vSize;
    lastFrame->surface.getSize(hSize, vSize);
    unsigned int hbSize = hSize / map->nbBlocks2Unwarp; // might be rounded off
    double fDec = (double)hbSize / (double)hSize; // = 1/realNbBlocks
    double fac = 1.0;
    DMatrix R; DVector t;
    for (int col = 0; col < hSize; ++col) {
      if ((col % hbSize == 0)) {
        // change current R/t
        //cout << "  f=" << fac << flush;
        DMatrix currT = lastFrame->positionHTM2v; // first transform point back to vehicle CS
        currT = ublas::prod(YawPitchRollXYZ_2_HTM(fac*poseDiff), currT); // then, apply motion compensation
        currT = ublas::prod(lastFrame->positionHTM2w,currT); // finally, transform back to world CS
        HTM_2_Rt(currT, R, t);
        fac = max(0.0,fac-fDec); // decrease influence of last position while approaching end position
      }
      for (int row = 0; row < vSize; ++row) {
        Surface::SPtr s = lastFrame->surface.get(col,row);
        if (s.get() == NULL) continue;
        s->position = ublas::prod(R,s->position) + t;
        s->origPosition = ublas::prod(R,s->origPosition) + t;
        s->normal = ublas::prod(R,s->normal);
      }
    }
    // force re-initilization of search as map was unwarped
    map->lastNbAddedScans = 0;
  }
  map->addToMap(currentFrame, adept, unwarp, remOldMap, remFarPts);

  ASSA::TimeVal neededT = ASSA::TimeVal::gettimeofday() - startT;
  cout << " [" << neededT.millisec() << "ms]" << flush;
}


void Mapper::mapRemovedPoint(double x, double y, double z, unsigned char intensity, unsigned int hitCnt)
{
  if ((remPtsFile) && (remPtsUseIntensity || (hitCnt >= 2))) {
    if (remPtsFileIsBinary) { 
      remPtsFile->write(reinterpret_cast<char *>(&x),sizeof(double)); // binary output
      remPtsFile->write(reinterpret_cast<char *>(&y),sizeof(double)); // binary output
      remPtsFile->write(reinterpret_cast<char *>(&z),sizeof(double)); // binary output
    } else
      (*remPtsFile) << x << " " << y << " " << z << " ";
    if (remPtsUseIntensity) {
      if (remPtsFileIsBinary)
        remPtsFile->write(reinterpret_cast<char *>(&intensity),sizeof(unsigned char)); // binary output
      else
        (*remPtsFile) << (int)intensity << " ";
    }
  }
}
