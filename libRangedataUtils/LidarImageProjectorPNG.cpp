#include "LidarImageProjectorPNG.hpp"

#include <cfloat>
#include <iostream>
#include <boost/filesystem.hpp>

#include <MatrixDefs.hpp>
#include <HomogeneousTransformationMatrix.hpp>
#include <PngDistanceImage.hpp>
#include <CsvReader.hpp>

using namespace std;
using namespace matrixTools;
using namespace HomogeneousTransformationMatrix;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////         PNGImageProjector         ///////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/*!
 * \brief constructor
 * \param cfgFile file holding all needed configuration settings.
 *                either: [dist-multiplier]; imgHSize; imgVSize; horizStartAngle(grad); horizStopAngle(grad); vertStartAngle(grad); vertStopAngle(grad)
 *                or:     [dist-multiplier]; imgHSize; imgVSize; horizStartAngle(grad); horizStopAngle(grad); vertAngle_1; vertAngle_2; ...; vertAngle_n
 * \throws runtime_error if the config file is invalid
 */
PNGImageProjector::PNGImageProjector(string cfgFile)
{
  CsvReader config(cfgFile, " ;", "#", false, false); //std::string filename, std::string separators = " ;", std::string commentIndicators = "", bool extractEmtpyLines = false, bool extractEmptyFields = true
  if (config.lineCount() < 1)
    throw runtime_error("PNGImageProjector: configuration file has no valid line");
  unsigned int nbEntries = config.fieldCount(0);
  //cout << endl << "found " << nbEntries << " entries in file " << cfgFile << flush;

  // several possiblitities for nbEntries:
  // 6: version 1 without dist-multiplier
  // 7: version 1 with dist-multiplier
  // 4+[1]: version 2 without dist-multiplier
  // 5+[2]: version 2 with dist-multiplier
  if (nbEntries < 6)
    throw runtime_error("PNGImageProjector: configuration file has not enough valid entries");
//  for (unsigned int i=0; i<6; ++i)
//    cout << endl << i << ": " << config.get(0,i) << " = " << config.get<unsigned int>(0,i) << flush;
  const bool isVersion1 = (nbEntries == 6) || (nbEntries == 7);
  const bool hasDistMult = ((nbEntries == 7) || (nbEntries == 5+config.get<unsigned int>(0,2)));
  const unsigned int offset = hasDistMult ? 1 : 0;
  const bool isVersion2 = (nbEntries == offset+4+config.get<unsigned int>(0,1+offset));
  if (!isVersion1 && !isVersion2)
    throw runtime_error("PNGImageProjector: configuration file has unexpected number of entries");
  cout << "recognized image config file as version " << (isVersion1 ? 1 : 2);
  if (hasDistMult) cout << " with distMult";
  cout << endl;

  const double toRAD = M_PI/180.0;
  imgHSize = config.get<unsigned int>(0,0+offset);
  imgVSize = config.get<unsigned int>(0,1+offset);
  hSpacing = NULL;
  vSpacing = NULL;
//  cout << "entries: " << nbEntries;
//  for (unsigned int i=0; i<nbEntries; ++i) cout << " " << i << ":" << config.get(0,i) << flush;
  if (isVersion1) {
    double horizStartAngleRAD_ = config.get<double>(0,2+offset)*toRAD;
    double horizStopAngleRAD_ = config.get<double>(0,3+offset)*toRAD;
    if (horizStopAngleRAD_ >= horizStartAngleRAD_) throw invalid_argument("PNGImageProjector: horizStopAngleRAD_ >= horizStartAngleRAD_");
    hSpacing = new PNGImageProjector::RegularSpacing(horizStartAngleRAD_, horizStopAngleRAD_, imgHSize);
    double vertStartAngleRAD_ = config.get<double>(0,4+offset)*toRAD;
    double vertStopAngleRAD_ = config.get<double>(0,5+offset)*toRAD;
    if (vertStopAngleRAD_ <= vertStartAngleRAD_) throw invalid_argument("PNGImageProjector: vertStopAngleRAD_ <= vertStartAngleRAD_");
    vSpacing = new PNGImageProjector::RegularSpacing(vertStartAngleRAD_, vertStopAngleRAD_, imgVSize);
  }
  if (isVersion2) {
    double horizStartAngleRAD_ = config.get<double>(0,2+offset)*toRAD;
    double horizStopAngleRAD_ = config.get<double>(0,3+offset)*toRAD;
    if (horizStopAngleRAD_ >= horizStartAngleRAD_)
      throw invalid_argument("PNGImageProjector: horizStopAngleRAD_ >= horizStartAngleRAD_");
    horizOpeningAngleRAD = horizStartAngleRAD_ - horizStopAngleRAD_;
    hSpacing = new PNGImageProjector::RegularSpacing(horizStartAngleRAD_, horizStopAngleRAD_, imgHSize);
    vector<double> vertAngles; vertAngles.reserve(imgVSize);
    for (unsigned int i=4+offset; i<nbEntries; ++i)
      vertAngles.push_back(config.get<double>(0,i)*toRAD);
    vertOpeningAngleRAD = vertAngles[0] - vertAngles[vertAngles.size()-1];
    vSpacing = new PNGImageProjector::IrregularSpacing(vertAngles.begin(), vertAngles.end(), imgVSize);
  }

  setupSphere();
//  horizStartAngleRAD = config.get<double>(0,2)*toRAD;
  //horizStartAngle(grad); horizStopAngle(grad); vertStartAngle(grad); vertStopAngle(grad)
}

/*!
 * \brief constructor for a regular lattice.

    horizontal angle increases to the left, i.e. NOT according to horizontal indexing => horizStopAngleRAD_ < horizStartAngleRAD
    vertical angle increases downward, i.e. according to vertical indexing => vertStopAngleRAD_ > vertStartAngleRAD

 * \param horizStartAngleRAD_ yaw-angle (mathematical ccw angle in radian) of the left border of the first=leftmost image column.
 * \param horizStopAngleRAD_ yaw-angle (mathematical ccw angle in radian) of the right border of the last=rightmost image column.
 * \param vertStartAngleRAD_ pitch-angle (mathematical ccw angle in radian) of the upper border of the first=uppermost image row.
 * \param vertStopAngleRAD_ pitch-angle (mathematical ccw angle in radian) of the lower border of the last=lowermost image row.
 * \param imgHSize_ horizontal image size in pixel
 * \param imgVSize_ vertical image size in pixel
 * \throws invalid_argument if (horizStopAngleRAD_ >= horizStartAngleRAD_) or (vertStopAngleRAD_ <= vertStartAngleRAD_)
 */
PNGImageProjector::PNGImageProjector(double horizStartAngleRAD_, double horizStopAngleRAD_, double vertStartAngleRAD_, double vertStopAngleRAD_, unsigned int imgHSize_, unsigned int imgVSize_)
  : imgHSize(imgHSize_)
   ,imgVSize(imgVSize_)
{
  if (horizStopAngleRAD_ >= horizStartAngleRAD_) throw invalid_argument("PNGImageProjector: horizStopAngleRAD_ >= horizStartAngleRAD_");
  if (vertStopAngleRAD_ <= vertStartAngleRAD_) throw invalid_argument("PNGImageProjector: vertStopAngleRAD_ <= vertStartAngleRAD_");

  horizOpeningAngleRAD = horizStartAngleRAD_ - horizStopAngleRAD_;
  vertOpeningAngleRAD = vertStartAngleRAD_ - vertStopAngleRAD_;
  hSpacing = new PNGImageProjector::RegularSpacing(horizStartAngleRAD_, horizStopAngleRAD_, imgHSize);
  vSpacing = new PNGImageProjector::RegularSpacing(vertStartAngleRAD_, vertStopAngleRAD_, imgVSize);
  setupSphere();
}

void PNGImageProjector::setupSphere()
{
  cout << "  allocating 3x" << imgHSize*imgVSize*sizeof(double)/1024 << "kB..." << flush;
  xUnitSphere = new double[imgHSize*imgVSize];
  yUnitSphere = new double[imgHSize*imgVSize];
  zUnitSphere = new double[imgHSize*imgVSize];
  rotMat = new double[imgHSize*imgVSize*9];

  cout << "initializing..." << flush;
  DVector unity(3,0.); unity[0] = 1.;
  DMatrix Rot;
  DVector trans;
  DVector pt;
  double yaw,pitch;
  for (unsigned int row=0; row<imgVSize; ++row) {
    for (unsigned int col=0; col<imgHSize; ++col) {
      yaw = hSpacing->getAngleRAD(col);
      pitch = vSpacing->getAngleRAD(row);
      YawPitchRollXYZ_2_Rt(yaw,pitch,0.,0.,0.,0.,Rot,trans);
      pt = ublas::prod(Rot, unity);
      unsigned int idx = row*imgHSize + col;
      xUnitSphere[idx] = pt[0];
      yUnitSphere[idx] = pt[1];
      zUnitSphere[idx] = pt[2];
      memcpy( &(rotMat[9*idx]), &(Rot.data()[0]), 9*sizeof(double) );
    }
  }
  cout << "done" << endl;
}

PNGImageProjector::~PNGImageProjector()
{
  delete [] rotMat;
  delete [] zUnitSphere;
  delete [] yUnitSphere;
  delete [] xUnitSphere;
  delete vSpacing;
  delete hSpacing;
}

unsigned int PNGImageProjector::getImgHorizSize() const
{
  return imgHSize;
}

unsigned int PNGImageProjector::getImgVertSize() const
{
  return imgVSize;
}

double PNGImageProjector::getImgHorizAngleRAD() const
{
  return horizOpeningAngleRAD;
}

double PNGImageProjector::getImgVertAngleRAD() const
{
  return vertOpeningAngleRAD;
}

bool PNGImageProjector::getImageIndexYP(double yawRAD, double pitchRAD, int &hi, int &vi) const
{
  return (hSpacing->getImageIndex(yawRAD, hi) && vSpacing->getImageIndex(pitchRAD, vi));
}


void PNGImageProjector::get3DCoordRel(int hi, int vi, double dist, double &x, double &y, double &z) const
{
  unsigned int idx = vi*imgHSize + hi;
  x = dist * xUnitSphere[idx];
  y = dist * yUnitSphere[idx];
  z = dist * zUnitSphere[idx];
}

void PNGImageProjector::getRotMatrix(int hi, int vi, matrixTools::DMatrix &mat) const
{
  BOOST_ASSERT(mat.size1() == 3);
  BOOST_ASSERT(mat.size2() == 3);
  assert(hi < (int)imgHSize && "PNGImageProjector::getRotMatrix: hi is >= imgHSize");
  assert(vi < (int)imgVSize && "PNGImageProjector::getRotMatrix: vi is >= imgVSize");
  unsigned int idx = vi*imgHSize + hi;
  memcpy( &(mat.data()[0]), &(rotMat[9*idx]), 9*sizeof(double) );
}

//////////////////////////////////////////////////////////////////////////////
//////////////////  Inner Class: Regular Spacing  ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

PNGImageProjector::RegularSpacing::RegularSpacing(double startAngleRAD_, double stopAngleRAD_, unsigned int nbPix_)
{
  cout << "  set up regular spacing with " << startAngleRAD_ << ".." << stopAngleRAD_ << ", pix " << nbPix_ << endl;
  // normalize starting angle to [0,2pi)
  while (startAngleRAD_ >= 2*M_PI) { startAngleRAD_ -= 2*M_PI; stopAngleRAD_ -= 2*M_PI;}
  while (startAngleRAD_ < 0.0) { startAngleRAD_ += 2*M_PI; stopAngleRAD_ += 2*M_PI;}
  double rangeRAD_ = stopAngleRAD_-startAngleRAD_;

  startAngleRAD = startAngleRAD_;
  rangeRAD = fabs(rangeRAD_);
  rangeSign = (rangeRAD_ < 0) ? -1.0 : 1.0;
  pixSizeRAD = rangeRAD/(double)nbPix_;
}

bool PNGImageProjector::RegularSpacing::getImageIndex(double angleRAD, int &idx)
{
  // 1) make angle relative
  angleRAD -= startAngleRAD;
  angleRAD *= rangeSign;
  // 2) normalize angle
  while (angleRAD < 0) angleRAD += 2*M_PI;
  while (angleRAD>= 2*M_PI) angleRAD -= 2*M_PI;
  // 3) check angle
  if (angleRAD > rangeRAD) return false;
  // 4) calculate index
  idx = angleRAD / pixSizeRAD;
  return true;
}

double PNGImageProjector::RegularSpacing::getAngleRAD(int idx)
{
  return startAngleRAD + rangeSign*pixSizeRAD*((double)idx+0.5); // center of pixel
}

//////////////////////////////////////////////////////////////////////////////
//////////////////  Inner Class: Irregular Spacing  //////////////////////////
//////////////////////////////////////////////////////////////////////////////

bool PNGImageProjector::IrregularSpacing::getImageIndex(double angleRAD, int &idx)
{
  if (angleRAD*pixSizeSign < (vAnglesRAD.front()-avgPixSizeRAD)*pixSizeSign) return false;
  if (angleRAD*pixSizeSign > (vAnglesRAD.back()+avgPixSizeRAD)*pixSizeSign) return false;
  idx = 0;
  double oldDst = DBL_MAX;
  while (idx < (int)vAnglesRAD.size()) {
    double cDst = fabs(angleRAD-vAnglesRAD[idx]);
    if (cDst < oldDst) {
      oldDst = cDst;
      ++idx;
    }  else
      break;
  }
  idx--;
  return true;
}

double PNGImageProjector::IrregularSpacing::getAngleRAD(int idx)
{
  return vAnglesRAD[idx];
}

