#include "LidarImageFeatures.hpp"

#include <stdexcept>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/random.hpp>


LidarImageFeatures::NormalConfidenceLookup* LidarImageFeatures::NormalConfidenceLookup::uniqueInstance = 0;

using namespace std;
using namespace matrixTools;

LidarImageFeatures::NormalConfidenceLookup* LidarImageFeatures::NormalConfidenceLookup::get(double resolHAngRAD, double resolVAngRAD, double stdDevDist, double stdDevHAngRAD, double stdDevVAngRAD)
{
//  if (uniqueInstance == 0) {
//    boost::mutex::scoped_lock lock(instanceMutex);
    if (uniqueInstance == 0) { // double check in case multithreaded apps call it in parallel
      uniqueInstance = new LidarImageFeatures::NormalConfidenceLookup(resolHAngRAD, resolVAngRAD, stdDevDist, stdDevHAngRAD, stdDevVAngRAD);
    }
//  }
  return uniqueInstance;
}

LidarImageFeatures::NormalConfidenceLookup::NormalConfidenceLookup(double resolHAngRAD_, double resolVAngRAD_, double stdDevDist_, double stdDevHAngRAD_, double stdDevVAngRAD_)
  : center(3)
  , top(3)
  , left(3)
  , normal(3)
  , covarTop3(3,3)
  , covarLeft3(3,3)
  , Samples(3,12)
  , SamplesT(12,3)
  , add(3)
{
  resolHAngRAD = resolHAngRAD_;
  resolVAngRAD = resolVAngRAD_;
  stdDevDist = stdDevDist_;
  stdDevHAngRAD = stdDevHAngRAD_;
  stdDevVAngRAD = stdDevVAngRAD_;

  tableDistance = 10.0; // distance of center point for which the lookup table is calculated
  lrange = 1.0; // start investigating at 1m from the scanner
  urange = tableDistance + 2*(tableDistance-lrange); // stop investigating at double the lower range
  lookupsize = 80; // number of entries per dimension
  lookuptableAngle = new double [lookupsize*lookupsize];
  lookuptableCovar = new double [9*lookupsize*lookupsize];
  stepsize = (urange-lrange)/(double)lookupsize;

  // try to load from default file
  if (loadFromFile()) {
    cout << "lookup-table loaded successfully" << flush;
  } else {
    // loading was not successful -> generate lookup table and save to default file
    // investigate on (tableDistance+lrange .. tableDistance+urange)
    unsigned int nbMaxMCSamplings=1000; // number of samplings in MC-simulation to generate table
    cout << "generating lookup-table at " << tableDistance << "m in range [" << lrange << "," << urange << "]..." << flush;

    //######## calculate MC-estimate of angle variance ########
    // Create a Mersenne twister random number generator that is seeded once with #seconds since 1970
    static boost::mt19937 rng(static_cast<unsigned> (std::time(0))); // MUST BE STATIC!!!!!
    // select Gaussian probability distribution
    boost::normal_distribution<double> norm_dist(0, stdDevDist);
    boost::normal_distribution<double> norm_angH(0, stdDevHAngRAD);
    boost::normal_distribution<double> norm_angV(0, stdDevVAngRAD);
    // bind random number generator to distribution, forming a function
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >  dist_sampler(rng, norm_dist);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >  angH_sampler(rng, norm_angH);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >  angV_sampler(rng, norm_angV);
    // sample from the distribution
    DVector center(3), left(3), top(3), leftRel(3), topRel(3), normal(3);
    DVector leftnoise(3), topnoise(3), leftnoiseRel(3), topnoiseRel(3), normalnoise(3);
    DMatrix normalSamples(3,nbMaxMCSamplings);
    DMatrix covar(3,3);
    center(0) = tableDistance; center(1) = 0; center(2) = 0;
    unsigned int nbTotalSamplings = 0;
    for (ih = 0; ih < lookupsize; ++ih) {
      for (iv = 0; iv < lookupsize; ++iv) {
        double dh = lrange + stepsize*(double)ih;
        double dv = lrange + stepsize*(double)iv;
        top(0) = dv; top(1) = 0; top(2) = dv*sin(resolHAngRAD);
        left(0) = dh; left(1) = dh*sin(resolVAngRAD); left(2) = 0;
        topRel = top-center;
        leftRel = left-center;
        matrixTools::cross_product(topRel/ublas::norm_2(topRel), leftRel/ublas::norm_2(leftRel), normal);
        normal /= ublas::norm_2(normal);
        // assumption: mean = 0.0, which is justified as the expected normal is defined as the calculated one
        double angSqSumRAD = 0.0;
        double lastAngVar = 0.0;
        unsigned int nbSamplings = 0;
        while (nbSamplings < nbMaxMCSamplings) {
          topnoise(0) = top(0) + dist_sampler();
          topnoise(1) = top(0)*sin(angH_sampler());
          topnoise(2) = top(0)*sin(resolVAngRAD+angV_sampler());
          leftnoise(0) = left(0) + dist_sampler();
          leftnoise(1) = left(0)*sin(resolHAngRAD+angH_sampler());
          leftnoise(2) = left(0)*sin(angV_sampler());
          topnoiseRel = topnoise-center;
          leftnoiseRel = leftnoise-center;
          matrixTools::cross_product(topnoiseRel/ublas::norm_2(topnoiseRel), leftnoiseRel/ublas::norm_2(leftnoiseRel), normalnoise);
          normalnoise /= matrixTools::ublas::norm_2(normalnoise);
          DMatrixCol(normalSamples, nbSamplings) = normalnoise - normal;
          double angRAD = acos(matrixTools::ublas::inner_prod(normal,normalnoise));
          angSqSumRAD += angRAD*angRAD;
          ++nbSamplings;
          if (nbSamplings % 50 == 0) { // test every 50 iterations for loop-exit
            double currAngVar = angSqSumRAD/(double)nbSamplings;
            if ((currAngVar-lastAngVar)/currAngVar < 0.01) // variance didn't change too much (<1%)
              break; // exit loop
            lastAngVar = currAngVar;
          }
        }
        DMatrixRange validNSamples(normalSamples, ublas::range(0,3), ublas::range(0,nbSamplings));
        DMatrix validNSamplesT = ublas::trans(validNSamples);
        covar = ublas::prod(validNSamples, validNSamplesT);
        covar /= (double)(nbSamplings);
        lookuptableAngle[iv*lookupsize+ih] = sqrt(angSqSumRAD/(double)nbSamplings);
        for (unsigned int d=0; d<9; ++d)
          lookuptableCovar[(iv*lookupsize+ih)*9+d] = covar(d/3,d%3);
        nbTotalSamplings += nbSamplings;
      }
    }
    cout << "used " << nbTotalSamplings/(lookupsize*lookupsize) << " iterations on average" << flush;
    saveToFile();
    //saveToTextFile("normalConfidence.txt");
  }
}

LidarImageFeatures::NormalConfidenceLookup::~NormalConfidenceLookup() {
  delete [] lookuptableAngle;
  delete [] lookuptableCovar;
}

bool LidarImageFeatures::NormalConfidenceLookup::loadFromFile(std::string filename) {
  try {
    ifstream infile(filename.c_str(), ios::binary);
    double vald[9]; unsigned int valui;
    infile.read(reinterpret_cast<char *>(vald),sizeof(double)*9);
    infile.read(reinterpret_cast<char *>(&valui),sizeof(unsigned int));
    if ((vald[0] != resolHAngRAD) || (vald[1] != resolVAngRAD) || (vald[2] != stdDevDist)
         || (vald[3] != stdDevHAngRAD) || (vald[4] != stdDevVAngRAD) || (vald[5] != tableDistance)
         || (vald[6] != lrange) || (vald[7] != urange) || (vald[8] != stepsize) || (valui != lookupsize)) {
      return false; // parameters of file are not the same as desired
    }
    infile.read(reinterpret_cast<char *>(lookuptableAngle),sizeof(double)*lookupsize*lookupsize);
    infile.read(reinterpret_cast<char *>(lookuptableCovar),sizeof(double)*lookupsize*lookupsize*9);
    return true;
  } catch (const exception &e) {
    return false;
  }
}

bool LidarImageFeatures::NormalConfidenceLookup::saveToFile(std::string filename) {
  try {
    ofstream outfile(filename.c_str(), ios::out | ios::binary);
    outfile.write(reinterpret_cast<char *>(&resolHAngRAD),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&resolVAngRAD),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&stdDevDist),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&stdDevHAngRAD),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&stdDevVAngRAD),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&tableDistance),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&lrange),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&urange),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&stepsize),sizeof(double));
    outfile.write(reinterpret_cast<char *>(&lookupsize),sizeof(unsigned int));
    outfile.write(reinterpret_cast<char *>(lookuptableAngle),sizeof(double)*lookupsize*lookupsize);
    outfile.write(reinterpret_cast<char *>(lookuptableCovar),sizeof(double)*lookupsize*lookupsize*9);
    return true;
  } catch (const exception &e) {
    return false;
  }
}

void LidarImageFeatures::NormalConfidenceLookup::saveToTextFile(std::string filename) {
  try {
    ofstream outfile(filename.c_str());
    for (ih = 0; ih < lookupsize; ++ih) {
      for (iv = 0; iv < lookupsize; ++iv) {
        outfile << lrange+ih*stepsize << "\t" << lrange+iv*stepsize << "\t" << lookuptableAngle[iv*lookupsize+ih]*180.0/M_PI << endl;
      }
      outfile << endl;
    }
    for (ih = 0; ih < lookupsize; ++ih) {
      for (iv = 0; iv < lookupsize; ++iv) {
        outfile << lrange+ih*stepsize << "\t" << lrange+iv*stepsize << "\t";
        for (unsigned int d = 0; d < 9; ++d) {
          outfile << "\t" << lookuptableAngle[(iv*lookupsize+ih)*9+d];
        }
        outfile << endl;
      }
      outfile << endl;
    }
  } catch (const exception &e) {
    cerr << e.what() << endl;
  }
}

double LidarImageFeatures::NormalConfidenceLookup::getStdDevRAD(double dist, double distHoriz, double distVert) const {
  // 1st) scale values as table was generated for only 1 distance
  fac = tableDistance/dist;
  dhRel = distHoriz * fac;
  dvRel = distVert * fac;
  if ((dhRel < lrange) || (dhRel >= urange) || (dvRel < lrange) || (dvRel >= urange))
    return M_PI/2.0;
  ih = (unsigned int)((dhRel-lrange) / stepsize);
  iv = (unsigned int)((dvRel-lrange) / stepsize);
  return lookuptableAngle[iv*lookupsize+ih];
}

DMatrix LidarImageFeatures::NormalConfidenceLookup::getCovar(double dist, double distHoriz, double distVert) const {
  // 1st) scale values as table was generated for only 1 distance
//  fac = tableDistance/dist;
//  dhRel = distHoriz * fac;
//  dvRel = distVert * fac;
//  if ((dhRel < lrange) || (dhRel >= urange) || (dvRel < lrange) || (dvRel >= urange))
//    return DZeroMatrix(3,3);
//  ih = (unsigned int)((dhRel-lrange) / stepsize);
//  iv = (unsigned int)((dvRel-lrange) / stepsize);
//  DMatrix ret(3,3);
//  for (unsigned int d=0; d<9; ++d)
//    ret(d/3,d%3) = lookuptableCovar[(iv*lookupsize+ih)*9+d];
//  return ret;

  // TODO (5): optimize normal-covariance-calculation -> correct lookup-table generation
  // unscented transform
  center(0) = dist;
  center(1) = 0;
  center(2) = 0;
  top(0) = distVert; // can be modified
  top(1) = 0;
  top(2) = top(0)*sin(resolHAngRAD);
  left(0) = distHoriz; // can be modified
  left(1) = left(0)*sin(resolVAngRAD);
  left(2) = 0;
   // calculate normal and sample noise
  matrixTools::cross_product((top-center)/norm_2(top-center), (left-center)/norm_2(left-center), normal);
  normal /= matrixTools::ublas::norm_2(normal);
  covarTop3 = DZeroMatrix(3,3);
  covarTop3(0,0) = pow(stdDevDist,2); //0.015^2 = 0,00225
  covarTop3(1,1) = pow(top(0)*sin(stdDevHAngRAD),2); // (8sin(0.0015))^2 = 0,00014
  covarTop3(2,2) = pow(top(0)*sin(stdDevVAngRAD),2); // (8sin(0.00017))^2 = 1,8e-6
  covarLeft3 = DZeroMatrix(3,3);
  covarLeft3(0,0) = pow(stdDevDist,2);
  covarLeft3(1,1) = pow(left(0)*sin(stdDevHAngRAD),2);
  covarLeft3(2,2) = pow(left(0)*sin(stdDevVAngRAD),2);
  for (unsigned int i=0; i<3; ++i) {
    add = DZeroVector(3);
    add(i) = sqrt(covarTop3(i,i));
    DMatrixCol sampleP(Samples, i);
    matrixTools::cross_product((top+add-center)/norm_2(top+add-center), (left-center)/norm_2(left-center), sampleP);
    sampleP /= matrixTools::ublas::norm_2(sampleP);
    DMatrixCol sampleM(Samples, i+3);
    matrixTools::cross_product((top-add-center)/norm_2(top-add-center), (left-center)/norm_2(left-center), sampleM);
    sampleM /= matrixTools::ublas::norm_2(sampleM);
  }
  for (unsigned int i=0; i<3; ++i) {
    add = DZeroVector(3);
    add(i) = sqrt(covarLeft3(i,i));
    DMatrixCol sampleP(Samples, i+6);
    matrixTools::cross_product((top-center)/norm_2(top-center), (left-center+add)/norm_2(left-center+add), sampleP);
    sampleP /= matrixTools::ublas::norm_2(sampleP);
    DMatrixCol sampleM(Samples, i+9);
    matrixTools::cross_product((top-center)/norm_2(top-center), (left-center-add)/norm_2(left-center-add), sampleM);
    sampleM /= matrixTools::ublas::norm_2(sampleM);
  }
  for (unsigned int i=0; i<12; ++i)
    DMatrixCol(Samples, i) -= normal;
  SamplesT = ublas::trans(Samples);
  return ublas::prod(Samples,SamplesT) / 12.0;
}


void LidarImageFeatures::postprocessDistanceImage(LidarImage<double> &distanceImage, double maxDist,
    int horizInterpolPixThresh, int vertInterpolPixThresh, double horizInterpolDistDiffThresh, bool smooth,
    double smoothMaxDist, LidarImage<double> *tmpBuffer) {
  // determine image size, should be representative to take resolution information out of distance image
  int hsize, vsize;
  distanceImage.getSize(hsize, vsize);

  // remove far measurements
  for (int row = 0; row < vsize; ++row) {
    assert(distanceImage.get(-1, row) == DBL_MAX);
    assert(distanceImage.get(hsize, row) == DBL_MAX);
    for (int col = 0; col < hsize; ++col) {
      if (distanceImage.get(col, row) > maxDist)
        distanceImage.set(col, row, DBL_MAX);
    }
  }

  // interpolate missing pixels horizontally
  for (int row = 0; row < vsize; ++row) {
    double oldVal = DBL_MAX;
    for (int col = 0; col < hsize; ++col) {
      double currVal = distanceImage.get(col, row);
      if ((currVal == DBL_MAX) && (oldVal != DBL_MAX)) {
        // search next valid pixel
        int nextCol = col;
        while (((currVal = distanceImage.get(nextCol, row)) == DBL_MAX) && (nextCol < hsize - 1))
          nextCol++;
        // interpolate all pixels inbetween
        if (currVal != DBL_MAX) {
          double diffPerPix = (currVal - oldVal) / (double) (nextCol - col);
          if ((diffPerPix <= horizInterpolDistDiffThresh) && ((nextCol - col) < horizInterpolPixThresh)) {
            for (int i = 1; col < nextCol; ++col, ++i)
              distanceImage.set(col, row, oldVal + diffPerPix * i);
          }
        }
        col = nextCol;
        currVal = distanceImage.get(col, row);
      }
      oldVal = currVal;
    }
  }

  // interpolate missing pixels vertically
  for (int col = 0; col < hsize; ++col) {
    double oldVal = DBL_MAX;
    for (int row = 0; row < vsize; ++row) {
      double currVal = distanceImage.get(col, row);
      if ((currVal == DBL_MAX) && (oldVal != DBL_MAX)) {
        // search next valid pixel
        int nextRow = row;
        while (((currVal = distanceImage.get(col, nextRow)) == DBL_MAX) && (nextRow < vsize - 1))
          nextRow++;
        // interpolate all pixels inbetween
        if (currVal != DBL_MAX) {
          double diffPerPix = (currVal - oldVal) / (double) (nextRow - row);
          if ((diffPerPix <= horizInterpolDistDiffThresh) && ((nextRow - row) < vertInterpolPixThresh)) {
            for (int i = 1; row < nextRow; ++row, ++i)
              distanceImage.set(col, row, oldVal + diffPerPix * i);
          }
        }
        row = nextRow;
        currVal = distanceImage.get(col, row);
      }
      oldVal = currVal;
    }
  }

  // smooth all values
  if (smooth) {
    // create temporary buffer, if it does not exist
    LidarImage<double> *ownTmpBuffer = NULL;
    if (tmpBuffer == NULL) {
      ownTmpBuffer = new LidarImage<double> (hsize, vsize);
      tmpBuffer = ownTmpBuffer;
    }

    (*tmpBuffer) = distanceImage; // copy all values
    for (int row = 0; row < vsize; ++row) {
      for (int col = 0; col < hsize; ++col) {
        double oval = tmpBuffer->get(col, row);
        if (oval != DBL_MAX) {
          double acc = oval;
          double cnt = 1.0;
          double e = tmpBuffer->get(col - 1, row);
          if ((e != DBL_MAX) && (fabs(oval - e) < smoothMaxDist)) {
            acc += e;
            cnt += 1.0;
          }
          e = tmpBuffer->get(col + 1, row);
          if ((e != DBL_MAX) && (fabs(oval - e) < smoothMaxDist)) {
            acc += e;
            cnt += 1.0;
          }
          e = tmpBuffer->get(col, row - 1);
          if ((e != DBL_MAX) && (fabs(oval - e) < smoothMaxDist)) {
            acc += e;
            cnt += 1.0;
          }
          e = tmpBuffer->get(col, row + 1);
          if ((e != DBL_MAX) && (fabs(oval - e) < smoothMaxDist)) {
            acc += e;
            cnt += 1.0;
          }
          acc /= cnt;
          distanceImage.set(col, row, acc);
        }
      }
    }
    delete ownTmpBuffer;
  } // end: smoothing
}

void LidarImageFeatures::transformTo3D(
    LidarImage<double> &distanceImage,
    LidarImage<matrixTools::DVector> &point3D,
    const LidarImageProjector &projector,
    LidarImage<matrixTools::DMatrix> *pointVariance3D,
    double stdDevDist, double stdDevHAngRAD, double stdDevVAngRAD) {
//#warning ---------------- disabled normal covariance calculation --------------------
//  normalVariance3D = NULL;


  int hsize, vsize;
  distanceImage.getSize(hsize, vsize);
  // assume point3D as same size

  double d;
  DVector pt;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      d = distanceImage.get(col, row);
      if (d == DBL_MAX) {
        point3D.set(col, row, DVector(3, DBL_MAX));
        if (pointVariance3D)
          pointVariance3D->set(col, row, DMatrix(3,3, 0.0));
      } else {
        projector.get3DCoordRel(col, row, distanceImage.get(col, row), pt);
        if (isnan(pt(0)) || isnan(pt(1)) || isnan(pt(2))) {
          cerr << endl << "col " << col << " row " << row << " dist " << distanceImage.get(col, row)
              << " was transformed to NAN" << flush;
          point3D.set(col, row, DVector(3, DBL_MAX));
          distanceImage.set(col, row, DBL_MAX);
          if (pointVariance3D)
            pointVariance3D->set(col, row, DMatrix(3,3, 0.0));
        } else {
          point3D.set(col, row, pt);
          if (pointVariance3D) {
            DMatrix covar(3,3,0.0);
            covar(0,0) = pow(stdDevDist,2);
            covar(1,1) = pow(d*sin(stdDevHAngRAD),2);
            covar(2,2) = pow(d*sin(stdDevVAngRAD),2);
            DMatrix rot(3,3); // rotation matrix for bringing the covar from laser-beam-frame into the vehicle frame
            projector.getRotMatrix(col, row, rot);
            DMatrix rotT = ublas::trans(rot);
            covar = ublas::prod(rot,covar);
            covar = ublas::prod(covar,rotT);
            pointVariance3D->set(col, row, covar);
          }
        } // point valid
      } // distance valid
    } // loop over columns
  } // loop over rows
}

void LidarImageFeatures::difference(const LidarImage<double> &first, const LidarImage<double> &second, LidarImage<
    double> &result) {
  int hsize, vsize;
  result.getSize(hsize, vsize);

  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double dc = first.get(col, row);
      double dl = second.get(col, row);
      if ((dc != DBL_MAX) && (dl != DBL_MAX))
        result.set(col, row, dc - dl);
      else
        result.set(col, row, DBL_MAX);
    }
  }
}

void LidarImageFeatures::derivative(const LidarImage<double> &distanceImage, LidarImage<double> &derivH, LidarImage<
    double> &derivV) {
  int hsize, vsize;

  // horizontally
  derivH.getSize(hsize, vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double dc = distanceImage.get(col, row);
      double dl = distanceImage.get(col + 1, row);
      if ((dc != DBL_MAX) && (dl != DBL_MAX))
        derivH.set(col, row, dc - dl);
      else
        derivH.set(col, row, DBL_MAX);
      //      double diff = distanceImage.get(col, row) - distanceImage.get(col+1, row);
      //      derivH.set(col, row, diff);
    }
  }

  // vertically
  derivV.getSize(hsize, vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      double dc = distanceImage.get(col, row);
      double dl = distanceImage.get(col, row + 1);
      if ((dc != DBL_MAX) && (dl != DBL_MAX))
        derivV.set(col, row, dc - dl);
      else
        derivV.set(col, row, DBL_MAX);
      //      double diff = distanceImage.get(col, row) - distanceImage.get(col, row+1);
      //      derivV.set(col, row, diff);
    }
  }
}

void LidarImageFeatures::connectionWeigths(const LidarImage<double> &distanceImage, const LidarImage<double> &derivH,
    const LidarImage<double> &derivV, LidarImage<double> &connWeightsH, LidarImage<double> &connWeightsV,
    double MAX_DST_PM, double MAX_DST_FAC, double A, double B, double C, double D) {
  int hsize, vsize;

  // horizontally
  connWeightsH.getSize(hsize, vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      // connection between "col" and "col+1"
      double weight = LidarImageFeatures::connectivity(distanceImage.get(col, row), distanceImage.get(col + 1, row),
          derivH.get(col, row), derivH.get(col - 1, row), derivH.get(col + 1, row),
          MAX_DST_PM, MAX_DST_FAC, A, B, C, D);
      connWeightsH.set(col, row, weight);
    }
  }

  // vertically
  connWeightsV.getSize(hsize, vsize);
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      // connection between "row" and "row+1"
      double weight = LidarImageFeatures::connectivity(distanceImage.get(col, row), distanceImage.get(col, row + 1),
          derivV.get(col, row), derivV.get(col, row - 1), derivV.get(col, row + 1),
          MAX_DST_PM, MAX_DST_FAC, A, B, C, D);
      connWeightsV.set(col, row, weight);
    }
  }
}

void LidarImageFeatures::normals(const LidarImage<double> &distanceImage, const LidarImage<matrixTools::DVector> &points, LidarImage<matrixTools::DVector> &normals,
    const LidarImage<double> *connWeightsH, const LidarImage<double> *connWeightsV,
    LidarImage<double> *normalStdDevRAD, LidarImage<matrixTools::DMatrix> *normalVariance3D,
    LidarImage<double> *normalConfidence, LidarImage<matrixTools::DVector> *tmpVec3D, LidarImage<double> *tmpDbl,
    const LidarImageProjector *projector,
    double DIST_DIFF, double W_FAC_CREATE, /*double W_FAC_SMOOTH,*/ double resolHAngRAD,
    double resolVAngRAD, double stdDevDist, double stdDevHAngRAD, double stdDevVAngRAD,
    unsigned int nbNConfMedianPasses, Neighborhood nConfNeighb, bool nConfMin, bool verbose)
{
//#warning ---------------- disabled covariance calculation --------------------
//  normalVariance3D = NULL;

  int hsize, vsize;
  normals.getSize(hsize, vsize);

  NormalConfidenceLookup *ncLookup = NULL;
  if (projector == NULL)
    normalVariance3D = NULL;
  if ((normalStdDevRAD) || (normalVariance3D)) // normal stdDev / covariance desired
    ncLookup = NormalConfidenceLookup::get(resolHAngRAD, resolVAngRAD, stdDevDist, stdDevHAngRAD, stdDevVAngRAD);
  if (verbose) cout << "...";

  double wi[5]; // weight of connection [0.0..1.0]
  double di[5]; // distance (m)
  bool bi[5]; // indicator if connection is valid [true/false]
  DVector pi[5]; // relative Cartesian coordinates (m,m,m)
  DVector pin[5]; // relative Cartesian coordinates (m,m,m), normalized to ||1||
  for (unsigned int i=0; i<5; ++i) {
    pi[i].resize(3);
    pin[i].resize(3);
  }
  DMatrix rot(3,3);
  DMatrix rotT(3,3);

//  LidarImage<matrixTools::DVector> *ownTmpVec3D = NULL;
//  if (tmpVec3D == NULL) {
//    ownTmpVec3D = new LidarImage<DVector> (hsize, vsize);
//    tmpVec3D = ownTmpVec3D;
//  }
  LidarImage<double> *ownTmpDbl = NULL;
  if (tmpDbl == NULL) {
    ownTmpDbl = new LidarImage<double> (hsize, vsize);
    tmpDbl = ownTmpDbl;
  }

  tmpVec3D = &normals; // directly store result in "normals" image
  DVector cp(3); // cross product

  // first pass over image to calculate normals and confidence values
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      // get data of point itself
      double d = distanceImage.get(col, row);
      DVector &tmpN = tmpVec3D->get(col, row);
      tmpN[0] = 0.0; tmpN[1] = 0.0; tmpN[2] = 0.0; // faster than assigning zero_vector

      bool skip = (d == DBL_MAX);
      // get data of 4 neighbors. store first connection twice, so a simple for-loop can be used for calculation
      // use connection-weights: if connecting-point is more distant, increase weight (edge of object -> let normal point outwards -> good for ICP)
      double sumcpw = 0.0;
      if (!skip) {
        const DVector &p = points.get(col, row);
        di[0] = distanceImage.get(col + 1, row);
        di[1] = distanceImage.get(col, row - 1);
        di[2] = distanceImage.get(col - 1, row);
        di[3] = distanceImage.get(col, row + 1);
        di[4] = di[0];
        bi[0] = (di[0] != DBL_MAX);
        bi[1] = (di[1] != DBL_MAX);
        bi[2] = (di[2] != DBL_MAX);
        bi[3] = (di[3] != DBL_MAX);
        bi[4] = bi[0];
        if (bi[0]) {
  //        pi[0] = points.get(col + 1, row) - p; // fucking slow!
  //        pin[0] = pi[0] / norm_2(pi[0]); // fucking slow!
          const DVector &p2 = points.get(col + 1, row);
          pi[0][0] = p2[0] - p[0];
          pi[0][1] = p2[1] - p[1];
          pi[0][2] = p2[2] - p[2];
          double pinorm = norm_2(pi[0]);
          pin[0][0] = pi[0][0] / pinorm;
          pin[0][1] = pi[0][1] / pinorm;
          pin[0][2] = pi[0][2] / pinorm;
        }
        if (bi[1]) {
  //        pi[1] = points.get(col, row - 1) - p;
  //        pin[1] = pi[1] / norm_2(pi[1]);
          const DVector &p2 = points.get(col, row - 1);
          pi[1][0] = p2[0] - p[0];
          pi[1][1] = p2[1] - p[1];
          pi[1][2] = p2[2] - p[2];
          double pinorm = norm_2(pi[1]);
          pin[1][0] = pi[1][0] / pinorm;
          pin[1][1] = pi[1][1] / pinorm;
          pin[1][2] = pi[1][2] / pinorm;
        }
        if (bi[2]) {
  //        pi[2] = points.get(col - 1, row) - p;
  //        pin[2] = pi[2] / norm_2(pi[2]);
          const DVector &p2 = points.get(col - 1, row);
          pi[2][0] = p2[0] - p[0];
          pi[2][1] = p2[1] - p[1];
          pi[2][2] = p2[2] - p[2];
          double pinorm = norm_2(pi[2]);
          pin[2][0] = pi[2][0] / pinorm;
          pin[2][1] = pi[2][1] / pinorm;
          pin[2][2] = pi[2][2] / pinorm;
        }
        if (bi[3]) {
  //        pi[3] = points.get(col, row + 1) - p;
  //        pin[3] = pi[3] / norm_2(pi[3]);
          const DVector &p2 = points.get(col, row + 1);
          pi[3][0] = p2[0] - p[0];
          pi[3][1] = p2[1] - p[1];
          pi[3][2] = p2[2] - p[2];
          double pinorm = norm_2(pi[3]);
          pin[3][0] = pi[3][0] / pinorm;
          pin[3][1] = pi[3][1] / pinorm;
          pin[3][2] = pi[3][2] / pinorm;
        }
        if (bi[4]) {
  //        pi[4] = pi[0]; // fucking slow
  //        pin[4] = pin[0]; // fucking slow
          pi[4][0] = pi[0][0];
          pi[4][1] = pi[0][1];
          pi[4][2] = pi[0][2];
          pin[4][0] = pin[0][0];
          pin[4][1] = pin[0][1];
          pin[4][2] = pin[0][2];
        }
        if (connWeightsH && connWeightsV) {
          wi = { 0.0, 0.0, 0.0, 0.0, 0.0 };
          if (bi[0]) wi[0] = (di[0] > d + DIST_DIFF) ? W_FAC_CREATE : connWeightsH->get(col, row); // to the right
          if (bi[1]) wi[1] = (di[1] > d + DIST_DIFF) ? W_FAC_CREATE : connWeightsV->get(col, row - 1); // to the top
          if (bi[2]) wi[2] = (di[2] > d + DIST_DIFF) ? W_FAC_CREATE : connWeightsH->get(col - 1, row); // to the left
          if (bi[3]) wi[3] = (di[3] > d + DIST_DIFF) ? W_FAC_CREATE : connWeightsV->get(col, row); // to the bottom
          wi[4] = wi[0];
        } else {
          wi = { 1.0, 1.0, 1.0, 1.0, 1.0 };
        }
        sumcpw = wi[0] + wi[1] + wi[2] + wi[3];
      } // end: gather neighbor-data

      skip = skip || (sumcpw < 0.01);
      // calculate 4 cross products and sum them up
      double nstdDevRAD = 0.0001; // uncertainty of normal vector in terms of angle
      DMatrix nCovar = DZeroMatrix(3,3);
      double cw = 0.0; // cumulative weight
      double mcpw = 0.0; // max weight -> used as normal confidence
      if (!skip) {
        //      cout << endl << row << "\t" << vsize << "\t" << col << "\t" << hsize << "\t" << DBL_MAX << "\t" << "\t" << "\t" << flush;
        for (unsigned int i = 0; i < 4; ++i) { // calculates right-top, top-left, left-bottom and bottom-left cross products
        //        cout << di[i] << "\t" << bi[i] << "\t" << pi[i] << di[i+1] << "\t" << bi[i+1] << "\t" << pi[i+1] << flush;
          if (bi[i] && bi[i+1]) {
            double cpw = wi[i] * wi[i+1]; // weight by connection weights
            mcpw = max(mcpw, cpw);
            cw += cpw;
            cross_product(pin[i], pin[i+1], cp); // |cp| = |pi||pi+1|sin(ang) =>  normal length ~ sin(ang)
            tmpN += cpw * cp; // normal points towards scanner
            unsigned int ih = i+(i%2); // horizontal connection index (0->0, 1/2->2, 3->4)
            unsigned int iv = i+((i+1)%2); // vertical connection index (0/1->1, 2/3->3)
            if (normalStdDevRAD)
              nstdDevRAD += cpw * ncLookup->getStdDevRAD(d, di[ih], di[iv]);
            if (normalVariance3D)
              nCovar += cpw * ncLookup->getCovar(d, di[ih], di[iv]);
          }
        } // end: loop over 4 neighbors
      } // end: calculate cross products
      double nl = norm_2(tmpN);
      assert(!isnan(tmpN(0)) && "tmpVec3D is NAN");
      assert(!isnan(tmpN(1)) && "tmpVec3D is NAN");
      assert(!isnan(tmpN(2)) && "tmpVec3D is NAN");

      // check if invalid normal
      skip = skip || (nl < 0.01);
      if (!skip) {
        // set normal vector based on sum of 4 cross-product-normal-vectors
        tmpN /= nl; // normalize normal to length 1
        nstdDevRAD /= cw;

        // calculate standard deviation of angle of normal vector
        if (normalStdDevRAD)
          normalStdDevRAD->set(col, row, nstdDevRAD);
        if (normalVariance3D) {
          nCovar /= cw;
          projector->getRotMatrix(col, row, rot);
          rotT = ublas::trans(rot);
          nCovar = ublas::prod(rot,nCovar);
          nCovar = ublas::prod(nCovar,rotT);
          normalVariance3D->set(col, row, nCovar);
        }

        // calculate confidence in normal vector estimation
        if (normalConfidence) {
          double confidence = min(1.0,mcpw); // maximum confidence of one of the 4 cross products
          double phpl = 1.0; // probability that the plane-assumption holds horizontally
          double pvpl = 1.0; // probability that the plane-assumption holds vertically
          // TODO (5): optimize: asin and exp need a lot of time to calculate (10% and 11% of this function respectively)
          if (bi[0]) phpl *= exp(-0.5*pow(asin(abs(inner_prod_3d(pin[0],tmpN)))/nstdDevRAD,2)); // to the right
          if (bi[2]) phpl *= exp(-0.5*pow(asin(abs(inner_prod_3d(pin[2],tmpN)))/nstdDevRAD,2)); // to the right
          if (bi[1]) pvpl *= exp(-0.5*pow(asin(abs(inner_prod_3d(pin[1],tmpN)))/nstdDevRAD,2)); // to the right
          if (bi[3]) pvpl *= exp(-0.5*pow(asin(abs(inner_prod_3d(pin[3],tmpN)))/nstdDevRAD,2)); // to the right
          confidence = min(confidence,max(phpl,pvpl));
          normalConfidence->set(col, row, confidence);
        } // end: calculate normal confidence
      }

      if (skip) {
        tmpN[0] = 0.0; tmpN[1] = 0.0; tmpN[2] = 0.0; // faster than assigning zero
        if (normalStdDevRAD)
          normalStdDevRAD->set(col, row, M_PI/2);
        if (normalVariance3D)
          normalVariance3D->set(col, row, DZeroMatrix(3));
        if (normalConfidence)
          normalConfidence->set(col, row, 0.0);
      }
    } // end: column
  } // end: row

  // smooth normal vectors
  /* disabled -> result was directly stored in correct image
//  DVector p;
//  double w, d;
  for (int row = 0; row < vsize; ++row) {
    for (int col = 0; col < hsize; ++col) {
      DVector &ntmp = tmpVec3D->get(col, row);
      DVector &n = normals.get(col, row);
      n[0] = ntmp[0]; n[1] = ntmp[1]; n[2] = ntmp[2]; // this is much faster than n = ntmp;
//      p = tmpVec3D->get(col, row);
      d = distanceImage.get(col, row);
      if (d == DBL_MAX) { // if point itself is invalid, set 0.0 and skip
        normals.set(col, row, ublas::zero_vector<double>(3));
        continue;
      }
      di[0] = distanceImage.get(col + 1, row);
      di[1] = distanceImage.get(col, row - 1);
      di[2] = distanceImage.get(col - 1, row);
      di[3] = distanceImage.get(col, row + 1);
      if (di[0] != DBL_MAX) {
        w = ((di[0] > d + DIST_DIFF) || (!connWeightsH)) ? W_FAC_SMOOTH : connWeightsH->get(col, row); // to the right
        p += w * tmpVec3D->get(col + 1, row);
      }
      if (di[1] != DBL_MAX) {
        w = ((di[1] > d + DIST_DIFF) || (!connWeightsV)) ? W_FAC_SMOOTH : connWeightsV->get(col, row - 1); // to the top
        p += w * tmpVec3D->get(col, row - 1);
      }
      if (di[2] != DBL_MAX) {
        w = ((di[2] > d + DIST_DIFF) || (!connWeightsH)) ? W_FAC_SMOOTH : connWeightsH->get(col - 1, row); // to the left
        p += w * tmpVec3D->get(col - 1, row);
      }
      if (di[3] != DBL_MAX) {
        w = ((di[3] > d + DIST_DIFF) || (!connWeightsV)) ? W_FAC_SMOOTH : connWeightsV->get(col, row); // to the bottom
        p += w * tmpVec3D->get(col, row + 1);
      }
//      double length = norm_2(p);
//      if (length < 0.01) // will be small if weights are low
//        normals.set(col, row, ublas::zero_vector<double>(3));
//      else
//        normals.set(col, row, p * (1 / length));
    }
  }
*/

   // smooth normal confidence
  if (normalConfidence) {
    for (unsigned int i=0; i<nbNConfMedianPasses; i++) {
      median(*normalConfidence, *tmpDbl, nConfNeighb, 0.0); // stores result in *tmpDbl
      if (nConfMin)
        minimum(*normalConfidence, *tmpDbl, *normalConfidence);
      else
        normalConfidence->copyShallowFrom(*tmpDbl); // sufficient for type double
    }
  }

  delete ownTmpDbl; // temporary lidar image, null if one was provided
//  delete ownTmpVec3D; // temporary lidar image, null if one was provided
}

void LidarImageFeatures::median(const LidarImage<double> &data, LidarImage<double> &target, double invalidVal)
{
  // implementation-idea: simply iterate over pixels making use of iterator-functions to get 4 neighbors
  assert(data.getHorizSize() == target.getHorizSize());
  assert(data.getVertSize()  == target.getVertSize() );
  LidarImage<double>::const_iterator dIt = data.begin();
  LidarImage<double>::const_iterator dEnd = data.end();
  LidarImage<double>::iterator tbIt = target.begin();
  LidarImage<double>::iterator tbEnd = target.end();
  vector<double> elements;
  elements.reserve(5);
  while ((dIt != dEnd) && (tbIt != tbEnd)) {
    elements.clear();
    if (*dIt == invalidVal) {
      *tbIt = invalidVal;
    } else {
      elements.push_back(*dIt);
      double val;
      val = dIt.getLeftPixel();
      if (val != invalidVal) elements.push_back(val);
      val = dIt.getRightPixel();
      if (val != invalidVal) elements.push_back(val);
      val = dIt.getUpperPixel();
      if (val != invalidVal) elements.push_back(val);
      val = dIt.getLowerPixel();
      if (val != invalidVal) elements.push_back(val);
      sort(elements.begin(), elements.end());
      double median = elements[elements.size()/2]; // there is at least 1 element
      *tbIt = median;
    }
    ++dIt;
    ++tbIt;
  }
}

void LidarImageFeatures::median(const LidarImage<double> &data, LidarImage<double> &target, Neighborhood nbh, double invalidVal)
{
  // implementation-idea: iterate over columns/rows, to get neighbors process a local window
  assert(data.getHorizSize() == target.getHorizSize());
  assert(data.getVertSize()  == target.getVertSize() );
  int hSize, vSize;
  data.getSize(hSize, vSize);
  vector<double> elements;
  int radius = 0; // radius of local window in pixel
  int maxRadSum = 0; // maximum sum of absolute radii, used to accept/reject diagonals
  switch (nbh) {
  case NBH4:
    radius = 1;
    maxRadSum = 1;
    elements.reserve(5);
    break;
  case NBH8:
    radius = 1;
    maxRadSum = 2;
    elements.reserve(9);
    break;
  case NBH24:
    radius = 2;
    maxRadSum = 4;
    elements.reserve(25);
    break;
  }
  for (int h=0; h<hSize; ++h) {
    for (int v=0; v<vSize; ++v) {
      if (data(h,v) == invalidVal) {
        // if pixel is invalid, directly overtake this value
        target(h,v) = invalidVal;
      } else {
        // if pixel is valid calculate median across neighborhood
        elements.clear();
        for (int ho=max(-radius,-h); ho<=min(radius,hSize-1-h); ++ho) {
          for (int vo=max(-radius,-v); vo<=min(radius,vSize-1-v); ++vo) {
            if (abs(ho)+abs(vo)>maxRadSum)
              continue;
            const double &val = data(h+ho,v+vo);
            if (val != invalidVal)
              elements.push_back(val);
          }
        }
        sort(elements.begin(), elements.end());
        target(h,v) = elements[elements.size()/2]; // there is at least 1 element
      }
    }
  }
}

