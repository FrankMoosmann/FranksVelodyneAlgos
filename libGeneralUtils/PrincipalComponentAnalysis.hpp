/*! 
    \file   PCA.h
    \brief  Provides a class for doing Principal Component Analysis
    \author  Frank Moosmann (<moosmann@mrt.uka.de>),
    \date    2008

    Copyright: Institute for Measurement and Control Systems
               University of Karlsruhe. All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef PCA_H_
#define PCA_H_

#include "MatrixDefs.hpp"


/*!
  \class PCA
  \brief Class for doing Principal Component Analysis
  
  Given several data points, PCA determines a new coordinate system,
  where the first basis corresponds to the direction that has most variance within the data,
  the second basis to the direction with the second most variance, etc.
  It can be used to reduce multidimensional data sets to lower dimensions
  by doing PCA, then transforming each data point with the transform method to the new
  basis, and finally cutting the last few dimensions (that have low eigenvalues)

  Example:
  int main (int argc, char **argv)
  {
    int n=6;
    int dim=2;
    matrixTools::DMatrix data(dim,n); // 3 rows, 4 columns
    data(0,0) = 1; data(0,1) = 2; data(0,2) = 0; data(0,3) = 1; data(0,4) = 2; data(0,5) = 0;
    data(1,0) = 3; data(1,1) = 2; data(1,2) = 2; data(1,3) =-1; data(1,4) = 0; data(1,5) = 0;
    cout << "data: " << endl << data << endl;

    PCA pca(data);
    matrixTools::DVector eigenval = pca.getSortedEigenValues();
    matrixTools::DMatrix eigenvec = pca.getSortedEigenVectors();
    cout << endl << "eigenval: " << endl << eigenval << flush;
    cout << endl << "eigenvec: " << endl << eigenvec << flush;

    matrixTools::DVector testv(dim);
    testv(0) = 2;  testv(1) = 0;
    cout << endl << "transform " << testv;
    cout << " --> " << pca.transform(testv) << flush;

    cout << endl;
    return 0;
  }
*/
class PCA
{
public:
  PCA(const matrixTools::DMatrix &data); //!< Constructor \param data holds the data with the number of rows corresponding to dimensionality of data.
  PCA(const matrixTools::DSMatrix &covar); //!< Constructor \param covar holds the (square, symmetric) covariance matrix.
  virtual ~PCA();
  
  matrixTools::DVector getSortedEigenValues() const; //!< Returns the Eigenvalues (i.e. weights) to the corresponding eigenvectors. The values are sorted, with the highest eigenvalue first. Eigenvalues represent the estimated variances along the new axis, i.e. principal components
  matrixTools::DMatrix getSortedEigenVectors() const; //!< Returns the Eigenvectors. The first vector corresponds to the direction with the highest variance of the data (thus highest eigenvalue).
  matrixTools::DVector transform(const matrixTools::DVector &v) const; //!< transforms the vector v given in the old coordinate systen into the new coordinate system represented by the eigenvectors.
  matrixTools::DVector transformBack(const matrixTools::DVector &v) const; //!< transforms the vector v given in the new coordinate systen represented by the eigenvectors into the old coordinate system.

private:
  void process(const matrixTools::DSMatrix &covar); //!< computes everything from the covariance matrix
  unsigned int nDim;
  matrixTools::DVector mean; //!< mean of the original data
  matrixTools::DVector eigenval; //!< holds the sorted eigenvalues (biggest first)
  matrixTools::DMatrix eigenvec; //!< holds the eigenvectors. Can be used as rotation matrix
  matrixTools::DMatrix basechange; //!< holds the transposed eigenvectors. Can be used as rotation matrix
};


#endif /*PCA_H_*/
