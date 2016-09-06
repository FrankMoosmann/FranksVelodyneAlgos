#include "PrincipalComponentAnalysis.hpp"

#include <iostream>
#include <vector>
#include <assert.h>
//#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>

using namespace std;
using namespace matrixTools;
using namespace boost::numeric::ublas;


  
PCA::PCA(const DMatrix &odata)
{
  cerr << endl << "PCA from data!" << flush;
  // check given data
  nDim = odata.size1(); // number of dimensions
  unsigned int nData = odata.size2(); // number of data elements
  assert((nData>1) && "PCA requires a data matrix with at least 2 rows!");
//  cout << endl << "doing PCA with a matrix of size " << nData << "x" << nDim << flush;
//  cout << endl << "doing PCA with the following matrix: " << odata << flush;

  // do mean calculation
  DMatrix mdata = odata; // make a copy, so MatrixCol adaptors work
  mean = zero_vector<double>(nDim);
  for (unsigned int dat=0; dat<nData; ++dat) {
    mean += DMatrixCol(mdata,dat);
  }
  mean /= (double)nData;
//  cout << endl << "mean: " << mean << flush;

  // do mean subtraction
  for (unsigned int dat=0; dat<nData; ++dat) {
    DMatrixCol(mdata,dat) -= mean;
  }
//  cout << endl << "data after mean subtraction: " << mdata << flush;

  // calculate covariance matrix and process it
  DMatrix datat = trans(mdata);
  DSMatrix covar = prod(mdata,datat);
//  cout << endl << "covar before normalization" << covar << flush;
  covar /= (double)(nData-1);
  process(covar);
} //symmetric_adaptor<matrix<double>, lower> sal (m);

PCA::PCA(const DSMatrix &covar)
{
  nDim = covar.size1(); // number of dimensions
  mean = zero_vector<double>(nDim);
  process(covar);
}

PCA::~PCA()
{
}

void PCA::process(const DSMatrix &covar)
{
//  cout << endl << "processing covar " << covar << flush;
//  cout << endl << "mean=" << mean << flush;
  // do PCA
  DVector eigenvalTmp(nDim);
  DCMatrix eigenvecTmp(covar);
  boost::numeric::bindings::lapack::syev('V', 'U', eigenvecTmp, eigenvalTmp);
//  cout << endl << "EVD successful:" << eigenvalTmp << eigenvecTmp << flush;

  // store correctly sorted results, i.e. only invert the direction, as lapack already sorts the results.
  eigenval.resize(nDim);
  eigenvec.resize(nDim,nDim);
  basechange.resize(nDim,nDim);
  for (unsigned int dim=0; dim<nDim; ++dim) {
    DMatrixCol(eigenvec,dim) = DCMatrixCol(eigenvecTmp,nDim-dim-1);
    eigenval(dim) = eigenvalTmp(nDim-dim-1);
  }
  basechange = trans(eigenvec);
//  cout << endl << "EVD successful:" << eigenval << eigenvec << flush;
}

DVector PCA::getSortedEigenValues() const
{
  return eigenval;
}

DMatrix PCA::getSortedEigenVectors() const
{
  return eigenvec;
}

DVector PCA::transform(const DVector &v) const
{
  return prod(basechange,(v-mean));
}

DVector PCA::transformBack(const DVector &v) const
{
  return prod(eigenvec,v)+mean;
}
