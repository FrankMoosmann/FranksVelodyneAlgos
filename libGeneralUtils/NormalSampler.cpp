/*
 * NormalSampler.cpp
 *
 *  Created on: Apr 16, 2010
 *      Author: moosmann
 */

#include "NormalSampler.hpp"
#include "FMUtils.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace matrixTools;
using namespace std;

NormalSampler::NormalSampler(const DVector &mean_, const DMatrix &covariance)
  : dim(mean_.size())
   ,mean(mean_)
   ,L(dim)
{
//  cout << endl << mean;
//  cout << endl << covariance;
  if (covariance.size1() != dim) throw range_error("NormalSampler: specified covariance matrix has wrong dimensionality");
  if (covariance.size2() != dim) throw range_error("NormalSampler: specified covariance matrix has wrong dimensionality");
  DMatrix covar = covariance; // copy data
  double *values = &(covar.data()[0]);
  const cvm::srsmatrix covarCVM(values,dim); // doesn't copy data
//  cout << endl << covarCVM;
  L.cholesky(covarCVM); // delivers covar = Lt * L
  L.transpose(); // but we want: covar = L * Lt -> transpose
//  cout << " --> " << L << endl;
//  cout << " L*Lt = " << L * ~L << endl;
}

NormalSampler::~NormalSampler() {
}

void NormalSampler::sample(DVector &target)
{
  // cholesky decompose -> covar = L*Lt
  // random sample: draw vec_u from N(0,I) (thus each element independently from N(0,1)
  //                calculate vec_sample = mean + L * vec_u
  cvm::rvector u(dim);
  for (unsigned int j=0; j<dim; ++j)
    u(j+1) = randomGauss(0, 1);
  u = L * u;
  target.resize(dim);
  for (unsigned int j=0; j<dim; ++j)
    target(j) = mean(j) + u(j+1);
}

DVector NormalSampler::sample()
{
  DVector tmp;
  sample(tmp);
  return tmp;
}


void NormalSampler::test(bool storetofiles)
{
  // some given 2-dimensional data:
  unsigned int rows = 7;
  unsigned int dims = 2;
  double MA1[] = {-2.06667,-3.666667, 2.8,-0.666667, 5.666667,-7.333333, -10.06066,-10.06066, 20.96985,20.707107, 0.707107,0, -20.334,-18.234};
  DMatrix A1(rows,dims); for (unsigned int i=0;i<rows;++i)  for (unsigned int j=0;j<dims;++j) A1(i,j) = MA1[i*dims+j];
  if (storetofiles) {
    ofstream outfile("initialSamples.txt", ios::out);
    for (unsigned int i=0;i<rows;++i) {
      for (unsigned int j=0;j<dims;++j) {
        outfile << A1(i,j) << "\t";
      }
      outfile << endl;
    }
  }
  // calculate empirical mean and covariance:
  DVector mean = ublas::zero_vector<double>(dims);
  for (unsigned int j=0;j<dims;++j) {
    for (unsigned int i=0;i<rows;++i) {
      mean(j) += A1(i,j);
    }
    mean(j) /= (double)rows;
  }
  DMatrix centeredData = A1;
  for (unsigned int i=0;i<rows;++i) {
    for (unsigned int j=0;j<dims;++j) {
      centeredData(i,j) -= mean(j);
    }
  }
  DMatrix A1t = ublas::trans(centeredData);
  DMatrix covar = ublas::prod(A1t, centeredData)  / (double)(rows-1);

  // now sample from this distribution
  // TODO (9): calculate mean+covariance from samples and compare with initial mean+covariance
  ofstream *newsamplesout = NULL;
  if (storetofiles)
    newsamplesout = new ofstream("newSamples.txt", ios::out);
  cout << "sampling" << endl;
  NormalSampler ns(mean, covar);
  DVector sample(dims);
  cout << "now" << endl;
  for (unsigned int i=0; i<800; ++i) {
    ns.sample(sample);
    if (storetofiles) {
      for (unsigned int j=0;j<dims;++j) {
        *newsamplesout << sample(j) << "\t";
      }
      *newsamplesout << endl;
    } else {
      cout << sample << endl;
    }
  }
  if (storetofiles)
    delete newsamplesout;
}
