/*
 * NormalSampler.h
 *
 *  Created on: Apr 16, 2010
 *      Author: moosmann
 */

#ifndef NORMALSAMPLER_H_
#define NORMALSAMPLER_H_

#include <cvm.h>
#include "MatrixDefs.hpp"

class NormalSampler {
public:
  //template <class VectorT, class MatrixT> // could make it templated
  NormalSampler(const matrixTools::DVector &mean, const matrixTools::DMatrix &covariance);
  virtual ~NormalSampler();

  void sample(matrixTools::DVector &target);
  matrixTools::DVector sample();

  static void test(bool storetofiles = false);

private:
  const unsigned int dim;
  const matrixTools::DVector mean;
  cvm::srmatrix L;
};


#endif /* NORMALSAMPLER_H_ */
