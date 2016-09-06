#ifndef MATRIXDEFS_H_
#define MATRIXDEFS_H_

#define BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#include <sstream>
#include <cfloat>
#include <boost/io/ios_state.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>

namespace matrixTools {
  namespace ublas = boost::numeric::ublas;
  namespace lapack = boost::numeric::bindings::lapack;

  typedef ublas::vector<double>                                DVector;
  typedef ublas::matrix<double>                                DMatrix;
  typedef ublas::matrix<double,ublas::column_major>            DCMatrix; //-- lapack needs column major!
  typedef ublas::symmetric_matrix<double>                      DSMatrix;
  typedef ublas::symmetric_matrix<double,ublas::column_major>  DCSMatrix;
  typedef ublas::banded_matrix<double,ublas::column_major>    DCBMatrix;

  typedef ublas::vector_range< DVector >              DVectorRange;
  typedef ublas::vector_range< const DVector >        DVectorRangeConst;
  typedef ublas::matrix_range< DMatrix >              DMatrixRange;
  typedef ublas::matrix_range< const DMatrix >        DMatrixRangeConst;
  typedef ublas::matrix_range< DCMatrix >              DCMatrixRange;
  typedef ublas::matrix_range< const DCMatrix >       DCMatrixRangeConst;
  typedef ublas::matrix_range< DSMatrix >             DSMatrixRange;
  typedef ublas::matrix_range< const DSMatrix >       DSMatrixRangeConst;
  typedef ublas::matrix_column< DMatrix >              DMatrixCol;
  typedef ublas::matrix_column< const DMatrix >        DMatrixColConst;
  typedef ublas::matrix_column< DCMatrix >            DCMatrixCol;
  typedef ublas::matrix_column< const DCMatrix >      DCMatrixColConst;
  typedef ublas::matrix_row< DMatrix >                DMatrixRow;
  typedef ublas::matrix_row< const DMatrix >          DMatrixRowConst;
  typedef ublas::matrix_row< DCMatrix >                DCMatrixRow;
  typedef ublas::matrix_column< DCMatrixRange >        DCMatrixRangeCol;

  typedef ublas::zero_vector<double>                  DZeroVector;
  typedef ublas::zero_matrix<double>                  DZeroMatrix;
  typedef ublas::identity_matrix<double>              DIdMatrix;

  inline DVector DZeroHomVec(unsigned int dim) {DVector v = DZeroVector(dim+1); v[dim] = 1.0; return v;};
  //inline DVector DHomVec(const DVector &v) {DVector vh = v; vh.resize(v.size()+1); vh[v.size()] = 1.0; return vh;};

  // computes the square of the euclidean norm more efficiently by omitting the root
  inline double norm_2_square(DVector v) {
    double ret = 0;
    unsigned int s = v.size();
    for (unsigned int i=0; i<s; ++i)
      ret += v(i)*v(i);
    return ret;
  }

  // computes the 3D cross product of two vectors
  template <class VectExpr1, class VectExpr2, class VectExpr3>
  inline void cross_product(const VectExpr1 &first, const VectExpr2 &second, VectExpr3 &third) {
    BOOST_ASSERT((first.size() == 3) && "cross_product: first vector has size != 3");
    BOOST_ASSERT((second.size() == 3) && "cross_product: second vector has size != 3");
    BOOST_ASSERT((third.size() == 3) && "cross_product: third vector has size != 3");
    //third.resize(3);
    third[0] = first[1]*second[2]-first[2]*second[1];
    third[1] = first[2]*second[0]-first[0]*second[2];
    third[2] = first[0]*second[1]-first[1]*second[0];
  }

  // computes the inner product of two 3D vectors
  template <class VectExpr1, class VectExpr2>
  inline double inner_prod_3d(const VectExpr1 &first, const VectExpr2 &second) {
    BOOST_ASSERT((first.size() == 3) && "cross_product: first vector has size != 3");
    BOOST_ASSERT((second.size() == 3) && "cross_product: second vector has size != 3");
    //third.resize(3);
    return first[0]*second[0] + first[1]*second[1] + first[2]*second[2];
  }

  // makes an alligned output of a matrix
  template <class MatExpr>
  inline std::string aligned_write(const MatExpr &a) {
    //boost::io::ios_flags_saver coutstate(stream);
    std::ostringstream oss;
    oss << std::scientific;
    for (unsigned int i=0; i<a.size1(); ++i) {
      for (unsigned int j=0; j<a.size2(); ++j)
        oss << a(i,j) << "\t";
      if (i < a.size1()-1)
        oss << std::endl;
    }
    return oss.str();
  }

  template <class MatrixT> // templating the parameter makes if usable for all types of matrices: DMatrix, DCMatrix, DSMatrix,..
  inline DCMatrix inv(const MatrixT &mat, double regularizationCoeff = 0.0) {
    BOOST_ASSERT(mat.size1() == mat.size2());
    unsigned int n = mat.size1();
    DCMatrix inv = mat; // copy data, as it will be modified below
    if (regularizationCoeff != 0.0)
      inv += regularizationCoeff * ublas::identity_matrix<double>(n);
    std::vector<int> ipiv(n); // pivot vector, is first filled by trf, then used by tri to inverse matrix
    lapack::getrf(inv,ipiv); // inv and ipiv will both be modified
    lapack::getri(inv,ipiv); // afterwards, "inv" is the inverse
    return inv;
  }

  template <class MatrixT> // templating the parameter makes if usable for all types of matrices: DMatrix, DCMatrix, DSMatrix,..
  inline DCMatrix invSym(const MatrixT &mat, double regularizationCoeff = 0.0) {
    BOOST_ASSERT(mat.size1() == mat.size2());
    unsigned int n = mat.size1();
    DCMatrix inv = mat; // copy data, as it will be modified below
    if (regularizationCoeff != 0.0)
      inv += regularizationCoeff * ublas::identity_matrix<double>(n);
    std::vector<int> ipiv(n); // pivot vector, is first filled by trf, then used by tri to inverse matrix
    // TODO (9): use "po..." (po=positive definite matrix) instead if "sy..." (symmetric indefinite matrix) to make it faster
    lapack::sytrf('U',inv,ipiv); // inv and ipiv will both be modified
    lapack::sytri('U',inv,ipiv); // afterwards, "inv" is the real inverse, but only the upper elements are valid!!!
    ublas::symmetric_adaptor<DCMatrix, ublas::upper> iSym(inv);
    return iSym; // copies upper matrix to lower
  }

  template<class matrix_T>
  inline double determinant(ublas::matrix_expression<matrix_T> const& mat_r)
  {
    double det = 1.0;
    matrix_T mLu(mat_r() );
    ublas::permutation_matrix<std::size_t> pivots(mat_r().size1() );
    int is_singular = lu_factorize(mLu, pivots);
    if (!is_singular) {
      for (std::size_t i=0; i < pivots.size(); ++i) {
        if (pivots(i) != i)
          det *= -1.0;
        det *= mLu(i,i);
      }
    } else {
      det = 0.0;
    }
    return det;
  }

  inline DMatrix Rt_2_HTM(const DMatrix &R, const DVector &t)
  {
    DMatrix htm = DIdMatrix(4,4);
    DMatrixRange(htm, ublas::range(0,3), ublas::range(0,3)) = R;
    DMatrixCol htmTHom(htm, 3);
    ublas::vector_range< DMatrixCol > (htmTHom, ublas::range(0,3)) = t;
    return htm;
  }

  template<class matrix_T>
  DMatrix orthonormalize(const matrix_T& mat) {
    // after the method of Augustin A. Dubrulle:
    // "An Optimum Iteration for the Matrix Polar Decomposition", Elect. Trans. Num. Anal, 1999
    BOOST_ASSERT(mat.size1() == mat.size2());
//    std::cout << "fixing " << ublas::prod(mat,ublas::trans(mat)) << " with det " << determinant(mat) << std::flush;
    DMatrix A = mat;
    DMatrix W = A;
    double eps = 1e-15;
    double limit = (1+eps)*sqrt(mat.size1());
    A = inv(ublas::trans(W));
    double g = sqrt(norm_frobenius(A)/norm_frobenius(W)); // frobenius norm
    W = 0.5*(g*W+(1/g)*A);
    double f = norm_frobenius(W);
    double pf = DBL_MAX;
    while ((f > limit) && (f < pf)) {
//      std::cout << "." << std::flush;
      pf = f;
      A = inv(ublas::trans(W));
      g = sqrt(norm_frobenius(A)/f);
      W = 0.5*(g*W+(1/g)*A);
      f = norm_frobenius(W);
    }
    // now, W^T*W=I, i.e. W^T = inv(W)
    //W /= determinant(W); ?? Necessary??
//    std::cout << "to " << ublas::prod(mat,ublas::trans(W)) << " with det " << determinant(W) << std::flush;
    return W;
  }
  inline void testOrthonormalization() {
    // r = result from HTM::YawPitchRollXYZ_2_Rt(0.1, 0.4, 0.8, 0, 0, 0, R, t);
    double ra[9] = {0.91646,0.208401,0.341571,0.0919527,0.721115,-0.686686,-0.389418,0.660729,0.641709};
    DMatrix R(3,3);
    for (unsigned int r=0; r<3; ++r)
      for (unsigned int c=0; c<3; ++c)
        R(r,c) = ra[r*3+c];
    DMatrix Rt, invR; DVector t;
    Rt = ublas::trans(R);
    invR = inv(R);
    std::cout << "R =      " << R << std::endl;
    std::cout << "det(R) = " << determinant(R) << std::endl;
    std::cout << "Rt =     " << Rt << std::endl;
    std::cout << "RtR =    " << ublas::prod(Rt,R) << std::endl;
    std::cout << "inv(R) = " << invR << std::endl;
    std::cout << "inv(R)*R=" << ublas::prod(invR,R) << std::endl;

    std::cout << std::endl << "altering R a bit" << std::endl;
    R(0,0) += 0.001;
    R(1,0) += 0.004;
    R(2,0) -= 0.01;
    Rt = ublas::trans(R);
    invR = inv(R);
    std::cout << "R =      " << R << std::endl;
    std::cout << "det(R) = " << determinant(R) << std::endl;
    std::cout << "Rt =     " << Rt << std::endl;
    std::cout << "RtR =    " << ublas::prod(Rt,R) << std::endl;
    std::cout << "inv(R) = " << invR << std::endl;
    std::cout << "inv(R)*R=" << ublas::prod(invR,R) << std::endl;

    std::cout << std::endl << "fixing R" << std::endl;
    R = orthonormalize(R);
    Rt = ublas::trans(R);
    invR = inv(R);
    std::cout << "R =      " << R << std::endl;
    std::cout << "det(R) = " << determinant(R) << std::endl;
    std::cout << "Rt =     " << Rt << std::endl;
    std::cout << "RtR =    " << ublas::prod(Rt,R) << std::endl;
    std::cout << "inv(R) = " << invR << std::endl;
    std::cout << "inv(R)*R=" << ublas::prod(invR,R) << std::endl;
  }

  class DCMatrixColConstIterator : public std::iterator<std::forward_iterator_tag, const DCMatrixCol>
  {
  private:
    DCMatrix &m;
    unsigned int idx;
  public:
    DCMatrixColConstIterator(DCMatrix &matrix, unsigned int index) : m(matrix), idx(index) {};
    DCMatrixColConstIterator(const DCMatrixColConstIterator &other) : m(other.m), idx(other.idx) {};
    ~DCMatrixColConstIterator() {};

    static DCMatrixColConstIterator begin(DCMatrix &matrix) {return DCMatrixColConstIterator(matrix, 0);};
    static DCMatrixColConstIterator end(DCMatrix &matrix) {return DCMatrixColConstIterator(matrix, matrix.size2());};

    DCMatrixColConstIterator& operator=(const DCMatrixColConstIterator& other) {m = other.m; idx = other.idx; return *this;};
    bool operator==(const DCMatrixColConstIterator& other) {return (/*(&m == &other.m) &&*/ (idx == other.idx));};
    bool operator!=(const DCMatrixColConstIterator& other) {return (/*(&m != &other.m) ||*/ (idx != other.idx));};
    DCMatrixColConstIterator& operator++() {idx++; return *this;}; //pre-increment
    DCMatrixColConstIterator operator++(int) {DCMatrixColConstIterator tmp(*this); operator++(); return tmp;}; //post-increment
    DCMatrixCol operator*() {return DCMatrixCol(m,idx);};
    DCMatrixCol operator->() {return DCMatrixCol(m,idx);};
  };

} // end namespace matrixTools


#endif /*MATRIXDEFS_H_*/
