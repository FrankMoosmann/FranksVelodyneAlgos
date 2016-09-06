#ifndef HomogeneousTransformationMatrix_hpp__
#define HomogeneousTransformationMatrix_hpp__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <assert.h>

// deprecated: usage of out-of-namespace typed
// typedef boost::numeric::ublas::matrix<float> ublas_matrix_f;
// typedef boost::numeric::ublas::matrix<double> ublas_matrix_d;
// typedef boost::numeric::ublas::vector<float> ublas_vec_f;
// typedef boost::numeric::ublas::vector<double> ublas_vec_d;

namespace HomogeneousTransformationMatrix {  // HTM = 4x4 homogeneous transformation matrix

  typedef boost::numeric::ublas::matrix<float> ublas_matrix_f;
  typedef boost::numeric::ublas::matrix<double> ublas_matrix_d;
  typedef boost::numeric::ublas::vector<float> ublas_vec_f;
  typedef boost::numeric::ublas::vector<double> ublas_vec_d;

  // given a base coordinate system "base_cs" that is first translated by tx,ty,tz,
  // then rotated by yaw,pitch,roll (in rad!) results in the new coord-sys "new_cs"
  // the following functions will create HTMs for the two tasks:
  
  // with the resulting HTM, a point in the new CS p_new_cs = (x,y,z,1) can be transformed into p_base_cs by   p_base_cs = HTM * p_new_cs
  inline static ublas_matrix_d YawPitchRollXYZ_2_HTM( const double yawRAD, const double pitchRAD, const double rollRAD, const double tx, const double ty, const double tz );
  inline static ublas_matrix_d YawPitchRollXYZ_2_HTM( const ublas_vec_d &params );
  // with the resulting HTM, a point p_base_cs = (x,y,z,1) can be transformed into the transformed CS by   p_new_cs = HTM * p_base_cs
  inline static ublas_matrix_d YawPitchRollXYZ_2_HTM_inv( const double yawRAD, const double pitchRAD, const double rollRAD, const double tx, const double ty, const double tz );
  // inverts the transformation
  inline static ublas_matrix_d invert_HTM( const ublas_matrix_d& matrix );
  // extract single transformation parameters from the HTM
  inline static void HTM_2_YawPitchRollXYZ( const ublas_matrix_d &matrix, double &yawRAD, double &pitchRAD, double &rollRAD, double &tx, double &ty, double &tz );
  inline static ublas_vec_d HTM_2_YawPitchRollXYZ( const ublas_matrix_d &matrix);


  // the same for non-homogeneous coordinates:
  inline static void YawPitchRollXYZ_2_Rt( const double yawRAD, const double pitchRAD, const double rollRAD, const double tx, const double ty, const double tz, ublas_matrix_d &R, ublas_vec_d &t );
  inline static void YawPitchRollXYZ_2_Rt( const ublas_vec_d &params, ublas_matrix_d &R, ublas_vec_d &t );
  inline static void YawPitchRollXYZ_2_Rt_inv( const double yawRAD, const double pitchRAD, const double rollRAD, const double tx, const double ty, const double tz, ublas_matrix_d &R, ublas_vec_d &t );
  inline static void invert_Rt( ublas_matrix_d &R, ublas_vec_d &t );
  inline static void Rt_2_YawPitchRollXYZ( const ublas_matrix_d &R, const ublas_vec_d &t, double &yawRAD, double &pitchRAD, double &rollRAD, double &tx, double &ty, double &tz );

  // conversion:
  inline static ublas_matrix_d Rt_2_HTM(const ublas_matrix_d &R, const ublas_vec_d &t);
  inline static void HTM_2_Rt(const ublas_matrix_d &htm, ublas_matrix_d &R, ublas_vec_d &t);
};

HomogeneousTransformationMatrix::ublas_matrix_d HomogeneousTransformationMatrix::YawPitchRollXYZ_2_HTM( const double yaw, const double pitch, const double roll, const double tx, const double ty, const double tz )
{
  double c_a = cos( yaw );
  double s_a = sin( yaw );

  double c_b = cos( pitch );
  double s_b = sin( pitch );

  double c_c = cos( roll );
  double s_c = sin( roll );

  HomogeneousTransformationMatrix::ublas_matrix_d trafo( 4, 4 ); //(row,col)

  trafo(0, 0) = c_a * c_b;
  trafo(1, 0) = s_a * c_b;
  trafo(2, 0) = -s_b;
  trafo(3, 0) = 0;
  
  trafo(0, 1) = c_a * s_b * s_c - s_a * c_c;
  trafo(1, 1) = s_a * s_b * s_c + c_a * c_c;
  trafo(2, 1) = c_b * s_c;
  trafo(3, 1) = 0;
  
  trafo(0, 2) = c_a * s_b * c_c + s_a * s_c;
  trafo(1, 2) = s_a * s_b * c_c - c_a * s_c;
  trafo(2, 2) = c_b * c_c;
  trafo(3, 2) = 0;

  trafo(0, 3) = tx;
  trafo(1, 3) = ty;
  trafo(2, 3) = tz;
  trafo(3, 3) = 1.;

  return trafo;
}

HomogeneousTransformationMatrix::ublas_matrix_d HomogeneousTransformationMatrix::YawPitchRollXYZ_2_HTM( const HomogeneousTransformationMatrix::ublas_vec_d &v )
{
  assert((v.size() == 6) && "YawPitchRollXYZ_2_HTM: parameter vector is not 6D");
  return YawPitchRollXYZ_2_HTM( v[0], v[1], v[2], v[3], v[4], v[5] );
}

HomogeneousTransformationMatrix::ublas_matrix_d HomogeneousTransformationMatrix::YawPitchRollXYZ_2_HTM_inv( const double yaw, const double pitch, const double roll, const double tx, const double ty, const double tz )
{
  double c_a = cos( yaw );
  double s_a = sin( yaw );

  double c_b = cos( pitch );
  double s_b = sin( pitch );

  double c_c = cos( roll );
  double s_c = sin( roll );

  HomogeneousTransformationMatrix::ublas_matrix_d trafo( 4, 4 );

  trafo(0, 0) = c_a * c_b;
  trafo(0, 1) = s_a * c_b;
  trafo(0, 2) = -s_b;
  trafo(0, 3) = 0;
  
  trafo(1, 0) = c_a * s_b * s_c - s_a * c_c;
  trafo(1, 1) = s_a * s_b * s_c + c_a * c_c;
  trafo(1, 2) = c_b * s_c;
  trafo(1, 3) = 0;
  
  trafo(2, 0) = c_a * s_b * c_c + s_a * s_c;
  trafo(2, 1) = s_a * s_b * c_c - c_a * s_c;
  trafo(2, 2) = c_b * c_c;
  trafo(2, 3) = 0;

  trafo(0, 3) = -tx * trafo(0, 0) - ty * trafo(0, 1) - tz * trafo(0, 2);
  trafo(1, 3) = -tx * trafo(1, 0) - ty * trafo(1, 1) - tz * trafo(1, 2);
  trafo(2, 3) = -tx * trafo(2, 0) - ty * trafo(2, 1) - tz * trafo(2, 2);
  trafo(3, 3) = 1.;

  return trafo;
}


HomogeneousTransformationMatrix::ublas_matrix_d HomogeneousTransformationMatrix::invert_HTM( const HomogeneousTransformationMatrix::ublas_matrix_d& matrix )
{
  assert((matrix.size1() == 4) && (matrix.size2() == 4) && "invert_HTM: input matrix is not 4x4 matrix");
  HomogeneousTransformationMatrix::ublas_matrix_d inverse(4, 4);

  inverse(0, 0) = matrix(0, 0);
  inverse(1, 1) = matrix(1, 1);
  inverse(2, 2) = matrix(2, 2);

  inverse(1, 0) = matrix(0, 1);
  inverse(2, 0) = matrix(0, 2);

  inverse(0, 1) = matrix(1, 0);
  inverse(2, 1) = matrix(1, 2);

  inverse(0, 2) = matrix(2, 0);
  inverse(1, 2) = matrix(2, 1);

  inverse(0, 3) = -( inverse(0, 0) * matrix(0, 3) + inverse(0, 1) * matrix(1, 3) + inverse(0, 2) * matrix(2, 3) );
  inverse(1, 3) = -( inverse(1, 0) * matrix(0, 3) + inverse(1, 1) * matrix(1, 3) + inverse(1, 2) * matrix(2, 3) );
  inverse(2, 3) = -( inverse(2, 0) * matrix(0, 3) + inverse(2, 1) * matrix(1, 3) + inverse(2, 2) * matrix(2, 3) );
  inverse(3, 3) = 1.;

  inverse(3, 0) = 0.;
  inverse(3, 1) = 0.;
  inverse(3, 2) = 0.;
    
  //  std::cout << inverse << "\n";

  return inverse;
}


void HomogeneousTransformationMatrix::HTM_2_YawPitchRollXYZ( const HomogeneousTransformationMatrix::ublas_matrix_d &matrix, double &yaw, double &pitch, double &roll, double &tx, double &ty, double &tz )
{
  assert((matrix.size1() == 4) && (matrix.size2() == 4) && "HTM_2_YawPitchRollXYZ: input matrix is not 4x4 matrix");
// YawPitchRoll_2_HTM (a=yaw, b=pitch, c=roll):
// trafo(0, 0) = c_a * c_b;
// trafo(1, 0) = s_a * c_b;
// trafo(2, 0) = -s_b;
// trafo(2, 1) = c_b * s_c;
// trafo(2, 2) = c_b * c_c;
// in case |pitch| ~= pi --> s_b
// trafo(0, 1) = c_a * s_c - s_a * c_c;
// trafo(1, 1) = s_a * s_c + c_a * c_c;
// trafo(0, 2) = c_a * c_c + s_a * s_c;
// trafo(1, 2) = s_a * c_c - c_a * s_c;

  //pitch = atan2(sqrt(matrix(2,1)*matrix(2,1)+matrix(2,2)*matrix(2,2)), -matrix(2,0)); //always correct!!!
  //BUG CORRECTED: 07.08.08: enumerator and denumerator swaped (-->next line)
  pitch = atan2(-matrix(2,0), sqrt(matrix(2,1)*matrix(2,1)+matrix(2,2)*matrix(2,2)));
  if ((M_PI - fabs(pitch)) > 0.0001) { // no singularity!
    yaw = atan2(matrix(1,0), matrix(0,0));
    roll = atan2(matrix(2,1), matrix(2,2));
  } else { // singularity! in this case either yaw or roll can be arbitrarily selected
    yaw = 0;
    roll = atan2(matrix(0,1), matrix(1,1));
  }

  tx = matrix(0,3);
  ty = matrix(1,3);
  tz = matrix(2,3);

}

inline static HomogeneousTransformationMatrix::ublas_vec_d HomogeneousTransformationMatrix::HTM_2_YawPitchRollXYZ( const HomogeneousTransformationMatrix::ublas_matrix_d &matrix)
{
  ublas_vec_d v(6);
  HTM_2_YawPitchRollXYZ( matrix, v[0], v[1], v[2], v[3], v[4], v[5] );
  return v;
}

inline static void HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt( const double yaw, const double pitch, const double roll, const double tx, const double ty, const double tz, HomogeneousTransformationMatrix::ublas_matrix_d &R, HomogeneousTransformationMatrix::ublas_vec_d &t )
{
  double c_a = cos( yaw );
  double s_a = sin( yaw );
  double c_b = cos( pitch );
  double s_b = sin( pitch );
  double c_c = cos( roll );
  double s_c = sin( roll );

  R.resize( 3, 3 ); //(row,col)

  R(0, 0) = c_a * c_b;
  R(1, 0) = s_a * c_b;
  R(2, 0) = -s_b;
  
  R(0, 1) = c_a * s_b * s_c - s_a * c_c;
  R(1, 1) = s_a * s_b * s_c + c_a * c_c;
  R(2, 1) = c_b * s_c;
  
  R(0, 2) = c_a * s_b * c_c + s_a * s_c;
  R(1, 2) = s_a * s_b * c_c - c_a * s_c;
  R(2, 2) = c_b * c_c;

  t.resize( 3 );
  t(0) = tx;
  t(1) = ty;
  t(2) = tz;
}

inline static void HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt( const HomogeneousTransformationMatrix::ublas_vec_d &v, HomogeneousTransformationMatrix::ublas_matrix_d &R, HomogeneousTransformationMatrix::ublas_vec_d &t )
{
  assert((v.size() == 6) && "YawPitchRollXYZ_2_HTM: parameter vector is not 6D");
  return YawPitchRollXYZ_2_Rt( v[0], v[1], v[2], v[3], v[4], v[5], R, t );
}

inline static void HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt_inv( const double yaw, const double pitch, const double roll, const double tx, const double ty, const double tz, HomogeneousTransformationMatrix::ublas_matrix_d &R, HomogeneousTransformationMatrix::ublas_vec_d &t )
{
  double c_a = cos( yaw );
  double s_a = sin( yaw );
  double c_b = cos( pitch );
  double s_b = sin( pitch );
  double c_c = cos( roll );
  double s_c = sin( roll );

  R.resize( 3, 3 );

  R(0, 0) = c_a * c_b;
  R(0, 1) = s_a * c_b;
  R(0, 2) = -s_b;
  
  R(1, 0) = c_a * s_b * s_c - s_a * c_c;
  R(1, 1) = s_a * s_b * s_c + c_a * c_c;
  R(1, 2) = c_b * s_c;
  
  R(2, 0) = c_a * s_b * c_c + s_a * s_c;
  R(2, 1) = s_a * s_b * c_c - c_a * s_c;
  R(2, 2) = c_b * c_c;

  t.resize(3);
  t(0)=-tx;
  t(1)=-ty;
  t(2)=-tz;
  t = prod(R,t);
}

inline static void HomogeneousTransformationMatrix::invert_Rt( HomogeneousTransformationMatrix::ublas_matrix_d &R, HomogeneousTransformationMatrix::ublas_vec_d &t )
{
  assert((R.size1() == 3) && (R.size2() == 3) && "invert_Rt: input matrix is not 3x3 matrix");
  assert((t.size() == 3) && "invert_Rt: input translation vector is 3 dimensional");
  R = boost::numeric::ublas::trans(R);
  t = -boost::numeric::ublas::prod(R,t);
}

inline static void HomogeneousTransformationMatrix::Rt_2_YawPitchRollXYZ( const HomogeneousTransformationMatrix::ublas_matrix_d &R, const HomogeneousTransformationMatrix::ublas_vec_d &t, double &yaw, double &pitch, double &roll, double &tx, double &ty, double &tz )
{
  assert((R.size1() == 3) && (R.size2() == 3) && "invert_Rt: input matrix is not 3x3 matrix");
  assert((t.size() == 3) && "invert_Rt: input translation vector is 3 dimensional");
  pitch = atan2(-R(2,0), sqrt(R(2,1)*R(2,1)+R(2,2)*R(2,2)));
  if ((M_PI - fabs(pitch)) > 0.0001) { // no singularity!
    yaw = atan2(R(1,0), R(0,0));
    roll = atan2(R(2,1), R(2,2));
  } else { // singularity! in this case either yaw or roll can be arbitrarily selected
    yaw = 0;
    roll = atan2(R(0,1), R(1,1));
  }

  tx = t(0);
  ty = t(1);
  tz = t(2);
}

inline static HomogeneousTransformationMatrix::ublas_matrix_d HomogeneousTransformationMatrix::Rt_2_HTM(const HomogeneousTransformationMatrix::ublas_matrix_d &R, const HomogeneousTransformationMatrix::ublas_vec_d &t)
{
  using namespace boost::numeric::ublas;
  ublas_matrix_d htm = identity_matrix<double>(4,4);
  matrix_range< ublas_matrix_d > htmR(htm, range(0,3), range(0,3));
  matrix_column<ublas_matrix_d> htmTHom(htm, 3);
  vector_range< matrix_column<ublas_matrix_d> > htmT(htmTHom, range(0,3));
  htmR = R;
  htmT = t;
  return htm;
}

inline static void HomogeneousTransformationMatrix::HTM_2_Rt(const HomogeneousTransformationMatrix::ublas_matrix_d &htm, HomogeneousTransformationMatrix::ublas_matrix_d &R, HomogeneousTransformationMatrix::ublas_vec_d &t)
{
  using namespace boost::numeric::ublas;
  R = matrix_range< const ublas_matrix_d >(htm, range(0,3), range(0,3));
  matrix_column<const ublas_matrix_d> tHom(htm, 3);
  t = vector_range< matrix_column<const ublas_matrix_d> >(tHom, range(0,3));
}

#endif
