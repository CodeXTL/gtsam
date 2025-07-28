/**
 * @file  SL4.h
 * @brief Projective Special Linear Group (SL(4, R)) factor
 * @author: Hyungtae Lim
 */

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/base/OptionalJacobian.h>
#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam_unstable/dllexport.h>

// To use exp(), log()
#include <string>
#include <unsupported/Eigen/MatrixFunctions>

using SL4Jacobian = gtsam::OptionalJacobian<15, 15>;

inline Eigen::Matrix<double, 15, 15> I_15x15 =
    Eigen::Matrix<double, 15, 15>::Identity();

using Matrix15x15 = Eigen::Matrix<double, 15, 15>;
using Matrix16x16 = Eigen::Matrix<double, 16, 16>;

inline Eigen::Matrix<double, 16, 15> setVecToAlgMatrix() {
  Eigen::Matrix<double, 16, 15> alg = Eigen::Matrix<double, 16, 15>::Zero();

  // 12 Off-diagonal E_ij generators
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (i != j) {
        alg(i * 4 + j, k++) = 1.0;
      }
    }
  }

  // For Diagonal generators B1 = diag(1, -1, 0, 0)
  alg(0, 12) = 1.0;
  alg(5, 12) = -1.0;

  // For B2 = diag(0, 1, -1, 0)
  alg(5, 13) = 1.0;
  alg(10, 13) = -1.0;

  // For B3 = diag(0, 0, 1, -1)
  alg(10, 14) = 1.0;
  alg(15, 14) = -1.0;

  return alg;
}

inline Eigen::Matrix<double, 15, 16> setAlgtoVecMatrix() {
  Eigen::Matrix<double, 15, 16> mat;
  mat << 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
      1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
      0., 0., 0., 0., 0., -1.;
  return mat;
}

// ALG_TO_VEC * VEC_TO_ALG is equals to I_15x15
const Eigen::Matrix<double, 16, 15> VEC_TO_ALG = setVecToAlgMatrix();
const Eigen::Matrix<double, 15, 16> ALG_TO_VEC = setAlgtoVecMatrix();

namespace gtsam {
// NOTE(hlim): Strictly speaking, it should be expressed as SL(4, ℝ),
// but for simplicity, we omit ℝ, assuming our target is over the real numbers.
// And the variable `sl4` represents SL(4, ℝ).
class GTSAM_EXPORT SL4 : public LieGroup<SL4, 15> {
 public:
  static const size_t dimension = 15;

 protected:
  Matrix44 T_;

 public:
  /// @name Standard Constructors
  /// @{

  /// Default constructor initializes at origin
  SL4() : T_(Matrix44::Identity()) {}

  /// Copy constructor
  SL4(const Matrix44& pose);

  SL4(const SL4& pose);

  /** print with optional string */
  void print(const std::string& s = "") const;

  /** assert equality up to a tolerance */
  bool equals(const SL4& sl4, double tol = 1e-9) const;

  inline const Matrix44& matrix() const { return T_; }

  /// @}
  /// @name Manifold
  /// @{

  /// Dimensionality of tangent space = 15 DOF - used to autodetect sizes
  inline static size_t Dim() { return dimension; }

  /// Dimensionality of tangent space = 15 DOF
  inline size_t dim() const { return dimension; }

  SL4 retract(const Vector& v, SL4Jacobian Horigin, SL4Jacobian Hv) const;

  Vector localCoordinates(const SL4& sl4, SL4Jacobian Horigin,
                          SL4Jacobian Hp2) const;

  /// @}
  /// @name Group
  /// @{

  /// identity for group operation
  static SL4 Identity() { return SL4(); }

  /// inverse transformation with derivatives
  SL4 inverse(SL4Jacobian H1 = {}) const;

  SL4 compose(const SL4& sl4, SL4Jacobian H1 = {}, SL4Jacobian H2 = {}) const;

  /// compose syntactic sugar
  inline SL4 operator*(const SL4& sl4) const { return compose(sl4); }

  /**
   * Return relative pose between p1 and p2, in p1 coordinate frame
   * as well as optionally the derivatives
   */
  SL4 between(const SL4& sl4, SL4Jacobian H1 = {}, SL4Jacobian H2 = {}) const;

  static Matrix44 Hat(const Vector& xi) {
    assert(xi.size() == 15);
    Matrix44 mat;
    const double d11 = xi(12);
    const double d22 = -xi(12) + xi(13);
    const double d33 = -xi(13) + xi(14);
    const double d44 = -xi(14);

    mat << d11, xi(0), xi(1), xi(2), xi(3), d22, xi(4), xi(5), xi(6), xi(7),
        d33, xi(8), xi(9), xi(10), xi(11), d44;

    return mat;
  }

  // NOTE(hlim): Why 'X'? - I just follow the convetion of GTSAM
  static Vector Vee(const Matrix44& X) {
    Vector vec(15);
    const double x12 = X(0, 0);
    const double x13 = X(1, 1) + x12;
    const double x14 = -X(3, 3);
    vec << X(0, 1), X(0, 2), X(0, 3), X(1, 0), X(1, 2), X(1, 3), X(2, 0),
        X(2, 1), X(2, 3), X(3, 0), X(3, 1), X(3, 2), x12, x13, x14;
    return vec;
  }

  /// Exponential map at identity - create a rotation from canonical coordinates
  static SL4 Expmap(const Vector& xi);

  /// Log map at identity - return the canonical coordinates of this rotation
  static Vector Logmap(const Matrix44& T);

  static Vector Logmap(const SL4& p);

  /// @}
  ///
  Matrix15x15 AdjointMap() {
    Matrix44 H_inv_T = T_.inverse().transpose();
    Matrix16x16 C_H;

    // Kronecker product H ⊗ H^{-T}
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        C_H.block<4, 4>(i * 4, j * 4) = T_(i, j) * H_inv_T;

    return ALG_TO_VEC * C_H * VEC_TO_ALG;
  }

 private:
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
  // Serialization function
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_NVP(T_);
  }
#endif
};  // \class Pose4DoF

template <>
struct traits<SL4> : public internal::LieGroup<SL4> {
  using ChartJacobian = SL4Jacobian;

  static SL4 Identity() { return SL4(); }

  // NOTE(hlim) In PGO, this function is not used
  static SL4 Compose(const SL4& g, const SL4& h, SL4Jacobian H1 = {},
                     SL4Jacobian H2 = {}) {
    return g.compose(h, H1, H2);
  }

  static SL4 Between(const SL4& g, const SL4& h, SL4Jacobian H1 = {},
                     SL4Jacobian H2 = {}) {
    return g.between(h, H1, H2);
  }

  // NOTE(hlim) In PGO, this function is not used
  static SL4 Inverse(const SL4& g, SL4Jacobian H = {}) { return g.inverse(H); }

  static SL4 Expmap(const Vector& xi, SL4Jacobian H = {}) {
    return SL4::Expmap(xi);
  }

  static Vector Logmap(const SL4& g, SL4Jacobian H = {}) {
    return SL4::Logmap(g);
  }
};

template <>
struct traits<const SL4> : public traits<SL4> {};

}  // namespace gtsam
