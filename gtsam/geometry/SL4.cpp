/**
 * @file  SL4.cpp
 * @brief Projective Special Linear Group (PSL(4, R)) Pose
 * @author: Hyungtae Lim
 */

#include <gtsam/geometry/SL4.h>

using namespace std;

namespace gtsam {

SL4::SL4(const Matrix44& pose) {
  double det = pose.determinant();
  if (det <= 0.0) {
    throw std::runtime_error(
        "Matrix determinant must be positive for SL(4) normalization. Got det "
        "= " +
        std::to_string(det));
  }

  T_ = pose / std::pow(det, 1.0 / 4.0);
}

SL4::SL4(const SL4& pose) : T_(pose.T_) {}

/* ************************************************************************* */
void SL4::print(const std::string& s) const { cout << s << T_ << "\n"; }

/* ************************************************************************* */
bool SL4::equals(const SL4& sl4, double tol) const {
  return T_.isApprox(sl4.T_, tol);
}

/* ************************************************************************* */
// NOTE(hlim): In PGO, this function is not used
SL4 SL4::inverse(SL4Jacobian H1) const {
  if (!H1) {
    return SL4(T_.inverse());
  }

  // TODO(hlim): Might not affect the PGO quality at all,
  // but should be implemented for the complete implementation
  SL4 result(T_.inverse());
  throw std::runtime_error("H matrix is not implemented.");
  return result;
}

/* ************************************************************************* */
// NOTE(hlim): In PGO, this function is not used
// SL4 SL4::compose(const SL4& sl4) const {
//   return SL4(T_ * sl4.T_);
// }

/* ************************************************************************* */
// NOTE(hlim): In PGO, this function is not used
SL4 SL4::compose(const SL4& sl4, SL4Jacobian H1, SL4Jacobian H2) const {
  if (!H1 && !H2) return SL4(T_ * sl4.T_);

  // TODO(hlim): Might not affect the PGO quality at all,
  // but should be implemented for the complete implementation
  SL4 result(T_ * sl4.T_);
  throw std::runtime_error("H matrix is not implemented.");
  return result;
}

/* ************************************************************************* */
SL4 SL4::between(const SL4& sl4, SL4Jacobian H1, SL4Jacobian H2) const {
  if (!H1 && !H2) return SL4(T_.inverse() * sl4.T_);

  SL4 result(T_.inverse() * sl4.T_);
  if (H1) {
    *H1 = -result.inverse().AdjointMap();
  }
  if (H2) *H2 = I_15x15;
  return result;
}

SL4 SL4::retract(const Vector& v, SL4Jacobian Horigin, SL4Jacobian Hv) const {
  assert(v.size() == 15);

  SL4 retracted_pose = SL4(T_ * (I_4x4 + Hat(v)));

  if (Horigin) {
    // TODO(hlim) Should be implemented
    throw std::runtime_error("H matrix is not implemented.");
    *Horigin = I_15x15;
  }

  if (Hv) {
    // TODO(hlim) This matrix should be double checked
    throw std::runtime_error("H matrix is not implemented.");
    *Hv = I_15x15;
  }

  return retracted_pose;
}

Vector SL4::localCoordinates(const SL4& sl4, SL4Jacobian Horigin,
                             SL4Jacobian Hp2) const {
  OptionalJacobian<3, 3>::Jacobian H3x3_1, H3x3_2;
  Vector result = SL4::Logmap(T_.inverse() * sl4.T_);

  if (Horigin) {
    // TODO(hlim) This matrix should be double checked
    throw std::runtime_error("H matrix is not implemented.");
    *Horigin = -I_15x15;
  }

  if (Hp2) {
    // TODO(hlim) This matrix should be double checked
    throw std::runtime_error("H matrix is not implemented.");
    *Hp2 = I_15x15;
  }

  return result;
}

/* ************************************************************************* */
SL4 SL4::Expmap(const Vector& xi) {
  assert(xi.size() == 15);
  const auto& mat = Hat(xi);

  // NOTE(hlim):
  // The cost of the computation is approximately 20n^3 for matrices of size n.
  // The number 20 depends weakly on the norm of the matrix.
  // See
  // https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html

  // TODO(hlim): Approximate exp function? But it introduces non-negligible
  // numerical error. return SL4(Matrix44::Identity() + mat + 0.5 * mat * mat);

  return SL4(mat.exp());
}

/* ************************************************************************* */
Vector SL4::Logmap(const Matrix44& T) {
  const Matrix44& mat = T.log();

  return Vee(mat);
}

Vector SL4::Logmap(const SL4& p) { return Logmap(p.T_); }
/* ************************************************************************* */

}  // namespace gtsam
