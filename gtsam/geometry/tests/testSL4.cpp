/* ************************************************************************* */
/**
 * @file    testSL4.cpp
 * @brief   Unit tests for SL4 manifold
 * @author  Hyungtae Lim
 */
/* ************************************************************************* */

#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/TestableAssertions.h>
#include <gtsam/base/lieProxies.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/base/testLie.h>
#include <gtsam/geometry/SL4.h>

using namespace std;
using namespace gtsam;

GTSAM_CONCEPT_TESTABLE_INST(SL4)
GTSAM_CONCEPT_MATRIX_LIE_GROUP_INST(SL4)

// Common static variables for tests
static const Vector xi0 = (Vector(15) << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                           0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
                              .finished();

static const Vector xi1 = (Vector(15) << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                           0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
                              .finished();

static const Vector xi2 = (Vector(15) << 0.5, 0.4, 0.9, 0.2, 0.3, 0.8, 0.7, 0.6,
                           0.1, 1.0, 1.2, 1.1, 1.5, 1.3, 1.4)
                              .finished();

// Create a SL4
const Matrix4 T_matrix =
    (Matrix4() << 1, 0, 0, 1, 0, 1, 0, 2, 0, 0, 1, 3, 0, 0, 0, 1).finished();
const SL4 T1(T_matrix);

const Matrix4 T_matrix2 =
    (Matrix4() << 1, 0, 0, 4, 0, 1, 0, 5, 0, 0, 1, 6, 0, 0, 0, 1).finished();
const SL4 T2(T_matrix2);

/* ************************************************************************* */
TEST(SL4, Constructor) {
  SL4 sl4_default;
  EXPECT(assert_equal(SL4(), sl4_default));
}

/* ************************************************************************* */
TEST(SL4, Expmap) {
  SL4 expected(SL4::Expmap(xi0));
  EXPECT(assert_equal(expected, SL4::Expmap(xi0), 1e-8));
}

/* ************************************************************************* */
TEST(SL4, Logmap) {
  SL4 sl4_exp(SL4::Expmap(xi0));
  EXPECT(assert_equal(xi0, SL4::Logmap(sl4_exp), 1e-8));
}

/* ************************************************************************* */
TEST(SL4, Retract) {
  SL4 sl4_retracted = T1.retract(xi0);
  SL4 expected(T1.matrix() * (I_4x4 + SL4::Hat(xi0)));
  EXPECT(assert_equal(expected, sl4_retracted, 1e-8));
}

/* ************************************************************************* */
TEST(SL4, LocalCoordinates) {
  SL4 sl4_retracted = T1.retract(xi0);
  EXPECT(assert_equal(xi0, T1.localCoordinates(sl4_retracted), 1e-8));
}

/* ************************************************************************* */
TEST(SL4, Between) {
  SL4 sl4_1 = SL4::Expmap(xi1);
  SL4 sl4_2 = SL4::Expmap(xi2);

  Matrix H1, H2;
  SL4 expected(sl4_1.matrix().inverse() * sl4_2.matrix());
  SL4 actual = sl4_1.between(sl4_2, H1, H2);
  EXPECT(assert_equal(expected, actual, 1e-8));

  Matrix numericalH1 =
      numericalDerivative21(testing::between<SL4>, sl4_1, sl4_2);
  EXPECT(assert_equal(numericalH1, H1, 5e-3));

  Matrix numericalH2 =
      numericalDerivative22(testing::between<SL4>, sl4_1, sl4_2);
  EXPECT(assert_equal(numericalH2, H2, 1e-5));
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */