#include "FiniteSolenoidField.hh"

#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

FiniteSolenoidField::FiniteSolenoidField(G4double radius, G4double length,
                                         G4double entranceZ, G4double centerField,
                                         std::size_t radialPoints,
                                         std::size_t axialPoints,
                                         std::size_t currentSlices)
    : radius_(radius),
      length_(length),
      entranceZ_(entranceZ),
      centerField_(centerField),
      zMin_(std::min(0.0, entranceZ - 2.5 * radius)),
      zMax_(entranceZ + length + 2.5 * radius),
      rhoMax_(0.995 * radius),
      scale_(1.0),
      radialPoints_(std::max<std::size_t>(radialPoints, 3)),
      axialPoints_(std::max<std::size_t>(axialPoints, 3)),
      currentSlices_(std::max<std::size_t>(currentSlices, 8)),
      brGrid_(radialPoints_ * axialPoints_),
      bzGrid_(radialPoints_ * axialPoints_) {
  if (radius <= 0.0 || length <= 0.0) {
    G4Exception("FiniteSolenoidField", "InvalidGeometry", FatalException,
                "Finite solenoid radius and length must be positive.");
  }

  const auto unitCenter = UnitSheetField(0.0, entranceZ_ + 0.5 * length_).second;
  if (std::abs(unitCenter) < 1.0e-30 * tesla) {
    G4Exception("FiniteSolenoidField", "ZeroNormalization", FatalException,
                "Could not normalize the finite-solenoid field.");
  }
  scale_ = centerField_ / unitCenter;

  for (std::size_t ir = 0; ir < radialPoints_; ++ir) {
    const G4double rho = rhoMax_ * static_cast<G4double>(ir) /
                         static_cast<G4double>(radialPoints_ - 1);
    for (std::size_t iz = 0; iz < axialPoints_; ++iz) {
      const G4double z = zMin_ + (zMax_ - zMin_) * static_cast<G4double>(iz) /
                                    static_cast<G4double>(axialPoints_ - 1);
      const auto [br, bz] = UnitSheetField(rho, z);
      brGrid_[Index(ir, iz)] = scale_ * br;
      bzGrid_[Index(ir, iz)] = scale_ * bz;
    }
  }
}

std::pair<G4double, G4double> FiniteSolenoidField::LoopField(
    G4double rho, G4double deltaZ, G4double current) const {
  const G4double a = radius_;
  const G4double mu0 = 4.0 * pi * 1.0e-7 * tesla * m / ampere;
  if (rho < 1.0e-10 * m) {
    const G4double denom = std::pow(a * a + deltaZ * deltaZ, 1.5);
    return {0.0, mu0 * current * a * a / (2.0 * denom)};
  }

  const G4double sum2 = (a + rho) * (a + rho) + deltaZ * deltaZ;
  const G4double diff2 = (a - rho) * (a - rho) + deltaZ * deltaZ;
  const G4double k2 = std::clamp(4.0 * a * rho / sum2, 0.0, 1.0 - 1.0e-14);
  const G4double k = std::sqrt(k2);
  const G4double K = std::comp_ellint_1(k);
  const G4double E = std::comp_ellint_2(k);
  const G4double common = mu0 * current / (2.0 * pi * std::sqrt(sum2));
  const G4double br =
      common * deltaZ / rho *
      (-K + (a * a + rho * rho + deltaZ * deltaZ) / diff2 * E);
  const G4double bz =
      common * (K + (a * a - rho * rho - deltaZ * deltaZ) / diff2 * E);
  return {br, bz};
}

std::pair<G4double, G4double> FiniteSolenoidField::UnitSheetField(
    G4double rho, G4double z) const {
  const G4double dz = length_ / static_cast<G4double>(currentSlices_);
  G4double br = 0.0;
  G4double bz = 0.0;
  for (std::size_t i = 0; i < currentSlices_; ++i) {
    const G4double loopZ = entranceZ_ + (static_cast<G4double>(i) + 0.5) * dz;
    const auto [loopBr, loopBz] = LoopField(rho, z - loopZ, dz * ampere / m);
    br += loopBr;
    bz += loopBz;
  }
  return {br, bz};
}

std::pair<G4double, G4double> FiniteSolenoidField::CylindricalField(
    G4double rho, G4double z) const {
  const G4double r = std::clamp(std::abs(rho), 0.0, rhoMax_);
  if (z < zMin_ || z > zMax_) {
    const auto [br, bz] = UnitSheetField(r, z);
    return {scale_ * br, scale_ * bz};
  }
  const G4double zz = z;

  const G4double fr = r / rhoMax_ * static_cast<G4double>(radialPoints_ - 1);
  const G4double fz = (zz - zMin_) / (zMax_ - zMin_) *
                      static_cast<G4double>(axialPoints_ - 1);
  const auto ir0 = std::min<std::size_t>(static_cast<std::size_t>(fr), radialPoints_ - 2);
  const auto iz0 = std::min<std::size_t>(static_cast<std::size_t>(fz), axialPoints_ - 2);
  const G4double tr = fr - static_cast<G4double>(ir0);
  const G4double tz = fz - static_cast<G4double>(iz0);

  const auto interpolate = [&](const std::vector<G4double>& grid) {
    const G4double v00 = grid[Index(ir0, iz0)];
    const G4double v10 = grid[Index(ir0 + 1, iz0)];
    const G4double v01 = grid[Index(ir0, iz0 + 1)];
    const G4double v11 = grid[Index(ir0 + 1, iz0 + 1)];
    return (1.0 - tr) * (1.0 - tz) * v00 + tr * (1.0 - tz) * v10 +
           (1.0 - tr) * tz * v01 + tr * tz * v11;
  };
  return {interpolate(brGrid_), interpolate(bzGrid_)};
}

void FiniteSolenoidField::GetFieldValue(const G4double point[4],
                                        G4double* field) const {
  const G4double x = point[0];
  const G4double y = point[1];
  const G4double rho = std::hypot(x, y);
  const auto [br, bz] = CylindricalField(rho, point[2]);
  if (rho > 1.0e-12 * m) {
    field[0] = br * x / rho;
    field[1] = br * y / rho;
  } else {
    field[0] = 0.0;
    field[1] = 0.0;
  }
  field[2] = bz;
}

G4double FiniteSolenoidField::CenterFieldValue() const {
  return CylindricalField(0.0, entranceZ_ + 0.5 * length_).second;
}
