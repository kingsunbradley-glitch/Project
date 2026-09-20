#pragma once

#include "G4MagneticField.hh"
#include "globals.hh"

#include <cstddef>
#include <utility>
#include <vector>

class FiniteSolenoidField final : public G4MagneticField {
 public:
  FiniteSolenoidField(G4double radius, G4double length, G4double entranceZ,
                      G4double centerField, std::size_t radialPoints = 81,
                      std::size_t axialPoints = 301, std::size_t currentSlices = 96);

  void GetFieldValue(const G4double point[4], G4double* field) const override;

  std::pair<G4double, G4double> CylindricalField(G4double rho,
                                                 G4double z) const;
  G4double CenterFieldValue() const;
  G4double ZMin() const { return zMin_; }
  G4double ZMax() const { return zMax_; }

 private:
  std::pair<G4double, G4double> UnitSheetField(G4double rho, G4double z) const;
  std::pair<G4double, G4double> LoopField(G4double rho, G4double deltaZ,
                                         G4double current) const;
  std::size_t Index(std::size_t ir, std::size_t iz) const {
    return ir * axialPoints_ + iz;
  }

  G4double radius_;
  G4double length_;
  G4double entranceZ_;
  G4double centerField_;
  G4double zMin_;
  G4double zMax_;
  G4double rhoMax_;
  G4double scale_;
  std::size_t radialPoints_;
  std::size_t axialPoints_;
  std::size_t currentSlices_;
  std::vector<G4double> brGrid_;
  std::vector<G4double> bzGrid_;
};

