#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  const G4LogicalVolume* GetBoreVolume() const { return boreLogical_; }

 private:
  G4LogicalVolume* boreLogical_ = nullptr;
};
