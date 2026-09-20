#include "DetectorConstruction.hh"

#include "FiniteSolenoidField.hh"
#include "SimulationConfig.hh"

#include "G4Colour.hh"
#include "G4FieldBuilder.hh"
#include "G4FieldParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include <algorithm>

DetectorConstruction::DetectorConstruction() { G4FieldBuilder::Instance(); }

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto& cfg = SimulationConfig::Instance();
  cfg.Validate();

  auto* nist = G4NistManager::Instance();
  auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* transportMaterial = vacuum;
  if (cfg.gasEnabled) {
    const auto density = cfg.gasPressure * (4.002602 * g / mole) /
                         (CLHEP::Avogadro * CLHEP::k_Boltzmann * cfg.gasTemperature);
    transportMaterial = new G4Material("He_transport", 2.0, 4.002602 * g / mole,
                                       density, kStateGas, cfg.gasTemperature,
                                       cfg.gasPressure);
  }

  const auto worldRadius = std::max(2.0 * cfg.solenoidRadius, 0.75 * m);
  const auto downstreamMargin = std::max(0.75 * m, cfg.postExitTrackLength + 0.25 * m);
  const auto worldHalfZ = cfg.targetDistance + cfg.solenoidLength + downstreamMargin;
  auto* worldSolid = new G4Tubs("World", 0.0, worldRadius, worldHalfZ, 0.0, twopi);
  auto* worldLogical = new G4LogicalVolume(worldSolid, transportMaterial, "World");
  auto* worldPhysical = new G4PVPlacement(nullptr, {},
                                          worldLogical, "World", nullptr, false, 0, true);
  worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  worldLogical->SetUserLimits(new G4UserLimits(cfg.maxStep));

  auto* boreSolid = new G4Tubs("SolenoidBore", 0.0, cfg.solenoidRadius,
                               0.5 * cfg.solenoidLength, 0.0, twopi);
  boreLogical_ = new G4LogicalVolume(boreSolid, transportMaterial, "SolenoidBore");
  new G4PVPlacement(nullptr, {0.0, 0.0, cfg.targetDistance + 0.5 * cfg.solenoidLength},
                    boreLogical_, "SolenoidBore", worldLogical, false, 0, true);
  boreLogical_->SetUserLimits(new G4UserLimits(cfg.maxStep));
  auto* boreVis = new G4VisAttributes(G4Colour(0.15, 0.45, 1.0, 0.08));
  boreVis->SetForceSolid(true);
  boreLogical_->SetVisAttributes(boreVis);

  const auto shellThickness = std::max(10.0 * mm, 0.04 * cfg.solenoidRadius);
  auto* shellSolid = new G4Tubs("SolenoidCoil", cfg.solenoidRadius,
                                cfg.solenoidRadius + shellThickness,
                                0.5 * cfg.solenoidLength, 0.0, twopi);
  auto* shellLogical = new G4LogicalVolume(
      shellSolid, nist->FindOrBuildMaterial("G4_Cu"), "SolenoidCoil");
  new G4PVPlacement(nullptr, {0.0, 0.0, cfg.targetDistance + 0.5 * cfg.solenoidLength},
                    shellLogical, "SolenoidCoil", worldLogical, false, 0, true);
  auto* shellVis = new G4VisAttributes(G4Colour(0.85, 0.45, 0.12, 0.24));
  shellVis->SetForceSolid(true);
  shellLogical->SetVisAttributes(shellVis);

  auto* targetSolid = new G4Tubs("TargetMarker", 0.0, 7.0 * mm, 0.10 * mm, 0.0, twopi);
  auto* targetLogical = new G4LogicalVolume(targetSolid, transportMaterial, "TargetMarker");
  new G4PVPlacement(nullptr, {}, targetLogical, "TargetMarker", worldLogical, false, 0, true);
  auto* targetVis = new G4VisAttributes(G4Colour(0.95, 0.15, 0.15));
  targetVis->SetForceSolid(true);
  targetLogical->SetVisAttributes(targetVis);

  return worldPhysical;
}

void DetectorConstruction::ConstructSDandField() {
  const auto& cfg = SimulationConfig::Instance();
  auto* builder = G4FieldBuilder::Instance();
  auto* parameters = builder->GetFieldParameters();
  parameters->SetStepperType(kDormandPrince745);
  parameters->SetMinimumStep(0.01 * mm);
  parameters->SetDeltaChord(0.05 * mm);
  parameters->SetDeltaOneStep(0.005 * mm);
  parameters->SetDeltaIntersection(0.001 * mm);
  parameters->SetMinimumEpsilonStep(1.0e-7);
  parameters->SetMaximumEpsilonStep(1.0e-5);
  builder->SetGlobalField(new FiniteSolenoidField(
      cfg.solenoidRadius, cfg.solenoidLength, cfg.targetDistance, cfg.centerField));
  builder->ConstructFieldSetup();
}
