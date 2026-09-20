#include "ReactionKinematics.hh"

#include "G4Exception.hh"
#include "G4NucleiProperties.hh"

#include <algorithm>
#include <cmath>

namespace {
FusionKinematics Base(G4int projectileZ, G4int projectileA, G4int targetZ,
                      G4int targetA) {
  FusionKinematics result;
  result.compoundZ = projectileZ + targetZ;
  result.compoundA = projectileA + targetA;
  result.projectileMass = G4NucleiProperties::GetNuclearMass(projectileA, projectileZ);
  result.targetMass = G4NucleiProperties::GetNuclearMass(targetA, targetZ);
  result.compoundMass =
      G4NucleiProperties::GetNuclearMass(result.compoundA, result.compoundZ);
  return result;
}

void Check(const FusionKinematics& value) {
  if (value.beamKineticEnergy <= 0.0 || value.excitationEnergy < 0.0 ||
      value.projectileMass <= 0.0 || value.targetMass <= 0.0 || value.compoundMass <= 0.0) {
    G4Exception("ReactionKinematics", "InvalidFusionKinematics", FatalException,
                "The selected projectile/target and energy do not form a valid compound state.");
  }
}
}  // namespace

FusionKinematics ReactionKinematics::FromBeamEnergy(G4int projectileZ, G4int projectileA,
                                                    G4int targetZ, G4int targetA,
                                                    G4double beamKineticEnergy) {
  auto result = Base(projectileZ, projectileA, targetZ, targetA);
  const G4double projectileTotal = result.projectileMass + beamKineticEnergy;
  const G4double projectileMomentum =
      std::sqrt(std::max(0.0, projectileTotal * projectileTotal -
                                  result.projectileMass * result.projectileMass));
  const G4double totalEnergy = projectileTotal + result.targetMass;
  const G4double invariantMass =
      std::sqrt(std::max(0.0, totalEnergy * totalEnergy -
                                  projectileMomentum * projectileMomentum));
  result.beamKineticEnergy = beamKineticEnergy;
  result.excitationEnergy = invariantMass - result.compoundMass;
  result.compoundFourMomentum =
      G4LorentzVector(0.0, 0.0, projectileMomentum, totalEnergy);
  Check(result);
  return result;
}

FusionKinematics ReactionKinematics::FromExcitationEnergy(
    G4int projectileZ, G4int projectileA, G4int targetZ, G4int targetA,
    G4double excitationEnergy) {
  auto result = Base(projectileZ, projectileA, targetZ, targetA);
  const G4double invariantMass = result.compoundMass + excitationEnergy;
  const G4double projectileTotal =
      (invariantMass * invariantMass - result.projectileMass * result.projectileMass -
       result.targetMass * result.targetMass) /
      (2.0 * result.targetMass);
  const G4double projectileMomentum =
      std::sqrt(std::max(0.0, projectileTotal * projectileTotal -
                                  result.projectileMass * result.projectileMass));
  result.beamKineticEnergy = projectileTotal - result.projectileMass;
  result.excitationEnergy = excitationEnergy;
  result.compoundFourMomentum = G4LorentzVector(
      0.0, 0.0, projectileMomentum, projectileTotal + result.targetMass);
  Check(result);
  return result;
}

