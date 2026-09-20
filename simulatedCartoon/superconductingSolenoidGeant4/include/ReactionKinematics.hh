#pragma once

#include "G4LorentzVector.hh"
#include "globals.hh"

struct FusionKinematics {
  G4int compoundZ = 0;
  G4int compoundA = 0;
  G4double projectileMass = 0.0;
  G4double targetMass = 0.0;
  G4double compoundMass = 0.0;
  G4double beamKineticEnergy = 0.0;
  G4double excitationEnergy = 0.0;
  G4LorentzVector compoundFourMomentum;
};

class ReactionKinematics {
 public:
  static FusionKinematics FromBeamEnergy(G4int projectileZ, G4int projectileA,
                                         G4int targetZ, G4int targetA,
                                         G4double beamKineticEnergy);

  static FusionKinematics FromExcitationEnergy(G4int projectileZ, G4int projectileA,
                                               G4int targetZ, G4int targetA,
                                               G4double excitationEnergy);
};

