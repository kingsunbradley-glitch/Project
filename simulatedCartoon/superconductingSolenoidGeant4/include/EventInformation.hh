#pragma once

#include "G4VUserEventInformation.hh"
#include "globals.hh"

#include <limits>

enum class FinalStatus : G4int {
  Active = 0,
  Transmitted = 1,
  ApertureLoss = 2,
  NotEntered = 3,
  Stopped = 4,
  Backward = 5,
  WorldExit = 6
};

class EventInformation final : public G4VUserEventInformation {
 public:
  void Print() const override {}

  G4int channelNeutrons = 0;
  G4int residualZ = 0;
  G4int residualA = 0;
  G4double eventWeight = 1.0;
  G4bool absoluteCrossSection = false;
  G4double excitationMeV = 0.0;
  G4double beamEnergyMeV = 0.0;

  G4double initialThetaDeg = 0.0;
  G4double initialPhiDeg = 0.0;
  G4double initialEnergyMeV = 0.0;
  G4double initialCharge = 0.0;
  G4double initialVr = 0.0;
  G4double initialVphi = 0.0;
  G4double initialVz = 0.0;

  G4bool entered = false;
  G4double entryRhoM = std::numeric_limits<G4double>::quiet_NaN();
  G4double entryPhiDeg = std::numeric_limits<G4double>::quiet_NaN();
  G4double entryEnergyMeV = std::numeric_limits<G4double>::quiet_NaN();
  G4double entryCharge = std::numeric_limits<G4double>::quiet_NaN();

  G4double exitRhoM = std::numeric_limits<G4double>::quiet_NaN();
  G4double exitPhiDeg = std::numeric_limits<G4double>::quiet_NaN();
  G4double exitThetaDeg = std::numeric_limits<G4double>::quiet_NaN();
  G4double exitEnergyMeV = std::numeric_limits<G4double>::quiet_NaN();
  G4double exitCharge = std::numeric_limits<G4double>::quiet_NaN();
  FinalStatus status = FinalStatus::Active;
};

