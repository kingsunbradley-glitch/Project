#pragma once

#include "G4UserRunAction.hh"

#include <fstream>
#include <filesystem>
#include <vector>

class EventInformation;
class G4Run;

class RunAction final : public G4UserRunAction {
 public:
  RunAction();
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;
  void RecordEvent(G4int eventId, const EventInformation& info);
  void RecordTrack(G4int eventId, G4int step, G4double rhoM, G4double phiDeg,
                   G4double zM, G4double vrho, G4double vphi, G4double vz,
                   G4double energyMeV, G4double charge, G4double timeNs);

 private:
  struct BinCounter {
    G4int generated = 0;
    G4int entered = 0;
    G4int transmitted = 0;
    G4double generatedWeight = 0.0;
    G4double transmittedWeight = 0.0;
  };
  std::vector<BinCounter> bins_;
  G4int terminalCounts_[7] = {};
  G4double channelGeneratedWeight_[2] = {};
  G4double channelTransmittedWeight_[2] = {};
  G4double maxObservedTransmittedThetaDeg_ = -1.0;
  std::ofstream eventCsv_;
  std::ofstream trackCsv_;
  std::filesystem::path outputBase_;
};
