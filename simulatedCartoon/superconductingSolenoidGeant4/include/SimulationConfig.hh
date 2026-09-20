#pragma once

#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <memory>

class G4GenericMessenger;

class SimulationConfig {
 public:
  enum class GeneratorMode { Reaction, AngleScan };
  enum class KinematicsMode { BeamEnergy, ExcitationEnergy };

  static SimulationConfig& Instance();

  SimulationConfig(const SimulationConfig&) = delete;
  SimulationConfig& operator=(const SimulationConfig&) = delete;

  void Validate() const;
  void Print() const;

  void SetProjectileZA(G4int z, G4int a);
  void SetTargetZA(G4int z, G4int a);
  void SetGeneratorMode(const G4String& mode);
  void SetBeamEnergy(G4double value);
  void SetExcitationEnergy(G4double value);
  void SetRandomSeed(G4long seed);

  G4double solenoidRadius = 0.20 * m;
  G4double solenoidLength = 1.00 * m;
  G4double targetDistance = 0.30 * m;
  G4double centerField = 2.50 * tesla;
  G4bool gasEnabled = true;
  G4double gasPressure = 100.0 * pascal;
  G4double gasTemperature = 300.0 * kelvin;

  G4int projectileZ = 18;
  G4int projectileA = 40;
  G4int targetZ = 69;
  G4int targetA = 169;
  G4double beamEnergyLab = 200.0 * MeV;
  G4double excitationEnergy = 0.0;
  KinematicsMode kinematicsMode = KinematicsMode::BeamEnergy;
  G4bool beamEnergyWasSet = false;
  G4bool excitationEnergyWasSet = false;
  G4double sigma4nMicrobarn = 0.0;
  G4double sigma5nMicrobarn = 0.0;
  G4double initialChargeState = -1.0;
  G4int maxChannelTrials = 20000;

  GeneratorMode generatorMode = GeneratorMode::Reaction;
  G4double thetaMin = 0.0;
  G4double thetaMax = 45.0 * deg;
  G4int thetaBins = 60;
  G4int eventsPerBin = 1000;

  G4int trackSampleCount = 50;
  G4double postExitTrackLength = 0.0 * m;
  G4String outputFile = "output/solenoid_transport.root";
  G4long randomSeed = 20260810;
  G4double maxStep = 2.0 * mm;

 private:
  SimulationConfig();
  ~SimulationConfig();

  void BuildMessengers();
  std::unique_ptr<G4GenericMessenger> geometryMessenger_;
  std::unique_ptr<G4GenericMessenger> fieldMessenger_;
  std::unique_ptr<G4GenericMessenger> gasMessenger_;
  std::unique_ptr<G4GenericMessenger> reactionMessenger_;
  std::unique_ptr<G4GenericMessenger> scanMessenger_;
  std::unique_ptr<G4GenericMessenger> outputMessenger_;
};
