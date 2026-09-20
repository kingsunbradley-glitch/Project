#include "SimulationConfig.hh"

#include "G4Exception.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

namespace {
G4String Lower(G4String text) {
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return text;
}
}  // namespace

SimulationConfig& SimulationConfig::Instance() {
  static SimulationConfig config;
  return config;
}

SimulationConfig::SimulationConfig() { BuildMessengers(); }
SimulationConfig::~SimulationConfig() = default;

void SimulationConfig::BuildMessengers() {
  geometryMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/geometry/", "Geometry controls");
  geometryMessenger_->DeclarePropertyWithUnit("solenoidRadius", "m", solenoidRadius);
  geometryMessenger_->DeclarePropertyWithUnit("solenoidLength", "m", solenoidLength);
  geometryMessenger_->DeclarePropertyWithUnit("targetDistance", "m", targetDistance);
  geometryMessenger_->DeclarePropertyWithUnit("maxStep", "mm", maxStep);

  fieldMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/field/", "Magnetic field controls");
  fieldMessenger_->DeclarePropertyWithUnit("centerField", "tesla", centerField);

  gasMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/gas/", "Gas controls");
  gasMessenger_->DeclareProperty("enabled", gasEnabled);
  gasMessenger_->DeclarePropertyWithUnit("pressure", "pascal", gasPressure);
  gasMessenger_->DeclarePropertyWithUnit("temperature", "kelvin", gasTemperature);

  reactionMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/reaction/", "Reaction controls");
  reactionMessenger_->DeclareMethod("projectileZA", &SimulationConfig::SetProjectileZA,
                                    "Set projectile Z A");
  reactionMessenger_->DeclareMethod("targetZA", &SimulationConfig::SetTargetZA,
                                    "Set target Z A");
  reactionMessenger_->DeclareMethod("beamEnergyLab", &SimulationConfig::SetBeamEnergy,
                                    "Set lab beam kinetic energy and select beam-energy mode")
      .SetUnit("MeV");
  reactionMessenger_->DeclareMethod("excitationEnergy", &SimulationConfig::SetExcitationEnergy,
                                    "Set compound excitation energy and select E* mode")
      .SetUnit("MeV");
  reactionMessenger_->DeclareProperty("sigma4n", sigma4nMicrobarn)
      .SetGuidance("4n cross section in microbarn; <=0 means normalized-only mode");
  reactionMessenger_->DeclareProperty("sigma5n", sigma5nMicrobarn)
      .SetGuidance("5n cross section in microbarn; <=0 means normalized-only mode");
  reactionMessenger_->DeclareProperty("initialChargeState", initialChargeState)
      .SetGuidance("Initial ionic charge in e; negative means fully stripped at the target");
  reactionMessenger_->DeclareProperty("maxChannelTrials", maxChannelTrials);
  reactionMessenger_->DeclareMethod("mode", &SimulationConfig::SetGeneratorMode,
                                    "reaction or angleScan");

  scanMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/scan/", "Angle scan controls");
  scanMessenger_->DeclarePropertyWithUnit("thetaMin", "deg", thetaMin);
  scanMessenger_->DeclarePropertyWithUnit("thetaMax", "deg", thetaMax);
  scanMessenger_->DeclareProperty("bins", thetaBins);
  scanMessenger_->DeclareProperty("eventsPerBin", eventsPerBin);

  outputMessenger_ =
      std::make_unique<G4GenericMessenger>(this, "/sim/output/", "Output controls");
  outputMessenger_->DeclareProperty("trackSampleCount", trackSampleCount);
  outputMessenger_->DeclarePropertyWithUnit("postExitTrackLength", "m",
                                             postExitTrackLength)
      .SetGuidance("Continue sampled transmitted tracks this far beyond the exit");
  outputMessenger_->DeclareProperty("fileName", outputFile);
  outputMessenger_->DeclareMethod("randomSeed", &SimulationConfig::SetRandomSeed);
}

void SimulationConfig::SetProjectileZA(G4int z, G4int a) {
  projectileZ = z;
  projectileA = a;
}

void SimulationConfig::SetTargetZA(G4int z, G4int a) {
  targetZ = z;
  targetA = a;
}

void SimulationConfig::SetGeneratorMode(const G4String& mode) {
  const auto value = Lower(mode);
  if (value == "reaction") {
    generatorMode = GeneratorMode::Reaction;
  } else if (value == "anglescan" || value == "scan") {
    generatorMode = GeneratorMode::AngleScan;
  } else {
    G4Exception("SimulationConfig::SetGeneratorMode", "BadGeneratorMode", FatalException,
                "Generator mode must be reaction or angleScan.");
  }
}

void SimulationConfig::SetBeamEnergy(G4double value) {
  beamEnergyLab = value;
  kinematicsMode = KinematicsMode::BeamEnergy;
  beamEnergyWasSet = true;
}

void SimulationConfig::SetExcitationEnergy(G4double value) {
  excitationEnergy = value;
  kinematicsMode = KinematicsMode::ExcitationEnergy;
  excitationEnergyWasSet = true;
}

void SimulationConfig::SetRandomSeed(G4long seed) {
  if (seed <= 0) {
    G4Exception("SimulationConfig::SetRandomSeed", "BadSeed", FatalException,
                "Random seed must be positive.");
  }
  randomSeed = seed;
  G4Random::setTheSeed(seed);
}

void SimulationConfig::Validate() const {
  std::ostringstream errors;
  if (!(solenoidRadius > 0.0 && solenoidLength > 0.0 && targetDistance > 0.0)) {
    errors << "R, L and d must all be positive. ";
  }
  if (projectileZ <= 0 || projectileA < projectileZ || targetZ <= 0 || targetA < targetZ) {
    errors << "Projectile/target Z,A are invalid. ";
  }
  if (kinematicsMode == KinematicsMode::BeamEnergy && beamEnergyLab <= 0.0) {
    errors << "Beam kinetic energy must be positive. ";
  }
  if (kinematicsMode == KinematicsMode::ExcitationEnergy && excitationEnergy <= 0.0) {
    errors << "Excitation energy must be positive. ";
  }
  if (thetaBins <= 0 || eventsPerBin <= 0 || thetaMax <= thetaMin) {
    errors << "Angle scan range/bins/events are invalid. ";
  }
  if (gasEnabled && (gasPressure <= 0.0 || gasTemperature <= 0.0)) {
    errors << "Enabled He gas requires positive pressure and temperature. ";
  }
  if (maxStep <= 0.0 || maxChannelTrials <= 0 || trackSampleCount < 0 ||
      postExitTrackLength < 0.0) {
    errors << "Step size, channel trials, or track sample count is invalid. ";
  }
  if (!errors.str().empty()) {
    G4Exception("SimulationConfig::Validate", "InvalidConfiguration", FatalException,
                errors.str().c_str());
  }
}

void SimulationConfig::Print() const {
  G4cout << "\n=== superconducting solenoid configuration ===\n"
         << "R=" << solenoidRadius / m << " m, L=" << solenoidLength / m
         << " m, d=" << targetDistance / m << " m, B0=" << centerField / tesla
         << " T\nHe=" << (gasEnabled ? "on" : "off")
         << ", P=" << gasPressure / pascal << " Pa, T=" << gasTemperature / kelvin
         << " K\nprojectile=(Z=" << projectileZ << ",A=" << projectileA
         << "), target=(Z=" << targetZ << ",A=" << targetA << ")\nmode="
         << (generatorMode == GeneratorMode::Reaction ? "reaction" : "angleScan")
         << ", output=" << outputFile
         << ", post-exit track=" << postExitTrackLength / m << " m"
         << "\n==============================================\n"
         << G4endl;
}
