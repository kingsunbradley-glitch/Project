#include "PrimaryGeneratorAction.hh"

#include "EventInformation.hh"
#include "ReactionKinematics.hh"
#include "SimulationConfig.hh"

#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ReactionProduct.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <memory>

namespace {
FusionKinematics GetKinematics() {
  const auto& cfg = SimulationConfig::Instance();
  if (cfg.beamEnergyWasSet && cfg.excitationEnergyWasSet) {
    const auto fromBeam = ReactionKinematics::FromBeamEnergy(
        cfg.projectileZ, cfg.projectileA, cfg.targetZ, cfg.targetA,
        cfg.beamEnergyLab);
    if (std::abs(fromBeam.excitationEnergy - cfg.excitationEnergy) > 0.5 * MeV) {
      G4Exception("PrimaryGeneratorAction::GetKinematics", "InconsistentEnergyInput",
                  FatalException,
                  "beamEnergyLab and excitationEnergy imply E* values differing by more than 0.5 MeV.");
    }
  }
  if (cfg.kinematicsMode == SimulationConfig::KinematicsMode::ExcitationEnergy) {
    return ReactionKinematics::FromExcitationEnergy(
        cfg.projectileZ, cfg.projectileA, cfg.targetZ, cfg.targetA,
        cfg.excitationEnergy);
  }
  return ReactionKinematics::FromBeamEnergy(
      cfg.projectileZ, cfg.projectileA, cfg.targetZ, cfg.targetA,
      cfg.beamEnergyLab);
}

void DeleteProducts(G4ReactionProductVector* products, G4ReactionProduct* keep) {
  if (products == nullptr) return;
  for (auto* product : *products) {
    if (product != keep) delete product;
  }
  delete products;
}
}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : gun_(std::make_unique<G4ParticleGun>(1)),
      excitationHandler_(std::make_unique<G4ExcitationHandler>()) {
  excitationHandler_->Initialise();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

G4ReactionProduct* PrimaryGeneratorAction::SampleChannel(G4int desiredNeutrons) {
  const auto& cfg = SimulationConfig::Instance();
  const auto kin = GetKinematics();
  const auto desiredA = kin.compoundA - desiredNeutrons;
  for (G4int trial = 0; trial < cfg.maxChannelTrials; ++trial) {
    const G4Fragment compound(kin.compoundA, kin.compoundZ, kin.compoundFourMomentum);
    auto* products = excitationHandler_->BreakItUp(compound);
    G4ReactionProduct* selected = nullptr;
    if (products != nullptr) {
      for (auto* product : *products) {
        const auto* definition = product->GetDefinition();
        if (definition != nullptr && definition->GetAtomicNumber() == kin.compoundZ &&
            definition->GetAtomicMass() == desiredA) {
          selected = product;
          break;
        }
      }
    }
    if (selected != nullptr) {
      DeleteProducts(products, selected);
      return selected;
    }
    DeleteProducts(products, nullptr);
  }
  G4Exception("PrimaryGeneratorAction::SampleChannel", "ChannelSamplingFailed",
              FatalException,
              "Could not obtain the requested xn residue. Increase maxChannelTrials or use an E* where that channel is populated.");
  return nullptr;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  auto& cfg = SimulationConfig::Instance();
  const auto kin = GetKinematics();
  const G4int desiredNeutrons = (event->GetEventID() % 2 == 0) ? 4 : 5;
  std::unique_ptr<G4ReactionProduct> product(SampleChannel(desiredNeutrons));
  auto momentum = product->GetMomentum();
  auto direction = momentum.unit();

  if (cfg.generatorMode == SimulationConfig::GeneratorMode::AngleScan) {
    const auto bin = std::min(event->GetEventID() / cfg.eventsPerBin, cfg.thetaBins - 1);
    const auto theta = cfg.thetaMin + (static_cast<G4double>(bin) + 0.5) *
                                      (cfg.thetaMax - cfg.thetaMin) / cfg.thetaBins;
    const auto phi = twopi * G4UniformRand();
    direction = {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                 std::cos(theta)};
    momentum = momentum.mag() * direction;
  }

  auto* ion = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(
      kin.compoundZ, kin.compoundA - desiredNeutrons, 0.0);
  if (ion == nullptr) {
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "IonNotFound",
                FatalException, "Geant4 ion table could not create the evaporation residue.");
  }
  const auto charge = cfg.initialChargeState < 0.0 ? kin.compoundZ :
                                                    std::clamp(cfg.initialChargeState, 0.0,
                                                               static_cast<G4double>(kin.compoundZ));
  gun_->SetParticleDefinition(ion);
  gun_->SetParticlePosition({0.0, 0.0, 0.0});
  gun_->SetParticleMomentumDirection(direction);
  gun_->SetParticleEnergy(product->GetKineticEnergy());
  gun_->SetParticleCharge(charge * eplus);
  gun_->GeneratePrimaryVertex(event);

  auto* info = new EventInformation();
  info->channelNeutrons = desiredNeutrons;
  info->residualZ = kin.compoundZ;
  info->residualA = kin.compoundA - desiredNeutrons;
  const auto sigma = desiredNeutrons == 4 ? cfg.sigma4nMicrobarn : cfg.sigma5nMicrobarn;
  info->absoluteCrossSection = cfg.sigma4nMicrobarn > 0.0 && cfg.sigma5nMicrobarn > 0.0;
  info->eventWeight = info->absoluteCrossSection ? sigma : 1.0;
  info->excitationMeV = kin.excitationEnergy / MeV;
  info->beamEnergyMeV = kin.beamKineticEnergy / MeV;
  info->initialThetaDeg = direction.theta() / deg;
  info->initialPhiDeg = direction.phi() / deg;
  info->initialEnergyMeV = product->GetKineticEnergy() / MeV;
  info->initialCharge = charge;
  const auto totalEnergy = product->GetTotalEnergy();
  const auto speed = totalEnergy > 0.0 ? c_light * momentum.mag() / totalEnergy : 0.0;
  info->initialVr = speed * std::sin(direction.theta()) / (m / s);
  info->initialVphi = 0.0;
  info->initialVz = speed * direction.z() / (m / s);
  event->SetUserInformation(info);
}
