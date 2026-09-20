#include "SteppingAction.hh"

#include "EventInformation.hh"
#include "RunAction.hh"
#include "SimulationConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>

namespace {
G4ThreeVector InterpolateAtZ(const G4StepPoint* pre, const G4StepPoint* post,
                            G4double planeZ) {
  const auto z0 = pre->GetPosition().z();
  const auto dz = post->GetPosition().z() - z0;
  const auto fraction = std::abs(dz) > 1.0e-20 * mm
                            ? std::clamp((planeZ - z0) / dz, 0.0, 1.0)
                            : 0.0;
  return pre->GetPosition() + fraction * (post->GetPosition() - pre->GetPosition());
}

void CylindricalVelocity(const G4Track* track, G4double& vrho, G4double& vphi,
                         G4double& vz) {
  const auto& position = track->GetPosition();
  const auto rho = std::hypot(position.x(), position.y());
  const auto velocity = track->GetVelocity() * track->GetMomentumDirection();
  if (rho > 1.0e-15 * m) {
    vrho = (velocity.x() * position.x() + velocity.y() * position.y()) / rho;
    vphi = (-velocity.x() * position.y() + velocity.y() * position.x()) / rho;
  } else {
    vrho = velocity.x();
    vphi = velocity.y();
  }
  vz = velocity.z();
}
}  // namespace

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto* track = step->GetTrack();
  if (track->GetTrackID() != 1 || track->GetParentID() != 0) return;
  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  auto* info = dynamic_cast<EventInformation*>(event->GetUserInformation());
  if (info == nullptr) return;
  const auto& cfg = SimulationConfig::Instance();
  const auto sampled = event->GetEventID() < cfg.trackSampleCount;
  const auto extendingTransmittedTrack =
      info->status == FinalStatus::Transmitted && sampled &&
      cfg.postExitTrackLength > 0.0;
  if (info->status != FinalStatus::Active && !extendingTransmittedTrack) return;
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto& prePosition = pre->GetPosition();
  const auto& postPosition = post->GetPosition();
  const auto postRho = std::hypot(postPosition.x(), postPosition.y());
  const auto charge = track->GetDynamicParticle()->GetCharge() / eplus;

  if (sampled) {
    G4double vrho = 0.0;
    G4double vphi = 0.0;
    G4double vz = 0.0;
    CylindricalVelocity(track, vrho, vphi, vz);
    auto* analysis = G4AnalysisManager::Instance();
    analysis->FillNtupleIColumn(1, 0, event->GetEventID());
    analysis->FillNtupleIColumn(1, 1, track->GetCurrentStepNumber());
    analysis->FillNtupleDColumn(1, 2, postRho / m);
    analysis->FillNtupleDColumn(1, 3, postPosition.phi() / deg);
    analysis->FillNtupleDColumn(1, 4, postPosition.z() / m);
    analysis->FillNtupleDColumn(1, 5, vrho / (m / s));
    analysis->FillNtupleDColumn(1, 6, vphi / (m / s));
    analysis->FillNtupleDColumn(1, 7, vz / (m / s));
    analysis->FillNtupleDColumn(1, 8, track->GetKineticEnergy() / MeV);
    analysis->FillNtupleDColumn(1, 9, charge);
    analysis->FillNtupleDColumn(1, 10, track->GetGlobalTime() / ns);
    analysis->AddNtupleRow(1);
    runAction_->RecordTrack(event->GetEventID(), track->GetCurrentStepNumber(), postRho / m,
                            postPosition.phi() / deg, postPosition.z() / m,
                            vrho / (m / s), vphi / (m / s), vz / (m / s),
                            track->GetKineticEnergy() / MeV, charge,
                            track->GetGlobalTime() / ns);
  }

  const auto entryZ = cfg.targetDistance;
  const auto exitZ = cfg.targetDistance + cfg.solenoidLength;
  if (extendingTransmittedTrack) {
    if (postPosition.z() >= exitZ + cfg.postExitTrackLength ||
        post->GetKineticEnergy() <= 1.0 * eV ||
        post->GetMomentumDirection().z() <= 0.0) {
      track->SetTrackStatus(fStopAndKill);
    }
    return;
  }
  if (!info->entered && prePosition.z() < entryZ && postPosition.z() >= entryZ) {
    const auto crossing = InterpolateAtZ(pre, post, entryZ);
    const auto rho = std::hypot(crossing.x(), crossing.y());
    if (rho >= cfg.solenoidRadius) {
      info->status = FinalStatus::NotEntered;
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    info->entered = true;
    info->entryRhoM = rho / m;
    info->entryPhiDeg = crossing.phi() / deg;
    info->entryEnergyMeV = post->GetKineticEnergy() / MeV;
    info->entryCharge = charge;
  }

  if (info->entered && prePosition.z() < exitZ && postPosition.z() >= exitZ) {
    const auto crossing = InterpolateAtZ(pre, post, exitZ);
    const auto rho = std::hypot(crossing.x(), crossing.y());
    if (rho < cfg.solenoidRadius) {
      info->exitRhoM = rho / m;
      info->exitPhiDeg = crossing.phi() / deg;
      info->exitThetaDeg = post->GetMomentumDirection().theta() / deg;
      info->exitEnergyMeV = post->GetKineticEnergy() / MeV;
      info->exitCharge = charge;
      info->status = FinalStatus::Transmitted;
      if (sampled && cfg.postExitTrackLength > 0.0) return;
    } else {
      info->status = FinalStatus::ApertureLoss;
    }
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  if (info->entered && postPosition.z() >= entryZ && postPosition.z() <= exitZ &&
      postRho >= cfg.solenoidRadius) {
    info->status = FinalStatus::ApertureLoss;
    track->SetTrackStatus(fStopAndKill);
    return;
  }
  if (post->GetKineticEnergy() <= 1.0 * eV) {
    info->status = FinalStatus::Stopped;
    track->SetTrackStatus(fStopAndKill);
    return;
  }
  if (post->GetMomentumDirection().z() <= 0.0) {
    info->status = FinalStatus::Backward;
    track->SetTrackStatus(fStopAndKill);
  }
}
