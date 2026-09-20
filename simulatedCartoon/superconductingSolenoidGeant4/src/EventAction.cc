#include "EventAction.hh"

#include "EventInformation.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"

void EventAction::EndOfEventAction(const G4Event* event) {
  auto* info = dynamic_cast<EventInformation*>(event->GetUserInformation());
  if (info == nullptr) return;
  if (info->status == FinalStatus::Active) {
    info->status = info->entered ? FinalStatus::WorldExit : FinalStatus::NotEntered;
  }
  auto* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(0, 0, event->GetEventID());
  analysis->FillNtupleIColumn(0, 1, info->channelNeutrons);
  analysis->FillNtupleIColumn(0, 2, info->residualZ);
  analysis->FillNtupleIColumn(0, 3, info->residualA);
  analysis->FillNtupleDColumn(0, 4, info->eventWeight);
  analysis->FillNtupleIColumn(0, 5, info->absoluteCrossSection ? 1 : 0);
  analysis->FillNtupleDColumn(0, 6, info->beamEnergyMeV);
  analysis->FillNtupleDColumn(0, 7, info->excitationMeV);
  analysis->FillNtupleDColumn(0, 8, info->initialThetaDeg);
  analysis->FillNtupleDColumn(0, 9, info->initialPhiDeg);
  analysis->FillNtupleDColumn(0, 10, info->initialEnergyMeV);
  analysis->FillNtupleDColumn(0, 11, info->initialCharge);
  analysis->FillNtupleDColumn(0, 12, info->initialVr);
  analysis->FillNtupleDColumn(0, 13, info->initialVphi);
  analysis->FillNtupleDColumn(0, 14, info->initialVz);
  analysis->FillNtupleIColumn(0, 15, info->entered ? 1 : 0);
  analysis->FillNtupleDColumn(0, 16, info->entryRhoM);
  analysis->FillNtupleDColumn(0, 17, info->entryPhiDeg);
  analysis->FillNtupleDColumn(0, 18, info->entryEnergyMeV);
  analysis->FillNtupleDColumn(0, 19, info->entryCharge);
  analysis->FillNtupleDColumn(0, 20, info->exitRhoM);
  analysis->FillNtupleDColumn(0, 21, info->exitPhiDeg);
  analysis->FillNtupleDColumn(0, 22, info->exitThetaDeg);
  analysis->FillNtupleDColumn(0, 23, info->exitEnergyMeV);
  analysis->FillNtupleDColumn(0, 24, info->exitCharge);
  analysis->FillNtupleIColumn(0, 25, static_cast<G4int>(info->status));
  analysis->AddNtupleRow(0);
  runAction_->RecordEvent(event->GetEventID(), *info);
}
