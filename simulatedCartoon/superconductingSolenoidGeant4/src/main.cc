#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "SimulationConfig.hh"

#include "FTFP_BERT.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronicParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4RunManagerFactory.hh"
#include "Randomize.hh"

#include <memory>

int main(int argc, char** argv) {
  auto& cfg = SimulationConfig::Instance();
  G4Random::setTheSeed(cfg.randomSeed);
  G4EmParameters::Instance()->SetVerbose(0);
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  G4NuclearLevelData::GetInstance()->GetParameters()->SetVerbose(0);
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new DetectorConstruction());
  auto* physics = new FTFP_BERT();
  physics->SetVerboseLevel(0);
  physics->ReplacePhysics(new G4EmStandardPhysics_option4());
  physics->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new ActionInitialization());

  auto visManager = std::make_unique<G4VisExecutive>("quiet");
  visManager->Initialize();
  auto* uiManager = G4UImanager::GetUIpointer();
  G4int status = 0;
  if (argc > 1) {
    status = uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
  } else {
    auto ui = std::make_unique<G4UIExecutive>(argc, argv);
    status = uiManager->ApplyCommand("/control/execute macros/vis.mac");
    if (status == 0) ui->SessionStart();
  }
  delete runManager;
  return status == 0 ? 0 : 1;
}
