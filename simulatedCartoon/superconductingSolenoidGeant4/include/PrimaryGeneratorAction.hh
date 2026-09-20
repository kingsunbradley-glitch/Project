#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <memory>

class G4Event;
class G4ExcitationHandler;
class G4ParticleGun;
class G4ReactionProduct;

class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ReactionProduct* SampleChannel(G4int desiredNeutrons);
  std::unique_ptr<G4ParticleGun> gun_;
  std::unique_ptr<G4ExcitationHandler> excitationHandler_;
};
