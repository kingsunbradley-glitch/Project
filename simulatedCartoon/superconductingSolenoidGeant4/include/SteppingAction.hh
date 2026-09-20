#pragma once

#include "G4UserSteppingAction.hh"

class G4Step;
class RunAction;

class SteppingAction final : public G4UserSteppingAction {
 public:
  explicit SteppingAction(RunAction* runAction) : runAction_(runAction) {}
  void UserSteppingAction(const G4Step* step) override;

 private:
  RunAction* runAction_;
};
