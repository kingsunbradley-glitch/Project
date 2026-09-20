#pragma once

#include "G4UserEventAction.hh"

class G4Event;
class RunAction;

class EventAction final : public G4UserEventAction {
 public:
  explicit EventAction(RunAction* runAction) : runAction_(runAction) {}
  void EndOfEventAction(const G4Event* event) override;

 private:
  RunAction* runAction_;
};
