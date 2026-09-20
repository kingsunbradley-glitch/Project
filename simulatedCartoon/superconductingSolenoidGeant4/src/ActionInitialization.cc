#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGeneratorAction());
  auto* runAction = new RunAction();
  SetUserAction(runAction);
  SetUserAction(new EventAction(runAction));
  SetUserAction(new SteppingAction(runAction));
}
