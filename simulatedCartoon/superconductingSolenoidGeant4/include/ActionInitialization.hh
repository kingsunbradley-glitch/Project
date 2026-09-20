#pragma once

#include "G4VUserActionInitialization.hh"

class ActionInitialization final : public G4VUserActionInitialization {
 public:
  void Build() const override;
};
