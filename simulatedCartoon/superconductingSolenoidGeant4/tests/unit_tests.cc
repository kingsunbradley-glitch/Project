#include "FiniteSolenoidField.hh"
#include "ReactionKinematics.hh"

#include "G4ExcitationHandler.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4Fragment.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ReactionProduct.hh"
#include "G4SystemOfUnits.hh"
#include "G4ShortLivedConstructor.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void Require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

bool Close(double a, double b, double relative, double absolute = 0.0) {
  return std::abs(a - b) <= std::max(absolute, relative * std::max(std::abs(a), std::abs(b)));
}
}  // namespace

int main() {
  try {
    const auto fromBeam = ReactionKinematics::FromBeamEnergy(18, 40, 69, 169, 200.0 * MeV);
    Require(fromBeam.compoundZ == 87 && fromBeam.compoundA == 209,
            "40Ar+169Tm must form 209Fr*.");
    const auto invariant = fromBeam.compoundFourMomentum.mag();
    Require(Close(invariant, fromBeam.compoundMass + fromBeam.excitationEnergy, 1.0e-12),
            "Compound four-momentum invariant is inconsistent.");
    const auto fromExcitation = ReactionKinematics::FromExcitationEnergy(
        18, 40, 69, 169, fromBeam.excitationEnergy);
    Require(Close(fromExcitation.beamKineticEnergy, 200.0 * MeV, 1.0e-10),
            "E* to laboratory-energy inversion failed.");

    G4Random::setTheSeed(20260811);
    G4BosonConstructor().ConstructParticle();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4BaryonConstructor().ConstructParticle();
    G4IonConstructor().ConstructParticle();
    G4ShortLivedConstructor().ConstructParticle();
    auto* particleTable = G4ParticleTable::GetParticleTable();
    auto* genericIon = G4GenericIon::GenericIonDefinition();
    genericIon->SetProcessManager(new G4ProcessManager(genericIon));
    particleTable->SetGenericIon(genericIon);
    particleTable->SetReadiness();
    G4ExcitationHandler excitationHandler;
    excitationHandler.Initialise();
    const G4Fragment compound(fromBeam.compoundA, fromBeam.compoundZ,
                              fromBeam.compoundFourMomentum);
    auto* products = excitationHandler.BreakItUp(compound);
    Require(products != nullptr && !products->empty(),
            "Excitation handler returned no products.");
    G4LorentzVector productSum;
    G4int baryonSum = 0;
    G4int nuclearChargeSum = 0;
    for (auto* product : *products) {
      productSum += G4LorentzVector(product->GetMomentum(), product->GetTotalEnergy());
      const auto* definition = product->GetDefinition();
      baryonSum += definition->GetBaryonNumber();
      nuclearChargeSum += definition->GetAtomicNumber();
      delete product;
    }
    delete products;
    Require(baryonSum == fromBeam.compoundA && nuclearChargeSum == fromBeam.compoundZ,
            "De-excitation products violate A/Z conservation.");
    const auto momentumResidual =
        (productSum.vect() - fromBeam.compoundFourMomentum.vect()).mag();
    const auto energyResidual = std::abs(productSum.e() - fromBeam.compoundFourMomentum.e());
    if (momentumResidual >= 1.0e-5 * MeV || energyResidual >= 1.0e-6 * MeV) {
      std::cerr << "De-excitation residuals: |dp|=" << momentumResidual / MeV
                << " MeV/c, |dE|=" << energyResidual / MeV << " MeV\n";
    }
    Require(momentumResidual < 1.0e-5 * MeV && energyResidual < 1.0e-6 * MeV,
            "De-excitation products violate four-momentum conservation.");

    FiniteSolenoidField field(0.20 * m, 1.00 * m, 0.30 * m, 2.50 * tesla,
                              41, 151, 64);
    Require(Close(field.CenterFieldValue(), 2.50 * tesla, 1.0e-4),
            "Center field normalization is outside 1e-4.");
    G4double p1[4] = {0.04 * m, 0.0, 0.80 * m, 0.0};
    G4double p2[4] = {0.0, 0.04 * m, 0.80 * m, 0.0};
    G4double b1[3] = {};
    G4double b2[3] = {};
    field.GetFieldValue(p1, b1);
    field.GetFieldValue(p2, b2);
    Require(Close(b1[0], b2[1], 1.0e-10, 1.0e-12 * tesla) &&
                Close(b1[1], -b2[0], 1.0e-10, 1.0e-12 * tesla) &&
                Close(b1[2], b2[2], 1.0e-10),
            "Finite field is not axisymmetric.");
    const auto before = field.CylindricalField(0.02 * m, 0.30 * m - 0.01 * mm);
    const auto after = field.CylindricalField(0.02 * m, 0.30 * m + 0.01 * mm);
    Require(std::abs(after.second - before.second) < 1.0e-3 * tesla,
            "Entrance fringe field is discontinuous.");
    const auto targetField = field.CylindricalField(0.0, 0.0);
    Require(targetField.second > 0.0 && targetField.second < 2.50 * tesla,
            "Upstream axial fringe field is missing.");
    const auto entranceField = field.CylindricalField(0.10 * m, 0.30 * m);
    const auto exitField = field.CylindricalField(0.10 * m, 1.30 * m);
    Require(entranceField.first < 0.0 && exitField.first > 0.0,
            "Radial fringe field does not focus at entry and defocus at exit.");
    const auto tableEdge = field.CylindricalField(0.0, field.ZMax());
    const auto farDownstream = field.CylindricalField(0.0, field.ZMax() + 20.0 * m);
    Require(std::abs(farDownstream.second) < 1.0e-4 * std::abs(tableEdge.second),
            "Far fringe field does not decay outside the interpolation table.");

    const auto theta = 12.0 * deg;
    Require(Close(0.30 * m * std::tan(theta), 63.766968 * mm, 2.0e-6),
            "Straight drift geometry reference failed.");
  } catch (const std::exception& error) {
    std::cerr << "FAILED: " << error.what() << '\n';
    return 1;
  }
  std::cout << "All solenoid unit checks passed.\n";
  return 0;
}
