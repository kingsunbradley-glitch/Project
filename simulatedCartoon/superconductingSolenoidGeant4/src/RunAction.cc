#include "RunAction.hh"

#include "EventInformation.hh"
#include "SimulationConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>

namespace {
std::pair<G4double, G4double> WilsonInterval(G4int passed, G4int total) {
  if (total <= 0) return {0.0, 0.0};
  constexpr G4double z = 1.959963984540054;
  const auto n = static_cast<G4double>(total);
  const auto p = static_cast<G4double>(passed) / n;
  const auto denominator = 1.0 + z * z / n;
  const auto center = (p + z * z / (2.0 * n)) / denominator;
  const auto half = z * std::sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) /
                    denominator;
  return {std::max(0.0, center - half), std::min(1.0, center + half)};
}
}  // namespace

RunAction::RunAction() {
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(0);
  analysis->SetNtupleMerging(false);

  analysis->CreateNtuple("events", "Evaporation-residue transport events");
  analysis->CreateNtupleIColumn("event");
  analysis->CreateNtupleIColumn("channel_n");
  analysis->CreateNtupleIColumn("residual_Z");
  analysis->CreateNtupleIColumn("residual_A");
  analysis->CreateNtupleDColumn("weight");
  analysis->CreateNtupleIColumn("absolute_xs");
  analysis->CreateNtupleDColumn("beam_energy_MeV");
  analysis->CreateNtupleDColumn("excitation_MeV");
  analysis->CreateNtupleDColumn("initial_theta_deg");
  analysis->CreateNtupleDColumn("initial_phi_deg");
  analysis->CreateNtupleDColumn("initial_energy_MeV");
  analysis->CreateNtupleDColumn("initial_charge_e");
  analysis->CreateNtupleDColumn("initial_vrho_m_s");
  analysis->CreateNtupleDColumn("initial_vphi_m_s");
  analysis->CreateNtupleDColumn("initial_vz_m_s");
  analysis->CreateNtupleIColumn("entered");
  analysis->CreateNtupleDColumn("entry_rho_m");
  analysis->CreateNtupleDColumn("entry_phi_deg");
  analysis->CreateNtupleDColumn("entry_energy_MeV");
  analysis->CreateNtupleDColumn("entry_charge_e");
  analysis->CreateNtupleDColumn("exit_rho_m");
  analysis->CreateNtupleDColumn("exit_phi_deg");
  analysis->CreateNtupleDColumn("exit_theta_deg");
  analysis->CreateNtupleDColumn("exit_energy_MeV");
  analysis->CreateNtupleDColumn("exit_charge_e");
  analysis->CreateNtupleIColumn("status");
  analysis->FinishNtuple();

  analysis->CreateNtuple("tracks", "Sampled trajectories in cylindrical coordinates");
  analysis->CreateNtupleIColumn("event");
  analysis->CreateNtupleIColumn("step");
  analysis->CreateNtupleDColumn("rho_m");
  analysis->CreateNtupleDColumn("phi_deg");
  analysis->CreateNtupleDColumn("z_m");
  analysis->CreateNtupleDColumn("vrho_m_s");
  analysis->CreateNtupleDColumn("vphi_m_s");
  analysis->CreateNtupleDColumn("vz_m_s");
  analysis->CreateNtupleDColumn("energy_MeV");
  analysis->CreateNtupleDColumn("charge_e");
  analysis->CreateNtupleDColumn("time_ns");
  analysis->FinishNtuple();

  analysis->CreateNtuple("summary", "Angular acceptance summary");
  analysis->CreateNtupleIColumn("bin");
  analysis->CreateNtupleDColumn("theta_low_deg");
  analysis->CreateNtupleDColumn("theta_high_deg");
  analysis->CreateNtupleDColumn("theta_center_deg");
  analysis->CreateNtupleIColumn("generated");
  analysis->CreateNtupleIColumn("entered");
  analysis->CreateNtupleIColumn("transmitted");
  analysis->CreateNtupleDColumn("efficiency");
  analysis->CreateNtupleDColumn("ci95_low");
  analysis->CreateNtupleDColumn("ci95_high");
  analysis->CreateNtupleDColumn("weighted_efficiency");
  analysis->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
  auto& cfg = SimulationConfig::Instance();
  cfg.Validate();
  cfg.Print();
  std::fill(std::begin(terminalCounts_), std::end(terminalCounts_), 0);
  std::fill(std::begin(channelGeneratedWeight_), std::end(channelGeneratedWeight_), 0.0);
  std::fill(std::begin(channelTransmittedWeight_), std::end(channelTransmittedWeight_), 0.0);
  maxObservedTransmittedThetaDeg_ = -1.0;
  bins_.assign(cfg.generatorMode == SimulationConfig::GeneratorMode::AngleScan ? cfg.thetaBins
                                                                               : 1,
               {});
  const std::filesystem::path output(cfg.outputFile.c_str());
  if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
  outputBase_ = output;
  outputBase_.replace_extension();
  eventCsv_.open(outputBase_.string() + "_events.csv");
  trackCsv_.open(outputBase_.string() + "_tracks.csv");
  eventCsv_ << std::setprecision(12)
            << "event,channel_n,residual_Z,residual_A,weight,absolute_xs,beam_energy_MeV,"
               "excitation_MeV,initial_theta_deg,initial_phi_deg,initial_energy_MeV,"
               "initial_charge_e,initial_vrho_m_s,initial_vphi_m_s,initial_vz_m_s,entered,"
               "entry_rho_m,entry_phi_deg,entry_energy_MeV,entry_charge_e,exit_rho_m,"
               "exit_phi_deg,exit_theta_deg,exit_energy_MeV,exit_charge_e,status\n";
  trackCsv_ << std::setprecision(12)
            << "event,step,rho_m,phi_deg,z_m,vrho_m_s,vphi_m_s,vz_m_s,energy_MeV,charge_e,time_ns\n";
  G4AnalysisManager::Instance()->OpenFile(cfg.outputFile);
}

void RunAction::RecordEvent(G4int eventId, const EventInformation& info) {
  const auto& cfg = SimulationConfig::Instance();
  G4int bin = 0;
  if (cfg.generatorMode == SimulationConfig::GeneratorMode::AngleScan) {
    const auto relative = (info.initialThetaDeg * deg - cfg.thetaMin) /
                          (cfg.thetaMax - cfg.thetaMin);
    bin = std::clamp(static_cast<G4int>(relative * cfg.thetaBins), 0, cfg.thetaBins - 1);
  }
  auto& counter = bins_.at(bin);
  ++counter.generated;
  counter.generatedWeight += info.eventWeight;
  if (info.entered) ++counter.entered;
  if (info.status == FinalStatus::Transmitted) {
    ++counter.transmitted;
    counter.transmittedWeight += info.eventWeight;
    maxObservedTransmittedThetaDeg_ =
        std::max(maxObservedTransmittedThetaDeg_, info.initialThetaDeg);
  }
  const auto status = static_cast<G4int>(info.status);
  if (status >= 0 && status < 7) ++terminalCounts_[status];
  const auto channel = info.channelNeutrons == 4 ? 0 : 1;
  channelGeneratedWeight_[channel] += info.eventWeight;
  if (info.status == FinalStatus::Transmitted) {
    channelTransmittedWeight_[channel] += info.eventWeight;
  }
  eventCsv_ << eventId << ',' << info.channelNeutrons << ',' << info.residualZ << ','
            << info.residualA << ',' << info.eventWeight << ','
            << (info.absoluteCrossSection ? 1 : 0) << ',' << info.beamEnergyMeV << ','
            << info.excitationMeV << ',' << info.initialThetaDeg << ',' << info.initialPhiDeg
            << ',' << info.initialEnergyMeV << ',' << info.initialCharge << ','
            << info.initialVr << ',' << info.initialVphi << ',' << info.initialVz << ','
            << (info.entered ? 1 : 0) << ',' << info.entryRhoM << ',' << info.entryPhiDeg
            << ',' << info.entryEnergyMeV << ',' << info.entryCharge << ',' << info.exitRhoM
            << ',' << info.exitPhiDeg << ',' << info.exitThetaDeg << ',' << info.exitEnergyMeV
            << ',' << info.exitCharge << ',' << static_cast<G4int>(info.status) << '\n';
}

void RunAction::RecordTrack(G4int eventId, G4int step, G4double rhoM,
                            G4double phiDeg, G4double zM, G4double vrho,
                            G4double vphi, G4double vz, G4double energyMeV,
                            G4double charge, G4double timeNs) {
  trackCsv_ << eventId << ',' << step << ',' << rhoM << ',' << phiDeg << ',' << zM << ','
            << vrho << ',' << vphi << ',' << vz << ',' << energyMeV << ',' << charge << ','
            << timeNs << '\n';
}

void RunAction::EndOfRunAction(const G4Run*) {
  auto* analysis = G4AnalysisManager::Instance();
  const auto& cfg = SimulationConfig::Instance();
  std::ofstream summaryCsv(outputBase_.string() + "_summary.csv");
  summaryCsv << std::setprecision(12)
             << "bin,theta_low_deg,theta_high_deg,theta_center_deg,generated,entered,"
                "transmitted,efficiency,ci95_low,ci95_high,weighted_efficiency\n";
  G4double outerOnePercentLow = -1.0;
  G4double outerOnePercentHigh = -1.0;
  for (std::size_t i = 0; i < bins_.size(); ++i) {
    const auto& counter = bins_[i];
    G4double low = 0.0;
    G4double high = pi;
    if (cfg.generatorMode == SimulationConfig::GeneratorMode::AngleScan) {
      low = cfg.thetaMin + i * (cfg.thetaMax - cfg.thetaMin) / cfg.thetaBins;
      high = cfg.thetaMin + (i + 1) * (cfg.thetaMax - cfg.thetaMin) / cfg.thetaBins;
    }
    const auto efficiency = counter.generated > 0
                                ? static_cast<G4double>(counter.transmitted) / counter.generated
                                : 0.0;
    const auto weighted = counter.generatedWeight > 0.0
                              ? counter.transmittedWeight / counter.generatedWeight
                              : 0.0;
    const auto [ciLow, ciHigh] = WilsonInterval(counter.transmitted, counter.generated);
    analysis->FillNtupleIColumn(2, 0, static_cast<G4int>(i));
    analysis->FillNtupleDColumn(2, 1, low / deg);
    analysis->FillNtupleDColumn(2, 2, high / deg);
    analysis->FillNtupleDColumn(2, 3, 0.5 * (low + high) / deg);
    analysis->FillNtupleIColumn(2, 4, counter.generated);
    analysis->FillNtupleIColumn(2, 5, counter.entered);
    analysis->FillNtupleIColumn(2, 6, counter.transmitted);
    analysis->FillNtupleDColumn(2, 7, efficiency);
    analysis->FillNtupleDColumn(2, 8, ciLow);
    analysis->FillNtupleDColumn(2, 9, ciHigh);
    analysis->FillNtupleDColumn(2, 10, weighted);
    analysis->AddNtupleRow(2);
    summaryCsv << i << ',' << low / deg << ',' << high / deg << ','
               << 0.5 * (low + high) / deg << ',' << counter.generated << ','
               << counter.entered << ',' << counter.transmitted << ',' << efficiency << ','
               << ciLow << ',' << ciHigh << ',' << weighted << '\n';
    if (efficiency >= 0.01) {
      if (outerOnePercentLow < 0.0) outerOnePercentLow = low / deg;
      outerOnePercentHigh = high / deg;
    }
  }

  const auto geometryAngle = std::atan2(cfg.solenoidRadius, cfg.targetDistance) / deg;
  G4cout << "\nRun summary: theta_geom=" << geometryAngle << " deg\n";
  G4double sigmaOutTotal = 0.0;
  for (G4int channel = 0; channel < 2; ++channel) {
    const auto efficiency = channelGeneratedWeight_[channel] > 0.0
                                ? channelTransmittedWeight_[channel] /
                                      channelGeneratedWeight_[channel]
                                : 0.0;
    G4cout << (channel == 0 ? "4n" : "5n") << " efficiency=" << efficiency;
    if (cfg.sigma4nMicrobarn > 0.0 && cfg.sigma5nMicrobarn > 0.0) {
      const auto sigma = channel == 0 ? cfg.sigma4nMicrobarn : cfg.sigma5nMicrobarn;
      sigmaOutTotal += sigma * efficiency;
      G4cout << ", sigma_out=" << sigma * efficiency << " microbarn";
    }
    G4cout << '\n';
  }
  if (cfg.sigma4nMicrobarn > 0.0 && cfg.sigma5nMicrobarn > 0.0) {
    G4cout << "sigma_out_total=" << sigmaOutTotal << " microbarn\n";
  }
  G4cout << "terminal counts [transmitted, aperture, not-entered, stopped, backward, world] = ["
         << terminalCounts_[1] << ", " << terminalCounts_[2] << ", "
         << terminalCounts_[3] << ", " << terminalCounts_[4] << ", "
         << terminalCounts_[5] << ", " << terminalCounts_[6] << "]\n" << G4endl;

  std::ofstream summaryJson(outputBase_.string() + "_summary.json");
  summaryJson << std::setprecision(12)
              << "{\n  \"theta_geom_deg\": " << geometryAngle
              << ",\n  \"max_observed_pass_angle_deg\": "
              << maxObservedTransmittedThetaDeg_ << ",\n  \"outer_T_ge_1pct_deg\": ";
  if (cfg.generatorMode == SimulationConfig::GeneratorMode::AngleScan) {
    summaryJson << '[' << outerOnePercentLow << ", " << outerOnePercentHigh << ']';
  } else {
    summaryJson << "null";
  }
  summaryJson << ',';
  if (cfg.sigma4nMicrobarn > 0.0 && cfg.sigma5nMicrobarn > 0.0) {
    summaryJson << "\n  \"sigma_out_total_microbarn\": " << sigmaOutTotal << ',';
  }
  summaryJson << "\n  \"channels\": {\n";
  for (G4int channel = 0; channel < 2; ++channel) {
    const auto efficiency = channelGeneratedWeight_[channel] > 0.0
                                ? channelTransmittedWeight_[channel] /
                                      channelGeneratedWeight_[channel]
                                : 0.0;
    const auto sigma = channel == 0 ? cfg.sigma4nMicrobarn : cfg.sigma5nMicrobarn;
    summaryJson << "    \"" << (channel == 0 ? "4n" : "5n")
                << "\": {\"efficiency\": " << efficiency;
    if (cfg.sigma4nMicrobarn > 0.0 && cfg.sigma5nMicrobarn > 0.0) {
      summaryJson << ", \"sigma_out_microbarn\": " << sigma * efficiency;
    }
    summaryJson << "}" << (channel == 0 ? "," : "") << '\n';
  }
  summaryJson << "  },\n  \"terminal_counts\": {\"transmitted\": " << terminalCounts_[1]
              << ", \"aperture\": " << terminalCounts_[2]
              << ", \"not_entered\": " << terminalCounts_[3]
              << ", \"stopped\": " << terminalCounts_[4]
              << ", \"backward\": " << terminalCounts_[5]
              << ", \"world_exit\": " << terminalCounts_[6] << "}\n}\n";
  eventCsv_.close();
  trackCsv_.close();
  summaryCsv.close();
  summaryJson.close();

  analysis->Write();
  analysis->CloseFile();
}
