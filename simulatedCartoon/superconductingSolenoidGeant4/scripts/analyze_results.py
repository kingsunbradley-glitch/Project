#!/usr/bin/env python3
"""Create ROOT-only diagnostics (PNG/PDF/ROOT); Matplotlib is intentionally unused."""

from __future__ import annotations

import argparse
import pathlib
from collections import defaultdict

import ROOT


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=pathlib.Path, help="Geant4 ROOT output")
    parser.add_argument("--out-dir", type=pathlib.Path)
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kViridis)

    source = ROOT.TFile.Open(str(args.input))
    if not source or source.IsZombie():
        raise OSError(f"Cannot open {args.input}")
    events = source.Get("events")
    tracks = source.Get("tracks")
    summary = source.Get("summary")
    if not events or not tracks or not summary:
        raise RuntimeError("Input is missing events, tracks, or summary ntuple")

    out_dir = args.out_dir or args.input.parent / f"{args.input.stem}_plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    diagnostics = ROOT.TFile(str(out_dir / "diagnostics.root"), "RECREATE")

    theta_max = max(1.0, float(events.GetMaximum("initial_theta_deg")) * 1.08)
    h_theta = ROOT.TH1D("initial_theta", "Physical product angles;#theta_{lab} (deg);weighted counts",
                        100, 0.0, theta_max)
    rho_max = max(0.01, float(events.GetMaximum("entry_rho_m")),
                  float(events.GetMaximum("exit_rho_m"))) * 1.08
    h_entry = ROOT.TH1D("entry_rho", "Entrance/exit radius;#rho (m);events", 100, 0.0, rho_max)
    h_exit = ROOT.TH1D("exit_rho", "Entrance/exit radius;#rho (m);events", 100, 0.0, rho_max)
    h_entry.SetLineColor(ROOT.kAzure + 2)
    h_exit.SetLineColor(ROOT.kOrange + 7)
    h_entry.SetLineWidth(2)
    h_exit.SetLineWidth(2)
    for row in events:
        h_theta.Fill(row.initial_theta_deg, row.weight)
        if row.entered and row.entry_rho_m == row.entry_rho_m:
            h_entry.Fill(row.entry_rho_m)
        if row.status == 1 and row.exit_rho_m == row.exit_rho_m:
            h_exit.Fill(row.exit_rho_m)

    graph = ROOT.TGraphAsymmErrors(summary.GetEntries())
    graph.SetName("acceptance")
    graph.SetTitle("Angular transmission;#theta_{lab} (deg);T(#theta)")
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.75)
    graph.SetLineWidth(2)
    for index, row in enumerate(summary):
        x = row.theta_center_deg
        y = row.efficiency
        graph.SetPoint(index, x, y)
        graph.SetPointError(index, x - row.theta_low_deg, row.theta_high_deg - x,
                            y - row.ci95_low, row.ci95_high - y)

    track_points: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for row in tracks:
        track_points[int(row.event)].append((row.z_m, row.rho_m))
    multigraph = ROOT.TMultiGraph("rho_z_tracks", "Sampled trajectories;z (m);#rho (m)")
    owned_graphs = []
    colors = [ROOT.kAzure + 2, ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kRed + 1, ROOT.kCyan + 2]
    for event_id, points in sorted(track_points.items()):
        trajectory = ROOT.TGraph(len(points))
        trajectory.SetName(f"track_{event_id}")
        trajectory.SetLineColor(colors[event_id % len(colors)])
        trajectory.SetLineWidth(1)
        for point_index, (z_value, rho_value) in enumerate(points):
            trajectory.SetPoint(point_index, z_value, rho_value)
        multigraph.Add(trajectory, "L")
        owned_graphs.append(trajectory)

    canvas = ROOT.TCanvas("transport_diagnostics", "Solenoid transport", 1400, 1000)
    canvas.Divide(2, 2)
    canvas.cd(1)
    h_theta.Draw("HIST")
    canvas.cd(2)
    graph.Draw("AP")
    graph.GetYaxis().SetRangeUser(0.0, 1.05)
    canvas.cd(3)
    if owned_graphs:
        multigraph.Draw("AL")
    canvas.cd(4)
    h_entry.Draw("HIST")
    h_exit.Draw("HIST SAME")
    legend = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
    legend.AddEntry(h_entry, "entrance", "l")
    legend.AddEntry(h_exit, "exit", "l")
    legend.Draw()
    canvas.SaveAs(str(out_dir / "transport_diagnostics.png"))
    canvas.SaveAs(str(out_dir / "transport_diagnostics.pdf"))

    acceptance_canvas = ROOT.TCanvas("acceptance_canvas", "Angular acceptance", 1000, 720)
    graph.Draw("AP")
    graph.GetYaxis().SetRangeUser(0.0, 1.05)
    acceptance_canvas.SetGrid()
    acceptance_canvas.SaveAs(str(out_dir / "angular_acceptance.png"))
    acceptance_canvas.SaveAs(str(out_dir / "angular_acceptance.pdf"))

    diagnostics.cd()
    h_theta.Write()
    h_entry.Write()
    h_exit.Write()
    graph.Write()
    multigraph.Write()
    canvas.Write()
    acceptance_canvas.Write()
    diagnostics.Close()
    source.Close()
    print(f"ROOT diagnostics written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
