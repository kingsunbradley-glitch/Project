#!/usr/bin/env python3
"""Small deterministic transport checks run by CTest; no third-party Python modules."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import subprocess
import tempfile


BASE_MACRO = """
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/sim/geometry/solenoidRadius {radius} m
/sim/geometry/solenoidLength {length} m
/sim/geometry/targetDistance {distance} m
/sim/geometry/maxStep 1 mm
/sim/field/centerField {field} tesla
/sim/gas/enabled {gas}
/sim/gas/pressure {pressure} pascal
/sim/gas/temperature 300 kelvin
/sim/reaction/projectileZA 18 40
/sim/reaction/targetZA 69 169
/sim/reaction/beamEnergyLab 200 MeV
/sim/reaction/mode angleScan
/sim/scan/thetaMin 10 deg
/sim/scan/thetaMax 10.001 deg
/sim/scan/bins 1
/sim/scan/eventsPerBin {events}
/sim/output/fileName {output}
/sim/output/trackSampleCount {track_sample}
/sim/output/postExitTrackLength {post_exit} m
/sim/output/randomSeed 314159
/run/initialize
/run/beamOn {events}
"""


def run_case(executable: pathlib.Path, directory: pathlib.Path, name: str, *, gas: bool,
             pressure: float = 100.0, events: int = 6, track_sample: int = 0,
             post_exit: float = 0.0) -> list[dict[str, str]]:
    output = directory / f"{name}.root"
    macro = directory / f"{name}.mac"
    macro.write_text(BASE_MACRO.format(radius=0.5, length=0.3, distance=0.2,
                                       field=0.0, gas=str(gas).lower(),
                                       pressure=pressure, events=events, output=output,
                                       track_sample=track_sample, post_exit=post_exit),
                     encoding="utf-8")
    subprocess.run([str(executable), str(macro)], check=True, cwd=directory,
                   stdout=subprocess.DEVNULL)
    with (directory / f"{name}_events.csv").open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, type=pathlib.Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="solenoid-ctest-") as temporary:
        directory = pathlib.Path(temporary)
        first = run_case(args.exe, directory, "repeat_a", gas=False)
        second = run_case(args.exe, directory, "repeat_b", gas=False)
        fields = [key for key in first[0] if key != "event"]
        if [[row[key] for key in fields] for row in first] != [
            [row[key] for key in fields] for row in second
        ]:
            raise AssertionError("Fixed-seed serial runs are not reproducible")

        expected_rho = 0.2 * math.tan(math.radians(10.0005))
        transmitted = [row for row in first if int(row["status"]) == 1]
        if not transmitted:
            raise AssertionError("No event transmitted in zero-field geometry test")
        for row in first:
            if int(row["residual_Z"]) != 87 or int(row["residual_A"]) != 209 - int(row["channel_n"]):
                raise AssertionError("Conditioned 4n/5n residue has the wrong A/Z")
        for row in transmitted:
            if not math.isclose(float(row["entry_rho_m"]), expected_rho,
                                rel_tol=0.0, abs_tol=3.0e-5):
                raise AssertionError("rho_entry != d*tan(theta) in zero field")
            if not math.isclose(float(row["entry_energy_MeV"]),
                                float(row["exit_energy_MeV"]), abs_tol=2.0e-5):
                raise AssertionError("Vacuum zero-field run changed kinetic energy")

        gas_rows = run_case(args.exe, directory, "helium", gas=True, pressure=10000.0)
        gas_transmitted = [row for row in gas_rows if int(row["status"]) == 1]
        if gas_transmitted and not any(float(row["exit_energy_MeV"]) <
                                       float(row["entry_energy_MeV"]) for row in gas_transmitted):
            raise AssertionError("Helium did not reduce mean residue energy")
        for row in gas_rows:
            q = float(row["entry_charge_e"])
            if math.isfinite(q) and not 0.0 <= q <= float(row["residual_Z"]):
                raise AssertionError("Dynamic charge is outside 0 <= q <= Z")

        post_rows = run_case(args.exe, directory, "post_exit", gas=False, events=1,
                             track_sample=1, post_exit=0.1)
        if int(post_rows[0]["status"]) != 1:
            raise AssertionError("Post-exit sampling changed the transmitted classification")
        with (directory / "post_exit_tracks.csv").open(newline="", encoding="utf-8") as stream:
            track_rows = list(csv.DictReader(stream))
        if not track_rows or max(float(row["z_m"]) for row in track_rows) < 0.59:
            raise AssertionError("Sampled trajectory did not continue beyond the exit plane")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
