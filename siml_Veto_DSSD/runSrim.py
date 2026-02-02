#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WSL + SRIM-2013 TRIM batch runner (stable template)

User request:
- Ions: 1H, 2H, 3H, 3He, 4He
- Normal incidence (Angle = 0)
- Injected energy: 0.0–10.0 MeV, step 0.1 MeV
- 1000 ions per energy point
- Use the SAME TRIM.IN style that you already confirmed can run (2-layer Si 300um + 300um),
  because your new 4-layer (dead/active split) triggered TRIM "Run-time error 13: Type mismatch".

Notes:
- SRIM-2013 may write IONIZ.txt either in:
    C:\srim\IONIZ.txt
    C:\srim\SRIM Outputs\IONIZ.txt
  We monitor BOTH and archive the detected one.
- 0.0 MeV is not physically meaningful for TRIM. We map 0.0 MeV -> MIN_KEV (default 1 keV).
  If your TRIM doesn't like such low energy, set MIN_KEV=10 or 100.

Outputs:
- ./EnergyScan_Results/<ion>/IONIZ_E?.?MeV.txt (plus TRIMOUT / RANGE_3D / BACKSCAT / TRANSMIT / SPUTTER if present)
- ./EnergyScan_Results/<ion>/_logs/log_*.txt
- ./EnergyScan_Results/<ion>/failures.txt (if any)
"""

import os
import shutil
import subprocess
import time
from datetime import datetime
from typing import Optional, Tuple, List


# =========================
# 0) USER CONFIG
# =========================

# SRIM install path on Windows (WSL mount)
SRIM_DIR_WSL = "/mnt/c/srim"
TRIM_EXE_WIN = r"C:\srim\TRIM.exe"

# Energy scan
E_START_MEV = 0.0
E_STOP_MEV  = 10.0
E_STEP_MEV  = 0.1

# TRIM statistics per energy
N_IONS_PER_RUN = 1000

# If e_mev == 0.0, map to this keV to avoid TRIM parsing issues
MIN_KEV_FOR_ZERO = 1  # change to 10 or 100 if needed

# Target (stable template you can run): two Si layers, each 300 µm
# 300 µm = 3,000,000 Angstroms
LAYER_WIDTH_ANGSTROMS = 3_000_000
SI_DENSITY = 2.321

# Ions to simulate
IONS = {
    "1H":  (1, 1.008),
    "2H":  (1, 2.014),
    "3H":  (1, 3.016),
    "3He": (2, 3.016),
    "4He": (2, 4.003),
}

# Resume behavior: if output already exists, skip that energy
SKIP_IF_EXISTS = True


# =========================
# 1) INTERNAL CONFIG
# =========================

TRIM_IN_PATH = os.path.join(SRIM_DIR_WSL, "TRIM.IN")
TRIMAUTO_PATH = os.path.join(SRIM_DIR_WSL, "TRIMAUTO")
OUTPUT_DIR_WSL = os.path.join(SRIM_DIR_WSL, "SRIM Outputs")

IONIZ_CANDIDATES = [
    os.path.join(SRIM_DIR_WSL, "IONIZ.txt"),
    os.path.join(OUTPUT_DIR_WSL, "IONIZ.txt"),
]

ARCHIVE_FILES = [
    ("TRIMOUT.txt", os.path.join(OUTPUT_DIR_WSL, "TRIMOUT.txt")),
    ("RANGE_3D.txt", os.path.join(OUTPUT_DIR_WSL, "RANGE_3D.txt")),
    ("BACKSCAT.txt", os.path.join(OUTPUT_DIR_WSL, "BACKSCAT.txt")),
    ("TRANSMIT.txt", os.path.join(OUTPUT_DIR_WSL, "TRANSMIT.txt")),
    ("SPUTTER.txt", os.path.join(OUTPUT_DIR_WSL, "SPUTTER.txt")),
]

# Robustness knobs
OUTPUT_WAIT_TIMEOUT_S = 300
PROCESS_HARD_TIMEOUT_S = 1200
OUTPUT_STABLE_S = 2.0
POST_EXIT_GRAB_S = 8.0
EXIT_GRACE_S = 6
MAX_RETRIES = 1
SLEEP_BETWEEN_RUNS_S = 0.2

CURRENT_DIR = os.getcwd()
BASE_RESULTS_DIR = os.path.join(CURRENT_DIR, "EnergyScan_Results")


# =========================
# 2) Helpers
# =========================

def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dirs(path: str):
    os.makedirs(path, exist_ok=True)

def ensure_trimauto():
    # Windows CRLF required by SRIM's TRIMAUTO parsing
    with open(TRIMAUTO_PATH, "w", encoding="ascii", newline="\r\n") as f:
        f.write("1\r\n")

def kill_trim_processes():
    subprocess.run(["taskkill.exe", "/IM", "TRIM.exe", "/F", "/T"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["taskkill.exe", "/IM", "SRIM.exe", "/F", "/T"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def cleanup_old_outputs():
    # Remove old IONIZ candidates (overwritten per run)
    for p in IONIZ_CANDIDATES:
        if os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                kill_trim_processes()
                time.sleep(0.5)
                try:
                    os.remove(p)
                except OSError:
                    pass

def _ioniz_has_marker(path: str) -> bool:
    # Marker present in valid IONIZ files
    try:
        with open(path, "rb") as f:
            head = f.read(8192)
        return b"Total Ions calculated" in head
    except OSError:
        return False

def _best_ioniz_path() -> Optional[str]:
    good = []
    for p in IONIZ_CANDIDATES:
        if os.path.exists(p):
            try:
                size = os.path.getsize(p)
                mtime = os.path.getmtime(p)
            except OSError:
                continue
            if size > 0 and _ioniz_has_marker(p):
                good.append((mtime, p))
    if not good:
        return None
    good.sort(reverse=True)
    return good[0][1]

def fmt_energy(e_mev: float) -> str:
    return f"{e_mev:.1f}".replace(".", "p")

def energies_mev() -> List[float]:
    n_steps = int(round((E_STOP_MEV - E_START_MEV) / E_STEP_MEV))
    return [round(E_START_MEV + i * E_STEP_MEV, 10) for i in range(n_steps + 1)]

def mev_to_kev_int(e_mev: float) -> int:
    keV = int(round(e_mev * 1000.0))
    return MIN_KEV_FOR_ZERO if keV <= 0 else keV


# =========================
# 3) TRIM.IN (stable 2-layer template)
# =========================

def generate_trim_in(ion_z: int, ion_m: float, energy_kev: int):
    """
    This is the 2-layer Si(300um)+Si(300um) template that you already ran successfully.
    Normal incidence: Angle=0
    PlotType=5: no plot (reduces GUI issues)
    """
    plot_max = LAYER_WIDTH_ANGSTROMS * 2

    content = f"""==> SRIM-2013 This file controls TRIM Calculations.
Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.
     {ion_z}      {ion_m}      {energy_kev}       0    {N_IONS_PER_RUN}        0    10000
Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders
                      1                                   0       0
Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file
                          1       1           1       1               0                               0
Target material : Number of Elements & Layers
"Si(300um) + Si(300um)                "       1               2
PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]
       5                         0            {plot_max}
Target Elements:    Z   Mass(amu)
Atom 1 = Si =      14      28.09
Layer   Layer Name /               Width Density     Si(14)
Numb.   Description                (Ang) (g/cm3)    Stoich
 1      "Si_Top    "           {LAYER_WIDTH_ANGSTROMS}  {SI_DENSITY}       1
 2      "Si_Bottom "           {LAYER_WIDTH_ANGSTROMS}  {SI_DENSITY}       1
0  Target layer phases (0=Solid, 1=Gas)
0 0
Target Compound Corrections (Bragg)
 1   1
Individual target atom displacement energies (eV)
      15
Individual target atom lattice binding energies (eV)
       2
Individual target atom surface binding energies (eV)
     4.7
Stopping Power Version (1=2011, 0=2011)
 0
"""
    with open(TRIM_IN_PATH, "w", encoding="latin-1", newline="\r\n") as f:
        f.write(content)


# =========================
# 4) Run + detect outputs
# =========================

def wait_for_output(proc: subprocess.Popen, logf) -> Tuple[bool, str, Optional[str]]:
    t0 = time.monotonic()
    last_size = {p: None for p in IONIZ_CANDIDATES}
    stable_since = {p: None for p in IONIZ_CANDIDATES}
    last_report = 0.0

    def update_stable(p: str, now: float) -> bool:
        if not os.path.exists(p):
            last_size[p] = None
            stable_since[p] = None
            return False
        try:
            size = os.path.getsize(p)
        except OSError:
            return False
        if size <= 0:
            last_size[p] = size
            stable_since[p] = None
            return False
        if last_size[p] is None or size != last_size[p]:
            last_size[p] = size
            stable_since[p] = None
            return False
        if stable_since[p] is None:
            stable_since[p] = now
            return False
        if (now - stable_since[p]) >= OUTPUT_STABLE_S and _ioniz_has_marker(p):
            return True
        return False

    while True:
        now = time.monotonic()
        elapsed = now - t0

        if elapsed - last_report >= 5.0:
            last_report = elapsed
            sizes = []
            for p in IONIZ_CANDIDATES:
                s = 0
                if os.path.exists(p):
                    try:
                        s = os.path.getsize(p)
                    except OSError:
                        s = -1
                sizes.append(f"{os.path.basename(p)}={s}")
            msg = f"[{_ts()}]    ...waiting IONIZ ({', '.join(sizes)}, elapsed={int(elapsed)}s)\n"
            try:
                logf.write(msg.encode("utf-8", errors="ignore"))
                logf.flush()
            except Exception:
                pass

        for p in IONIZ_CANDIDATES:
            if update_stable(p, now):
                return True, f"output_ready({p})", p

        rc = proc.poll()
        if rc is not None:
            # after process exit, keep probing for a short while
            end_wait = time.monotonic() + POST_EXIT_GRAB_S
            while time.monotonic() < end_wait:
                pbest = _best_ioniz_path()
                if pbest is not None:
                    return True, f"proc_exit_rc={rc}_accept_output({pbest})", pbest
                time.sleep(0.3)
            return False, f"proc_exit_rc={rc}_no_output", None

        if elapsed > OUTPUT_WAIT_TIMEOUT_S:
            return False, f"output_wait_timeout>{OUTPUT_WAIT_TIMEOUT_S}s", None
        if elapsed > PROCESS_HARD_TIMEOUT_S:
            return False, f"hard_timeout>{PROCESS_HARD_TIMEOUT_S}s", None

        time.sleep(0.2)

def run_trim_once(log_path: str) -> Tuple[bool, str, Optional[str]]:
    ensure_dirs(OUTPUT_DIR_WSL)

    with open(log_path, "wb") as logf:
        proc = subprocess.Popen(
            ["cmd.exe", "/C", TRIM_EXE_WIN],
            cwd=SRIM_DIR_WSL,
            stdout=logf,
            stderr=logf
        )

        ok, reason, ioniz_path = wait_for_output(proc, logf)

        if ok:
            try:
                proc.wait(timeout=EXIT_GRACE_S)
            except subprocess.TimeoutExpired:
                kill_trim_processes()
                try:
                    proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    pass
            return True, reason, ioniz_path

        kill_trim_processes()
        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            pass
        return False, reason, None


# =========================
# 5) Archive outputs
# =========================

def archive_outputs(out_dir: str, ion_name: str, e_mev: float, ioniz_path: str):
    tag = f"E{fmt_energy(e_mev)}MeV"
    shutil.copy2(ioniz_path, os.path.join(out_dir, f"IONIZ_{ion_name}_{tag}.txt"))

    for shortname, src in ARCHIVE_FILES:
        if os.path.exists(src):
            base, ext = os.path.splitext(shortname)
            shutil.copy2(src, os.path.join(out_dir, f"{base}_{ion_name}_{tag}{ext}"))


# =========================
# 6) Main
# =========================

def main():
    ensure_dirs(BASE_RESULTS_DIR)
    ensure_trimauto()

    print(f"[{_ts()}] SRIM_DIR_WSL: {SRIM_DIR_WSL}")
    print(f"[{_ts()}] Results base: {BASE_RESULTS_DIR}")
    print(f"[{_ts()}] Energy scan: {E_START_MEV}–{E_STOP_MEV} MeV step {E_STEP_MEV} MeV")
    print(f"[{_ts()}] N ions per point: {N_IONS_PER_RUN}")
    print(f"[{_ts()}] IONIZ candidates:")
    for p in IONIZ_CANDIDATES:
        print(f"   - {p}")

    for ion_name, (z, m) in IONS.items():
        out_dir = os.path.join(BASE_RESULTS_DIR, ion_name)
        log_dir = os.path.join(out_dir, "_logs")
        ensure_dirs(out_dir)
        ensure_dirs(log_dir)

        print(f"\n========== Ion: {ion_name} (Z={z}, M={m}) ==========")

        failures = []
        for e_mev in energies_mev():
            tag = f"E{fmt_energy(e_mev)}MeV"
            out_ioniz = os.path.join(out_dir, f"IONIZ_{ion_name}_{tag}.txt")
            if SKIP_IF_EXISTS and os.path.exists(out_ioniz) and os.path.getsize(out_ioniz) > 0:
                print(f"[{_ts()}] -> {ion_name} {e_mev:.1f} MeV: skip (exists)")
                continue

            e_keV = mev_to_kev_int(e_mev)
            print(f"[{_ts()}] -> {ion_name} {e_mev:.1f} MeV ({e_keV} keV)")

            # write TRIM.IN using stable 2-layer template
            generate_trim_in(z, m, e_keV)

            success = False
            last_reason = "not_run"

            for attempt in range(MAX_RETRIES + 1):
                cleanup_old_outputs()
                ensure_trimauto()

                log_path = os.path.join(
                    log_dir,
                    f"log_{ion_name}_{tag}_try{attempt}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
                )

                ok, reason, ioniz_path = run_trim_once(log_path)
                last_reason = reason

                if ok and ioniz_path and os.path.exists(ioniz_path) and os.path.getsize(ioniz_path) > 0:
                    if _ioniz_has_marker(ioniz_path):
                        archive_outputs(out_dir, ion_name, e_mev, ioniz_path)
                        print(f"[{_ts()}]    ✅ done (reason={reason})")
                        success = True
                        break

                print(f"[{_ts()}]    ⚠️ fail/hang handled (reason={reason}) log={os.path.basename(log_path)}")
                time.sleep(1.0)

            if not success:
                failures.append((e_mev, last_reason))
                print(f"[{_ts()}]    ❌ skipped: {ion_name} {e_mev:.1f} MeV (last_reason={last_reason})")

            time.sleep(SLEEP_BETWEEN_RUNS_S)

        if failures:
            fail_path = os.path.join(out_dir, "failures.txt")
            with open(fail_path, "w", encoding="utf-8") as f:
                for e_mev, reason in failures:
                    f.write(f"{ion_name}\t{e_mev:.1f} MeV\t{reason}\n")
            print(f"[{_ts()}] Ion {ion_name} finished with failures -> {fail_path}")
        else:
            print(f"[{_ts()}] Ion {ion_name} finished successfully.")

    print(f"\n[{_ts()}] ALL DONE.")


if __name__ == "__main__":
    main()
