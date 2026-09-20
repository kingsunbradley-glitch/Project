#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate alpha-decay hindrance factors and build a nuclide-chart dashboard.

The experimental quantity in the HF numerator is read *directly* from the
database column ``T_alpha_half_s``:

    HF = T_alpha_half_s / Tcalc_s

``T_half_s`` (the total half-life) and ``b_alpha_percent`` are retained only
for provenance and consistency checks; they are never used to reconstruct the
partial alpha half-life.

The seven active models are the ones registered by ``method/all/run_all.py``.
Their pure calculation functions are imported directly so that all methods
share the same HF definition and data-quality policy.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DATABASE = BASE_DIR / "ackermann_ground_states.csv"
DEFAULT_METHOD_DIR = BASE_DIR / "method" / "all"
DEFAULT_OUTPUT_DIR = BASE_DIR / "results"

HF_DEFINITION = "HF = T_alpha_half_s / Tcalc_s"
Q_SOURCE = "Qalpha_MeV"
TALPHA_SOURCE = "T_alpha_half_s"

# Fixed, shared, discrete HF palette. Boundaries are intentionally explicit:
# [0.1, 4), [4, 10), [10, 100], with separate under/overflow colors.
HF_COLOR_BINS: tuple[dict[str, Any], ...] = (
    {"min": None, "max": 0.1, "include_max": False, "label": "HF < 0.1", "color": "#313695"},
    {"min": 0.1, "max": 4.0, "include_max": False, "label": "0.1 ≤ HF < 4", "color": "#74add1"},
    {"min": 4.0, "max": 10.0, "include_max": False, "label": "4 ≤ HF < 10", "color": "#fee090"},
    {"min": 10.0, "max": 100.0, "include_max": True, "label": "10 ≤ HF ≤ 100", "color": "#f46d43"},
    {"min": 100.0, "max": None, "include_max": False, "label": "HF > 100", "color": "#a50026"},
)


@dataclass(frozen=True)
class MethodSpec:
    """A normalized interface to one of the local half-life models."""

    name: str
    slug: str
    title: str
    note: str
    uses_l: bool
    calculate: Callable[[int, int, int, float, float], tuple[float, dict[str, Any]]]


def load_module(alias: str, path: Path) -> Any:
    """Import a Python file whose filename is not necessarily import-safe."""
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import method module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


def load_methods(method_dir: Path) -> list[MethodSpec]:
    """Load the seven methods used by method/all/run_all.py."""
    poenaru = load_module("hfdistr_poenaru", method_dir / "PoenaruMethod.py")
    dzr = load_module(
        "hfdistr_dzr", method_dir / "calc_dzr_poenaru_format.py"
    )
    ismail = load_module(
        "hfdistr_ismail", method_dir / "calc_ismail2022_formulaE.py"
    )
    qi = load_module("hfdistr_qi", method_dir / "calc_qi2009_udl.py")
    viola = load_module(
        "hfdistr_viola", method_dir / "calc_viola_seaborg.py"
    )
    xu = load_module("hfdistr_xu", method_dir / "calc_xu2022_unified.py")
    rurarz = load_module(
        "hfdistr_rurarz", method_dir / "calc_poenaru_E. Rurarz_format.py"
    )

    def calc_poenaru(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        result = poenaru.calc_poenaru_log10_t(A, Z, qalpha)
        return float(result["log10_Tcalc_s"]), {
            "parameter_set": "original1980",
            "x": result["x"],
            "N_interval": result["N_interval"],
            "Z_interval": result["Z_interval"],
        }

    def calc_dzr(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        parity = dzr.case_from_ZN(Z, N)
        return float(dzr.logT_dzr(A, Z, qalpha, ell, parity)), {
            "parity_case": parity,
            "l_used": ell,
        }

    def calc_ismail(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        parity = ismail.parity_case(Z, N)
        log10_t = ismail.ismail_log10T(
            A, Z, qalpha, ell=ell, Ac=4, Zc=2, param_set=None
        )
        return float(log10_t), {
            "parameter_set": parity,
            "parity_case": parity,
            "l_used": ell,
            "Ac_used": 4,
            "Zc_used": 2,
        }

    def calc_qi(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        log10_t, chi_prime, rho_prime = qi.qi_udl_log10T(
            A, Z, qalpha, Ac=4, Zc=2, coeff_set="alpha"
        )
        return float(log10_t), {
            "coefficient_set": "alpha",
            "chi_prime": chi_prime,
            "rho_prime": rho_prime,
            "Ac_used": 4,
            "Zc_used": 2,
        }

    def calc_viola(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        return float(viola.viola_seaborg_log10T(A, Z, qalpha)), {}

    def calc_xu(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        log10_t, x_value, bracket, c_value, c_status, h_value, parity = (
            xu.xu_log10T(A, Z, qalpha, ell)
        )
        return float(log10_t), {
            "X": x_value,
            "bracket": bracket,
            "C_ZN": c_value,
            "C_status": c_status,
            "h_blocking": h_value,
            "parity_case": parity,
            "l_used": ell,
        }

    def calc_rurarz(
        A: int, Z: int, N: int, qalpha: float, ell: float
    ) -> tuple[float, dict[str, Any]]:
        result = rurarz.calc_poenaru(
            A, Z, qalpha, rurarz.PARAM_SETS["target"]
        )
        return float(result["log10_T12a_Poenaru_s"]), {
            "parameter_set": "target",
            "x": result["x_calc"],
            "N_interval": result["N_interval"],
            "Z_interval": result["Z_interval"],
        }

    return [
        MethodSpec(
            "Poenaru",
            "poenaru",
            "Poenaru (original1980)",
            "Poenaru–Ivașcu–Mazilu parameters with B6 = −0.16820.",
            False,
            calc_poenaru,
        ),
        MethodSpec(
            "DZR",
            "dzr",
            "DZR",
            "Deng–Zhang–Royer model; l is not tabulated and is assumed.",
            True,
            calc_dzr,
        ),
        MethodSpec(
            "Ismail2022_FormulaE",
            "ismail2022_formula_e",
            "Ismail 2022 — Formula E",
            "Parity-specific Formula E coefficients; l is assumed.",
            True,
            calc_ismail,
        ),
        MethodSpec(
            "Qi2009_UDL",
            "qi2009_udl",
            "Qi 2009 UDL",
            "Universal Decay Law with the alpha-decay coefficient set.",
            False,
            calc_qi,
        ),
        MethodSpec(
            "Viola-Seaborg",
            "viola_seaborg",
            "Viola–Seaborg",
            "Viola–Seaborg implementation from method/all.",
            False,
            calc_viola,
        ),
        MethodSpec(
            "Xu2022_Unified",
            "xu2022_unified",
            "Xu 2022 Unified",
            (
                "C(Z,N) is defined by the local implementation only for Z≤90. "
                "For this database (Z=99–118), C=0 is an outside-domain extrapolation."
            ),
            True,
            calc_xu,
        ),
        MethodSpec(
            "E. Rurarz",
            "e_rurarz",
            "E. Rurarz / Poenaru target",
            "Poenaru-type target parameter set with B5 = B6 = −0.003033.",
            False,
            calc_rurarz,
        ),
    ]


def strip_table_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Strip fixed-width padding without changing numeric values."""
    df = df.copy()
    df.columns = [str(col).strip().replace("\ufeff", "") for col in df.columns]
    for col in df.columns:
        df[col] = df[col].map(
            lambda value: value.strip() if isinstance(value, str) else value
        )
        df[col] = df[col].replace("", pd.NA)
    return df


def parse_boolean(value: Any) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    return None


def prepare_database(
    database: Path, *, experimental_q_only: bool, ell: float
) -> pd.DataFrame:
    """Read, validate, and annotate the Ackermann CSV database."""
    # skipinitialspace is required because quoted fields follow padded commas.
    df = pd.read_csv(database, skipinitialspace=True)
    df = strip_table_strings(df)

    required = {
        "isotope",
        "symbol",
        "A",
        "Z",
        "N",
        "Qalpha_MeV",
        "T_alpha_half_s",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Database is missing required columns: {missing}")

    numeric_columns = [
        "A",
        "Z",
        "N",
        "Qalpha_keV",
        "Qalpha_MeV",
        "Qalpha_err_MeV",
        "Ealpha_MeV",
        "T_half_s",
        "T_half_err_plus_s",
        "T_half_err_minus_s",
        "T_half_err_s",
        "b_alpha_percent",
        "b_alpha_err_plus_percent",
        "b_alpha_err_minus_percent",
        "b_alpha_err_percent",
        "b_SF_percent",
        "b_SF_err_plus_percent",
        "b_SF_err_minus_percent",
        "b_SF_err_percent",
        "T_alpha_half_s",
        "T_SF_half_s",
        "T_alpha_half_err_plus_s",
        "T_alpha_half_err_minus_s",
        "T_alpha_half_err_s",
        "T_SF_half_err_plus_s",
        "T_SF_half_err_minus_s",
        "T_SF_half_err_s",
    ]
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ["A", "Z", "N"]:
        if df[col].isna().any():
            bad = (df.index[df[col].isna()] + 2).tolist()
            raise ValueError(f"Missing {col} at CSV line(s): {bad}")
        rounded = df[col].round().astype(int)
        if not np.allclose(df[col], rounded):
            raise ValueError(f"Non-integral values found in {col}")
        df[col] = rounded

    inconsistent = df["A"] != df["Z"] + df["N"]
    if inconsistent.any():
        lines = (df.index[inconsistent] + 2).tolist()
        raise ValueError(f"A != Z + N at CSV line(s): {lines}")

    if "Qalpha_keV" in df.columns:
        both = df["Qalpha_keV"].notna() & df["Qalpha_MeV"].notna()
        mismatch = both & ~np.isclose(
            df["Qalpha_keV"] / 1000.0,
            df["Qalpha_MeV"],
            rtol=0.0,
            atol=1.0e-12,
        )
        if mismatch.any():
            lines = (df.index[mismatch] + 2).tolist()
            raise ValueError(f"Qalpha_keV/1000 != Qalpha_MeV at lines: {lines}")

    for col in [
        "Qalpha_estimated",
        "Qalpha_is_experimental",
        "is_parenthesized",
        "mass_number_uncertain",
        "T_half_estimated",
        "b_alpha_estimated",
        "T_alpha_half_estimated",
        "b_SF_estimated",
        "T_SF_half_estimated",
        "has_estimated_marker",
    ]:
        if col in df.columns:
            df[col] = df[col].map(parse_boolean).astype("boolean")

    df.insert(0, "source_csv_line", df.index + 2)
    df.insert(
        0,
        "record_id",
        [f"source-row-{line:04d}" for line in df["source_csv_line"]],
    )

    q_estimated = pd.Series(False, index=df.index, dtype=bool)
    if "Qalpha_estimated" in df.columns:
        q_estimated |= df["Qalpha_estimated"].fillna(False).astype(bool)
    if "Qalpha_is_experimental" in df.columns:
        q_estimated |= ~df["Qalpha_is_experimental"].fillna(False).astype(bool)
    df["qalpha_estimated_used"] = q_estimated

    qualifier = df.get(
        "T_alpha_half_qualifier", pd.Series(pd.NA, index=df.index)
    )
    df["talpha_qualified"] = qualifier.notna()
    df["l_assumed"] = float(ell)

    count_by_cell = df.groupby(["Z", "N"])["record_id"].transform("size")
    df["cell_record_count"] = count_by_cell.astype(int)
    df["cell_record_index"] = (
        df.groupby(["Z", "N"], sort=False).cumcount() + 1
    ).astype(int)
    df["duplicate_nz_cell"] = df["cell_record_count"] > 1

    def reasons_for_row(row: pd.Series) -> str:
        reasons: list[str] = []
        qalpha = row["Qalpha_MeV"]
        talpha = row["T_alpha_half_s"]
        if pd.isna(qalpha):
            reasons.append("missing_Qalpha_MeV")
        elif not math.isfinite(float(qalpha)) or float(qalpha) <= 0.0:
            reasons.append("invalid_Qalpha_MeV")
        if pd.isna(talpha):
            reasons.append("missing_T_alpha_half_s")
        elif not math.isfinite(float(talpha)) or float(talpha) <= 0.0:
            reasons.append("invalid_T_alpha_half_s")
        if experimental_q_only and bool(row["qalpha_estimated_used"]):
            reasons.append("nonexperimental_Qalpha_excluded")
        return ";".join(reasons) if reasons else "ok"

    df["calculation_status"] = df.apply(reasons_for_row, axis=1)
    df["eligible_for_hf"] = df["calculation_status"] == "ok"

    # This is a diagnostic only. The reconstructed value is never used for HF.
    df["partial_half_life_reconstructed_s"] = np.nan
    branch_ok = (
        df.get("T_half_s", pd.Series(np.nan, index=df.index)).notna()
        & df.get("b_alpha_percent", pd.Series(np.nan, index=df.index)).notna()
        & (df["b_alpha_percent"] > 0.0)
    )
    df.loc[branch_ok, "partial_half_life_reconstructed_s"] = (
        df.loc[branch_ok, "T_half_s"]
        / (df.loc[branch_ok, "b_alpha_percent"] / 100.0)
    )
    comparable = branch_ok & df["T_alpha_half_s"].notna()
    df["partial_half_life_consistent"] = pd.Series(
        pd.NA, index=df.index, dtype="boolean"
    )
    df.loc[comparable, "partial_half_life_consistent"] = np.isclose(
        df.loc[comparable, "partial_half_life_reconstructed_s"],
        df.loc[comparable, "T_alpha_half_s"],
        rtol=1.0e-6,
        atol=0.0,
    )

    return df


def audit_database_consistency(database: pd.DataFrame) -> pd.DataFrame:
    """Return internal data-consistency issues that can affect interpretation.

    The audit distinguishes hard inconsistencies from source-level warnings. It
    never changes values automatically.
    """
    issues: list[dict[str, Any]] = []

    def add_issue(
        index: int,
        *,
        severity: str,
        check: str,
        observed: Any,
        expected: Any,
        details: str,
    ) -> None:
        row = database.loc[index]
        relative_difference = np.nan
        try:
            obs_float = float(observed)
            exp_float = float(expected)
            if math.isfinite(obs_float) and math.isfinite(exp_float) and exp_float != 0.0:
                relative_difference = (obs_float - exp_float) / exp_float
        except (TypeError, ValueError):
            pass
        issues.append(
            {
                "severity": severity,
                "check": check,
                "record_id": row["record_id"],
                "source_csv_line": int(row["source_csv_line"]),
                "isotope": row["isotope"],
                "observed": observed,
                "expected": expected,
                "relative_difference": relative_difference,
                "details": details,
            }
        )

    raw_missing = database["T_half_raw"].notna() & database["T_half_s"].isna()
    for index in database.index[raw_missing]:
        add_issue(
            index,
            severity="error",
            check="T_half_raw_not_parsed",
            observed=database.at[index, "T_half_raw"],
            expected="finite T_half_s",
            details="Non-empty raw total half-life was not converted to seconds.",
        )

    def exact_qualifier(column: str) -> pd.Series:
        if column not in database.columns:
            return pd.Series(True, index=database.index)
        values = database[column]
        return values.isna() | values.astype(str).str.strip().eq("")

    for branch_name, partial_name, qualifier_name, label in [
        ("b_alpha_percent", "T_alpha_half_s", "b_alpha_qualifier", "alpha"),
        ("b_SF_percent", "T_SF_half_s", "b_SF_qualifier", "SF"),
    ]:
        branch = pd.to_numeric(database[branch_name], errors="coerce")
        total = pd.to_numeric(database["T_half_s"], errors="coerce")
        partial = pd.to_numeric(database[partial_name], errors="coerce")
        positive_branch = branch.notna() & (branch > 0.0)
        missing_partial = positive_branch & total.notna() & partial.isna()
        for index in database.index[missing_partial]:
            add_issue(
                index,
                severity="error",
                check=f"missing_{label}_partial_half_life",
                observed=np.nan,
                expected=float(total[index] / (branch[index] / 100.0)),
                details=f"Positive {label} branch and total half-life are present.",
            )

        comparable = (
            positive_branch
            & total.notna()
            & partial.notna()
            & exact_qualifier(qualifier_name)
        )
        expected_partial = total / (branch / 100.0)
        mismatch = comparable & ~np.isclose(
            partial, expected_partial, rtol=1.0e-6, atol=1.0e-15
        )
        for index in database.index[mismatch]:
            add_issue(
                index,
                severity="error",
                check=f"{label}_partial_half_life_identity",
                observed=float(partial[index]),
                expected=float(expected_partial[index]),
                details=f"Expected T_{label}=T_total/(BR_{label}/100).",
            )

        for side in ["plus", "minus"]:
            total_error = pd.to_numeric(
                database[f"T_half_err_{side}_s"], errors="coerce"
            )
            partial_error = pd.to_numeric(
                database[f"{partial_name[:-2]}_err_{side}_s"], errors="coerce"
            )
            # T_partial = T_total / BR. A positive partial-life error receives
            # the negative BR error, and vice versa, because BR is inverted.
            opposite_side = "minus" if side == "plus" else "plus"
            branch_error_column = (
                f"{branch_name.removesuffix('percent')}err_"
                f"{opposite_side}_percent"
            )
            branch_error = pd.to_numeric(
                database.get(
                    branch_error_column,
                    pd.Series(np.nan, index=database.index),
                ),
                errors="coerce",
            )
            expected_error = expected_partial * np.sqrt(
                (total_error / total) ** 2
                + (branch_error.fillna(0.0) / branch) ** 2
            )
            error_comparable = (
                comparable & total_error.notna() & partial_error.notna()
            )
            error_mismatch = error_comparable & ~np.isclose(
                partial_error, expected_error, rtol=1.0e-6, atol=1.0e-15
            )
            for index in database.index[error_mismatch]:
                add_issue(
                    index,
                    severity="error",
                    check=f"{label}_partial_error_{side}_identity",
                    observed=float(partial_error[index]),
                    expected=float(expected_error[index]),
                    details=(
                        "Expected propagated uncertainty for "
                        "T_partial=T_total/(BR/100), including BR uncertainty "
                        "when available."
                    ),
                )

    for prefix in ["T_half", "T_alpha_half", "T_SF_half"]:
        plus = pd.to_numeric(database[f"{prefix}_err_plus_s"], errors="coerce")
        minus = pd.to_numeric(database[f"{prefix}_err_minus_s"], errors="coerce")
        symmetric = pd.to_numeric(database[f"{prefix}_err_s"], errors="coerce")
        expected = (plus + minus) / 2.0
        comparable = plus.notna() & minus.notna() & symmetric.notna()
        mismatch = comparable & ~np.isclose(
            symmetric, expected, rtol=1.0e-6, atol=1.0e-15
        )
        for index in database.index[mismatch]:
            add_issue(
                index,
                severity="error",
                check=f"{prefix}_symmetric_error_mean",
                observed=float(symmetric[index]),
                expected=float(expected[index]),
                details="Symmetric error should equal the mean of plus/minus errors.",
            )

    q_mev = pd.to_numeric(database["Qalpha_MeV"], errors="coerce")
    q_kev = pd.to_numeric(database["Qalpha_keV"], errors="coerce")
    comparable_q = q_mev.notna() & q_kev.notna()
    mismatch_q = comparable_q & ~np.isclose(
        q_mev, q_kev / 1000.0, rtol=0.0, atol=1.0e-12
    )
    for index in database.index[mismatch_q]:
        add_issue(
            index,
            severity="error",
            check="Qalpha_unit_conversion",
            observed=float(q_mev[index]),
            expected=float(q_kev[index] / 1000.0),
            details="Expected Qalpha_MeV=Qalpha_keV/1000.",
        )

    ealpha = pd.to_numeric(database.get("Ealpha_MeV"), errors="coerce")
    comparable_energy = ealpha.notna() & q_mev.notna()
    expected_q = ealpha * database["A"] / (database["A"] - 4.0)
    mismatch_energy = comparable_energy & ~np.isclose(
        q_mev, expected_q, rtol=1.0e-10, atol=1.0e-12
    )
    for index in database.index[mismatch_energy]:
        add_issue(
            index,
            severity="error",
            check="Ealpha_to_Qalpha_conversion",
            observed=float(q_mev[index]),
            expected=float(expected_q[index]),
            details="Expected Qalpha=Ealpha*A/(A-4).",
        )

    if {"Qalpha_estimated", "Qalpha_is_experimental"} <= set(database.columns):
        estimated = database["Qalpha_estimated"]
        experimental = database["Qalpha_is_experimental"]
        comparable_flags = estimated.notna() & experimental.notna()
        mismatch_flags = comparable_flags & (estimated.astype(bool) == experimental.astype(bool))
        for index in database.index[mismatch_flags]:
            add_issue(
                index,
                severity="error",
                check="Qalpha_quality_flags",
                observed=f"estimated={estimated[index]}, experimental={experimental[index]}",
                expected="opposite boolean values",
                details="Estimated and experimental Qalpha flags should be complementary.",
            )

    modes = (
        database["decay_modes_raw"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.replace(" ", "", regex=False)
    )
    alpha_branch = pd.to_numeric(database["b_alpha_percent"], errors="coerce")
    sf_branch = pd.to_numeric(database["b_SF_percent"], errors="coerce")
    exact_two_mode = (
        modes.isin(["alpha,sf", "sf,alpha"])
        & alpha_branch.notna()
        & sf_branch.notna()
        & exact_qualifier("b_alpha_qualifier")
        & exact_qualifier("b_SF_qualifier")
    )
    branch_sum = alpha_branch + sf_branch
    sum_warning = exact_two_mode & ((branch_sum - 100.0).abs() > 0.5)
    for index in database.index[sum_warning]:
        add_issue(
            index,
            severity="warning",
            check="two_mode_branch_sum",
            observed=float(branch_sum[index]),
            expected=100.0,
            details=(
                "Exact alpha+SF percentages do not sum to 100%; likely source "
                "rounding or unresolved uncertainty, not an automatic correction."
            ),
        )

    columns = [
        "severity",
        "check",
        "record_id",
        "source_csv_line",
        "isotope",
        "observed",
        "expected",
        "relative_difference",
        "details",
    ]
    return pd.DataFrame(issues, columns=columns)


def evaluate_methods(
    database: pd.DataFrame, methods: list[MethodSpec], ell: float
) -> pd.DataFrame:
    """Evaluate every method using the same Qalpha and partial-alpha half-life."""
    rows: list[dict[str, Any]] = []
    eligible = database.loc[database["eligible_for_hf"]].copy()

    for method in methods:
        for _, source in eligible.iterrows():
            A = int(source["A"])
            Z = int(source["Z"])
            N = int(source["N"])
            qalpha = float(source["Qalpha_MeV"])
            talpha = float(source["T_alpha_half_s"])
            if N != A - Z:
                raise ValueError(f"Inconsistent A/Z/N for {source['record_id']}")

            try:
                log10_tcalc, diagnostics = method.calculate(
                    A, Z, N, qalpha, ell
                )
                if not math.isfinite(log10_tcalc):
                    raise ValueError("non-finite log10(Tcalc)")
                tcalc = 10.0**log10_tcalc
                log10_hf = math.log10(talpha) - log10_tcalc
                hf = 10.0**log10_hf
                if not all(math.isfinite(value) and value > 0.0 for value in [tcalc, hf]):
                    raise ValueError("non-finite or non-positive Tcalc/HF")
                method_status = "ok"
            except (ValueError, OverflowError, ZeroDivisionError) as exc:
                raise RuntimeError(
                    f"{method.name} failed for {source['record_id']} "
                    f"({source['isotope']}): {exc}"
                ) from exc

            rows.append(
                {
                    "record_id": source["record_id"],
                    "source_csv_line": int(source["source_csv_line"]),
                    "isotope": source["isotope"],
                    "state_label": source.get("state_label", pd.NA),
                    "symbol": source["symbol"],
                    "A": A,
                    "Z": Z,
                    "N": N,
                    "method": method.name,
                    "method_slug": method.slug,
                    "Qalpha_MeV": qalpha,
                    "Ealpha_MeV": source.get("Ealpha_MeV", np.nan),
                    "Qalpha_err_MeV": source.get("Qalpha_err_MeV", np.nan),
                    "qalpha_estimated": bool(source["qalpha_estimated_used"]),
                    "Qalpha_qualifier": source.get("Qalpha_qualifier", pd.NA),
                    "Talpha_partial_exp_s": talpha,
                    "Talpha_err_plus_s": source.get(
                        "T_alpha_half_err_plus_s", np.nan
                    ),
                    "Talpha_err_minus_s": source.get(
                        "T_alpha_half_err_minus_s", np.nan
                    ),
                    "Talpha_qualifier": source.get(
                        "T_alpha_half_qualifier", pd.NA
                    ),
                    "decay_modes_raw": source.get("decay_modes_raw", pd.NA),
                    "b_alpha_percent": source.get("b_alpha_percent", np.nan),
                    "b_SF_percent": source.get("b_SF_percent", np.nan),
                    "alpha_daughter": source.get("alpha_daughter", pd.NA),
                    "alpha_daughter_population_raw": source.get(
                        "alpha_daughter_population_raw", pd.NA
                    ),
                    "log10_Tcalc_s": log10_tcalc,
                    "Tcalc_s": tcalc,
                    "log10_HF": log10_hf,
                    "HF": hf,
                    "l_used": float(ell) if method.uses_l else np.nan,
                    "l_assumed": method.uses_l,
                    "calculation_status": method_status,
                    "model_status": diagnostics.get(
                        "C_status", "within_implemented_formula"
                    ),
                    "diagnostics_json": json.dumps(
                        diagnostics, ensure_ascii=False, sort_keys=True
                    ),
                    "cell_record_index": int(source["cell_record_index"]),
                    "cell_record_count": int(source["cell_record_count"]),
                    "HF_definition": HF_DEFINITION,
                    "Qalpha_source_column": Q_SOURCE,
                    "Talpha_source_column": TALPHA_SOURCE,
                }
            )

    results = pd.DataFrame(rows)
    expected = int(database["eligible_for_hf"].sum()) * len(methods)
    if len(results) != expected:
        raise AssertionError(f"Expected {expected} results, got {len(results)}")
    if results[["Tcalc_s", "HF", "log10_HF"]].isna().any().any():
        raise AssertionError("NaN found in calculated results")

    expected_hf = results["Talpha_partial_exp_s"] / results["Tcalc_s"]
    if not np.allclose(results["HF"], expected_hf, rtol=2.0e-13, atol=0.0):
        raise AssertionError("HF identity check failed")
    expected_log = np.log10(results["Talpha_partial_exp_s"]) - results[
        "log10_Tcalc_s"
    ]
    if not np.allclose(results["log10_HF"], expected_log, rtol=0.0, atol=2e-13):
        raise AssertionError("log10(HF) identity check failed")

    counts = results.groupby("method")["record_id"].count()
    if counts.nunique() != 1:
        raise AssertionError(f"Method result counts differ: {counts.to_dict()}")
    return results


def build_wide_results(
    database: pd.DataFrame, results: pd.DataFrame, methods: list[MethodSpec]
) -> pd.DataFrame:
    base_columns = [
        "record_id",
        "source_csv_line",
        "isotope",
        "state_label",
        "symbol",
        "A",
        "Z",
        "N",
        "Qalpha_MeV",
        "Ealpha_MeV",
        "T_alpha_half_s",
        "b_alpha_percent",
        "b_SF_percent",
        "alpha_daughter",
        "alpha_daughter_population_raw",
        "Qalpha_estimated",
        "Qalpha_is_experimental",
        "T_alpha_half_qualifier",
        "cell_record_index",
        "cell_record_count",
    ]
    base_columns = [col for col in base_columns if col in database.columns]
    wide = database.loc[database["eligible_for_hf"], base_columns].copy()
    for method in methods:
        subset = results.loc[
            results["method_slug"] == method.slug,
            ["record_id", "Tcalc_s", "log10_Tcalc_s", "HF", "log10_HF"],
        ].rename(
            columns={
                "Tcalc_s": f"Tcalc_s__{method.slug}",
                "log10_Tcalc_s": f"log10_Tcalc_s__{method.slug}",
                "HF": f"HF__{method.slug}",
                "log10_HF": f"log10_HF__{method.slug}",
            }
        )
        wide = wide.merge(subset, on="record_id", how="left", validate="one_to_one")
    return wide


def method_summary(
    results: pd.DataFrame, methods: list[MethodSpec]
) -> pd.DataFrame:
    rows = []
    for method in methods:
        subset = results.loc[results["method_slug"] == method.slug]
        rows.append(
            {
                "method": method.name,
                "method_slug": method.slug,
                "n_records": len(subset),
                "n_unique_nz_cells": subset[["N", "Z"]].drop_duplicates().shape[0],
                "HF_min": subset["HF"].min(),
                "HF_median": subset["HF"].median(),
                "HF_max": subset["HF"].max(),
                "log10_HF_min": subset["log10_HF"].min(),
                "log10_HF_median": subset["log10_HF"].median(),
                "log10_HF_max": subset["log10_HF"].max(),
                "n_HF_lt_0p1": int((subset["HF"] < 0.1).sum()),
                "n_HF_0p1_to_lt4": int(
                    ((subset["HF"] >= 0.1) & (subset["HF"] < 4.0)).sum()
                ),
                "n_HF_4_to_lt10": int(
                    ((subset["HF"] >= 4.0) & (subset["HF"] < 10.0)).sum()
                ),
                "n_HF_10_to_100": int(
                    ((subset["HF"] >= 10.0) & (subset["HF"] <= 100.0)).sum()
                ),
                "n_HF_gt_100": int((subset["HF"] > 100.0).sum()),
                "note": method.note,
            }
        )
    return pd.DataFrame(rows)


def json_value(value: Any) -> Any:
    """Convert pandas/numpy values to strict JSON values."""
    if value is None or value is pd.NA:
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if pd.isna(value):
        return None
    return value


def make_dashboard_payload(
    database: pd.DataFrame,
    results: pd.DataFrame,
    methods: list[MethodSpec],
    summary: pd.DataFrame,
    *,
    ell: float,
    database_path: Path,
) -> dict[str, Any]:
    result_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for _, row in results.iterrows():
        result_lookup[(row["record_id"], row["method_slug"])] = {
            "Tcalc_s": float(row["Tcalc_s"]),
            "log10_Tcalc_s": float(row["log10_Tcalc_s"]),
            "HF": float(row["HF"]),
            "log10_HF": float(row["log10_HF"]),
            "status": row["calculation_status"],
            "model_status": row["model_status"],
        }

    cells: list[dict[str, Any]] = []
    for (Z, N), cell_frame in database.groupby(["Z", "N"], sort=True):
        records = []
        for _, row in cell_frame.sort_values("source_csv_line").iterrows():
            values = {
                method.slug: result_lookup.get((row["record_id"], method.slug))
                for method in methods
            }
            record = {
                "record_id": row["record_id"],
                "source_csv_line": int(row["source_csv_line"]),
                "isotope": json_value(row["isotope"]),
                "state_label": json_value(row.get("state_label")),
                "symbol": json_value(row["symbol"]),
                "A": int(row["A"]),
                "Z": int(row["Z"]),
                "N": int(row["N"]),
                "Qalpha_MeV": json_value(row["Qalpha_MeV"]),
                "Ealpha_MeV": json_value(row.get("Ealpha_MeV")),
                "Qalpha_err_MeV": json_value(row.get("Qalpha_err_MeV")),
                "qalpha_estimated": bool(row["qalpha_estimated_used"]),
                "Qalpha_qualifier": json_value(row.get("Qalpha_qualifier")),
                "Talpha_partial_exp_s": json_value(row["T_alpha_half_s"]),
                "Talpha_err_plus_s": json_value(
                    row.get("T_alpha_half_err_plus_s")
                ),
                "Talpha_err_minus_s": json_value(
                    row.get("T_alpha_half_err_minus_s")
                ),
                "Talpha_qualifier": json_value(
                    row.get("T_alpha_half_qualifier")
                ),
                "T_half_s": json_value(row.get("T_half_s")),
                "b_alpha_percent": json_value(row.get("b_alpha_percent")),
                "b_SF_percent": json_value(row.get("b_SF_percent")),
                "T_SF_half_s": json_value(row.get("T_SF_half_s")),
                "decay_modes_raw": json_value(row.get("decay_modes_raw")),
                "alpha_daughter": json_value(row.get("alpha_daughter")),
                "alpha_daughter_population_raw": json_value(
                    row.get("alpha_daughter_population_raw")
                ),
                "partial_half_life_consistent": json_value(
                    row.get("partial_half_life_consistent")
                ),
                "calculation_status": row["calculation_status"],
                "eligible_for_hf": bool(row["eligible_for_hf"]),
                "cell_record_index": int(row["cell_record_index"]),
                "cell_record_count": int(row["cell_record_count"]),
                "values": values,
            }
            records.append(record)
        cells.append({"key": f"{int(Z)}:{int(N)}", "Z": int(Z), "N": int(N), "records": records})

    method_payload = []
    summary_by_slug = summary.set_index("method_slug")
    for method in methods:
        stats = summary_by_slug.loc[method.slug]
        method_payload.append(
            {
                "id": method.slug,
                "name": method.name,
                "label": method.title,
                "note": method.note,
                "uses_l": method.uses_l,
                "stats": {
                    "n_records": int(stats["n_records"]),
                    "n_unique_nz_cells": int(stats["n_unique_nz_cells"]),
                    "HF_min": float(stats["HF_min"]),
                    "HF_median": float(stats["HF_median"]),
                    "HF_max": float(stats["HF_max"]),
                    "log10_HF_min": float(stats["log10_HF_min"]),
                    "log10_HF_max": float(stats["log10_HF_max"]),
                    "n_HF_lt_0p1": int(stats["n_HF_lt_0p1"]),
                    "n_HF_0p1_to_lt4": int(stats["n_HF_0p1_to_lt4"]),
                    "n_HF_4_to_lt10": int(stats["n_HF_4_to_lt10"]),
                    "n_HF_10_to_100": int(stats["n_HF_10_to_100"]),
                    "n_HF_gt_100": int(stats["n_HF_gt_100"]),
                },
            }
        )

    element_symbols = {
        str(int(z)): str(group["symbol"].dropna().iloc[0])
        for z, group in database.groupby("Z")
        if not group["symbol"].dropna().empty
    }

    return {
        "meta": {
            "database": database_path.name,
            "hf_definition": HF_DEFINITION,
            "talpha_source": TALPHA_SOURCE,
            "qalpha_source": Q_SOURCE,
            "l_assumed": ell,
            "n_database_records": len(database),
            "n_database_cells": int(
                database[["Z", "N"]].drop_duplicates().shape[0]
            ),
            "n_eligible_records": int(database["eligible_for_hf"].sum()),
            "n_eligible_cells": int(
                database.loc[database["eligible_for_hf"], ["Z", "N"]]
                .drop_duplicates()
                .shape[0]
            ),
            "n_excluded_records": int((~database["eligible_for_hf"]).sum()),
            "n_estimated_q_eligible": int(
                (
                    database["eligible_for_hf"]
                    & database["qalpha_estimated_used"]
                ).sum()
            ),
            "n_qualified_talpha_eligible": int(
                (database["eligible_for_hf"] & database["talpha_qualified"]).sum()
            ),
            "N_min": int(database["N"].min()),
            "N_max": int(database["N"].max()),
            "Z_min": int(database["Z"].min()),
            "Z_max": int(database["Z"].max()),
            "color_scale_type": "discrete_HF_bins",
            "color_bins": list(HF_COLOR_BINS),
            "element_symbols": element_symbols,
            "methods": method_payload,
        },
        "cells": cells,
    }


HTML_TEMPLATE = r'''<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>α 衰变阻碍因子 HF 核素图</title>
<style>
:root {
  color-scheme: light;
  --ink: #172033;
  --muted: #5c6577;
  --line: #d9dee8;
  --panel: #ffffff;
  --page: #f4f6fa;
  --accent: #2459a6;
  --warn-bg: #fff4dd;
  --warn-border: #e3a11a;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--page);
  color: var(--ink);
  font-family: Inter, "Segoe UI", "PingFang SC", "Microsoft YaHei", Arial, sans-serif;
  line-height: 1.45;
}
main { width: min(1540px, calc(100% - 28px)); margin: 18px auto 36px; }
header, .toolbar, .chart-shell, .notes { background: var(--panel); border: 1px solid var(--line); border-radius: 12px; }
header { padding: 18px 22px; }
h1 { margin: 0 0 7px; font-size: clamp(1.35rem, 2.3vw, 2rem); }
.formula { margin: 5px 0; font-size: 1.05rem; }
.formula code { padding: 2px 6px; border-radius: 5px; background: #edf2fa; color: #153b75; }
.subtle { color: var(--muted); margin: 4px 0 0; }
.toolbar { margin-top: 12px; padding: 12px 14px; display: flex; align-items: center; flex-wrap: wrap; gap: 10px 16px; }
.toolbar label { display: inline-flex; align-items: center; gap: 7px; font-size: .93rem; }
select, input[type="search"], button {
  min-height: 36px; border: 1px solid #bac3d2; border-radius: 7px; background: white; color: var(--ink); padding: 6px 10px; font: inherit;
}
select { min-width: 235px; }
input[type="search"] { width: 140px; }
button { cursor: pointer; }
button:hover { border-color: var(--accent); color: var(--accent); }
.status { margin-left: auto; color: var(--muted); font-variant-numeric: tabular-nums; }
.warning { display: none; margin-top: 12px; padding: 10px 14px; background: var(--warn-bg); border: 1px solid var(--warn-border); border-radius: 9px; color: #644300; }
.warning.visible { display: block; }
.cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px; margin-top: 12px; }
.card { background: var(--panel); border: 1px solid var(--line); border-radius: 10px; padding: 10px 13px; }
.card .label { display: block; color: var(--muted); font-size: .78rem; }
.card .value { display: block; margin-top: 2px; font-weight: 650; font-size: 1.02rem; font-variant-numeric: tabular-nums; }
.chart-shell { margin-top: 12px; padding: 10px; }
.chart-scroll { overflow: auto; overscroll-behavior: contain; border-radius: 8px; }
#chart { display: block; width: 100%; min-width: 1080px; height: auto; background: #fff; }
.cell { cursor: crosshair; }
.cell-band:focus .cell-border, .cell.highlight .cell-border { stroke: #00a1d6 !important; stroke-width: 2.4 !important; }
.cell-text { pointer-events: none; font-variant-numeric: tabular-nums; text-anchor: middle; dominant-baseline: central; paint-order: stroke; stroke-width: 1.4px; stroke-linejoin: round; }
.axis-label { fill: var(--ink); font-size: 16px; font-weight: 650; }
.tick-label { fill: #3d4658; font-size: 10px; }
.grid-line { stroke: #dfe3ea; stroke-width: .55; }
.plot-border { fill: none; stroke: #8993a4; stroke-width: 1; }
.notes { margin-top: 12px; padding: 14px 18px; display: grid; grid-template-columns: 1.2fr 1fr; gap: 18px; font-size: .91rem; }
.notes h2 { margin: 0 0 7px; font-size: 1rem; }
.notes ul { margin: 5px 0 0 20px; padding: 0; }
.legend-bins { display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 5px; }
.legend-bin { min-width: 0; }
.legend-swatch { height: 16px; border: 1px solid #737b89; border-radius: 3px; }
.legend-label { display: block; margin-top: 3px; color: var(--muted); font-size: .72rem; line-height: 1.2; text-align: center; }
.marker-demo { display: inline-block; width: 18px; height: 13px; border: 1px solid #202634; vertical-align: -2px; margin-right: 5px; background: #eee; }
.marker-demo.dashed { border-style: dashed; background: #fdb863; }
#tooltip {
  position: fixed; z-index: 20; width: min(360px, calc(100vw - 24px)); pointer-events: none; background: rgba(20,27,41,.97); color: #fff; border-radius: 9px; padding: 10px 12px; box-shadow: 0 7px 30px rgba(0,0,0,.22); font-size: .84rem;
}
#tooltip[hidden] { display: none; }
#tooltip .tip-title { font-weight: 700; font-size: 1rem; margin-bottom: 5px; }
#tooltip dl { display: grid; grid-template-columns: max-content 1fr; gap: 2px 10px; margin: 0; }
#tooltip dt { color: #bfc8d8; }
#tooltip dd { margin: 0; text-align: right; font-variant-numeric: tabular-nums; overflow-wrap: anywhere; }
@media (max-width: 850px) {
  .cards { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .notes { grid-template-columns: 1fr; }
  .status { width: 100%; margin-left: 0; }
}
@media (max-width: 520px) {
  main { width: calc(100% - 16px); margin-top: 8px; }
  header { padding: 14px; }
  .cards { grid-template-columns: 1fr; }
  .toolbar > label { width: 100%; justify-content: space-between; }
  select, input[type="search"] { flex: 1; }
}
@media print {
  body { background: white; }
  main { width: 100%; margin: 0; }
  .toolbar, #tooltip { display: none !important; }
  header, .chart-shell, .notes, .card { border-color: #bbb; box-shadow: none; }
  .chart-scroll { overflow: visible; }
  #chart { min-width: 0; width: 100%; }
}
</style>
</head>
<body>
<main>
  <header>
    <h1>α 衰变阻碍因子 HF 核素图</h1>
    <p class="formula"><code>HF = Tα,partial(exp) / Tα,calc</code></p>
    <p class="subtle">分子直接取自数据库字段 <strong>T_alpha_half_s</strong>（α 部分半衰期），不使用总半衰期 T_half_s，也不通过 α 分支比重算。</p>
  </header>

  <section class="toolbar" aria-label="图表控制">
    <label>计算方法 <select id="method-select"></select></label>
    <label><input id="show-labels" type="checkbox" checked>显示 HF 数值</label>
    <label><input id="show-na" type="checkbox" checked>显示不可计算核素</label>
    <label>查找核素 <input id="search" type="search" placeholder="如 273Ds"></label>
    <button id="reset" type="button">复位</button>
    <output id="status" class="status" aria-live="polite"></output>
  </section>
  <aside id="method-warning" class="warning" role="note"></aside>

  <section class="cards" aria-label="方法统计">
    <div class="card"><span class="label">当前方法</span><span id="card-method" class="value">—</span></div>
    <div class="card"><span class="label">可计算记录 / 核素格</span><span id="card-count" class="value">—</span></div>
    <div class="card"><span class="label">HF 中位数</span><span id="card-median" class="value">—</span></div>
    <div class="card"><span class="label">HF 范围</span><span id="card-range" class="value">—</span></div>
  </section>

  <section class="chart-shell">
    <div class="chart-scroll">
      <svg id="chart" role="img" aria-labelledby="chart-title chart-desc"></svg>
    </div>
  </section>

  <section class="notes">
    <div>
      <h2>数据与质量标记</h2>
      <ul>
        <li><span class="marker-demo dashed"></span>虚线边框和“*”：Qα 为估算值。</li>
        <li>HF 前的 ≈、&lt;、&gt;、≤、≥ 直接继承 T_alpha_half_s 的限定符。</li>
        <li>灰色“—”：数据库中存在该核素，但缺少正的 Qα 或 α 部分半衰期，HF 不定义。</li>
        <li>同一 N,Z 的多条记录拆成上下子格，未做平均。</li>
      </ul>
    </div>
    <div>
      <h2>共享颜色尺度</h2>
      <div id="legend-bins" class="legend-bins" aria-label="HF 离散分档色标"></div>
      <p id="color-note" class="subtle"></p>
    </div>
  </section>
</main>
<div id="tooltip" role="tooltip" hidden></div>
<script id="hf-data" type="application/json">__HF_PAYLOAD__</script>
<script>
"use strict";
const DATA = JSON.parse(document.getElementById("hf-data").textContent);
const SVG_NS = "http://www.w3.org/2000/svg";
const chart = document.getElementById("chart");
const methodSelect = document.getElementById("method-select");
const showLabels = document.getElementById("show-labels");
const showNA = document.getElementById("show-na");
const searchInput = document.getElementById("search");
const tooltip = document.getElementById("tooltip");
const warning = document.getElementById("method-warning");
const statusOutput = document.getElementById("status");

const CELL = 31;
const MARGIN = {left: 72, right: 34, top: 32, bottom: 66};
const nCount = DATA.meta.N_max - DATA.meta.N_min + 1;
const zCount = DATA.meta.Z_max - DATA.meta.Z_min + 1;
const WIDTH = MARGIN.left + MARGIN.right + nCount * CELL;
const HEIGHT = MARGIN.top + MARGIN.bottom + zCount * CELL;
chart.setAttribute("viewBox", `0 0 ${WIDTH} ${HEIGHT}`);
chart.setAttribute("preserveAspectRatio", "xMidYMid meet");

function svgEl(name, attrs = {}, text = null) {
  const el = document.createElementNS(SVG_NS, name);
  for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, String(value));
  if (text !== null) el.textContent = text;
  return el;
}

function hexRgb(hex) {
  const value = hex.replace("#", "");
  return [0,2,4].map(i => parseInt(value.slice(i, i + 2), 16));
}

function colorFor(hf) {
  if (!Number.isFinite(hf)) return "#e5e7eb";
  for (const bin of DATA.meta.color_bins) {
    const aboveMin = bin.min === null || hf >= bin.min;
    const belowMax = bin.max === null || hf < bin.max || (bin.include_max && hf <= bin.max);
    if (aboveMin && belowMax) return bin.color;
  }
  return "#e5e7eb";
}

function textColor(fill) {
  const srgb = hexRgb(fill).map(v => v / 255).map(v => v <= .04045 ? v / 12.92 : ((v + .055) / 1.055) ** 2.4);
  const lum = .2126 * srgb[0] + .7152 * srgb[1] + .0722 * srgb[2];
  const blackContrast = (lum + .05) / .05;
  const whiteContrast = 1.05 / (lum + .05);
  return blackContrast >= whiteContrast ? "#111827" : "#ffffff";
}

function compactNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "—";
  if (value === 0) return "0";
  const abs = Math.abs(value);
  if (abs < 1e-2 || abs >= 1e3) {
    const [mantissa, exponent] = value.toExponential(digits - 1).split("e");
    return `${Number(mantissa)}e${Number(exponent)}`;
  }
  return Number(value.toPrecision(digits)).toString();
}

function preciseNumber(value) {
  if (!Number.isFinite(value)) return "—";
  return value.toExponential(6).replace("e+", "e");
}

function qualifierPrefix(value) {
  const q = (value || "").trim();
  return ({"approx":"≈",">":">","<":"<",">=":"≥","<=":"≤","uncertain":"?"})[q] || "";
}

function hfLabel(record, result) {
  if (!result) return "—";
  const prefix = qualifierPrefix(record.Talpha_qualifier);
  const suffix = record.qalpha_estimated ? "*" : "";
  return `${prefix}${compactNumber(result.HF, 2)}${suffix}`;
}

function displayIsotope(name) {
  return String(name || "")
    .replace("^a", "ᵃ")
    .replace("^b", "ᵇ")
    .replace("^m", "ᵐ");
}

function stateTag(record, records) {
  if (records.length === 1) return "";
  const stateMatch = String(record.isotope).match(/\^([abm])$/);
  if (stateMatch) return stateMatch[1];
  if (records.some(item => String(item.isotope).includes("^m"))) return "g";
  return String(record.cell_record_index);
}

function drawAxes() {
  chart.replaceChildren();
  chart.append(svgEl("title", {id:"chart-title"}, "HF 核素图：横轴 N，纵轴 Z"));
  chart.append(svgEl("desc", {id:"chart-desc"}, "颜色按 HF 数值离散分档，格内文字表示 HF；可切换七种计算方法。"));
  const axes = svgEl("g", {id:"axes"});
  for (let i = 0; i <= nCount; i++) {
    const x = MARGIN.left + i * CELL;
    axes.append(svgEl("line", {x1:x,y1:MARGIN.top,x2:x,y2:MARGIN.top+zCount*CELL,class:"grid-line"}));
  }
  for (let i = 0; i <= zCount; i++) {
    const y = MARGIN.top + i * CELL;
    axes.append(svgEl("line", {x1:MARGIN.left,y1:y,x2:MARGIN.left+nCount*CELL,y2:y,class:"grid-line"}));
  }
  for (let N = DATA.meta.N_min; N <= DATA.meta.N_max; N++) {
    if ((N - DATA.meta.N_min) % 2 !== 0 && N !== DATA.meta.N_max) continue;
    const x = MARGIN.left + (N - DATA.meta.N_min + .5) * CELL;
    axes.append(svgEl("text", {x,y:HEIGHT-MARGIN.bottom+20,class:"tick-label","text-anchor":"middle"}, N));
  }
  for (let Z = DATA.meta.Z_min; Z <= DATA.meta.Z_max; Z++) {
    const y = MARGIN.top + (DATA.meta.Z_max - Z + .5) * CELL + 3;
    const symbol = DATA.meta.element_symbols[String(Z)] || "";
    axes.append(svgEl("text", {x:MARGIN.left-7,y,class:"tick-label","text-anchor":"end"}, `${Z} ${symbol}`));
  }
  axes.append(svgEl("rect", {x:MARGIN.left,y:MARGIN.top,width:nCount*CELL,height:zCount*CELL,class:"plot-border"}));
  axes.append(svgEl("text", {x:MARGIN.left+nCount*CELL/2,y:HEIGHT-13,class:"axis-label","text-anchor":"middle"}, "Neutron number N"));
  axes.append(svgEl("text", {x:17,y:MARGIN.top+zCount*CELL/2,class:"axis-label","text-anchor":"middle",transform:`rotate(-90 17 ${MARGIN.top+zCount*CELL/2})`}, "Proton number Z"));
  chart.append(axes);
  chart.append(svgEl("g", {id:"cells"}));
}

const recordMap = new Map();
for (const cell of DATA.cells) for (const record of cell.records) recordMap.set(record.record_id, record);

function render() {
  const methodId = methodSelect.value;
  const method = DATA.meta.methods.find(item => item.id === methodId);
  const cellsLayer = document.getElementById("cells");
  const fragment = document.createDocumentFragment();
  let visibleResults = 0;
  let visibleCells = 0;

  for (const cell of DATA.cells) {
    const validInCell = cell.records.filter(record => record.values[methodId]);
    if (!showNA.checked && validInCell.length === 0) continue;
    visibleCells++;
    const x = MARGIN.left + (cell.N - DATA.meta.N_min) * CELL + .7;
    const y = MARGIN.top + (DATA.meta.Z_max - cell.Z) * CELL + .7;
    const inner = CELL - 1.4;
    const outer = svgEl("g", {class:"cell", "data-cell-key":cell.key});
    const bandHeight = inner / cell.records.length;

    cell.records.forEach((record, index) => {
      const result = record.values[methodId];
      if (result) visibleResults++;
      const bandY = y + index * bandHeight;
      const fill = result ? colorFor(result.HF) : "#e5e7eb";
      const strokeDash = record.qalpha_estimated && result ? "3 2" : "none";
      const band = svgEl("g", {
        class:"cell-band",
        tabindex:"0",
        role:"button",
        "data-record-id":record.record_id,
        "aria-label":`${displayIsotope(record.isotope)}, N=${record.N}, Z=${record.Z}, ${method.label}, HF=${result ? preciseNumber(result.HF) : "不可计算"}`
      });
      band.append(svgEl("rect", {x,y:bandY,width:inner,height:bandHeight,fill,stroke:"#262b36","stroke-width":.75,"stroke-dasharray":strokeDash,class:"cell-border"}));
      if (record.qalpha_estimated && result) {
        band.append(svgEl("path", {d:`M ${x+inner-7} ${bandY} L ${x+inner} ${bandY} L ${x+inner} ${bandY+7} Z`,fill:"#111827","pointer-events":"none"}));
      }
      if (showLabels.checked) {
        const color = textColor(fill);
        if (cell.records.length === 1) {
          band.append(svgEl("text", {x:x+inner/2,y:bandY+bandHeight*.31,fill:color,"font-size":"7.4",class:"cell-text",stroke:color === "#ffffff" ? "#111827" : "#ffffff","stroke-opacity":.45}, displayIsotope(record.isotope)));
          band.append(svgEl("text", {x:x+inner/2,y:bandY+bandHeight*.69,fill:color,"font-size":"8.2","font-weight":"700",class:"cell-text",stroke:color === "#ffffff" ? "#111827" : "#ffffff","stroke-opacity":.45}, hfLabel(record,result)));
        } else {
          band.append(svgEl("text", {x:x+inner/2,y:bandY+bandHeight/2,fill:color,"font-size":"6.9","font-weight":"700",class:"cell-text",stroke:color === "#ffffff" ? "#111827" : "#ffffff","stroke-opacity":.45}, `${stateTag(record,cell.records)} ${hfLabel(record,result)}`));
        }
      }
      outer.append(band);
    });
    fragment.append(outer);
  }
  cellsLayer.replaceChildren(fragment);

  const stats = method.stats;
  document.getElementById("card-method").textContent = method.label;
  document.getElementById("card-count").textContent = `${stats.n_records} / ${stats.n_unique_nz_cells}`;
  document.getElementById("card-median").textContent = compactNumber(stats.HF_median, 4);
  document.getElementById("card-range").textContent = `${compactNumber(stats.HF_min,3)} – ${compactNumber(stats.HF_max,3)}`;
  statusOutput.textContent = `当前显示 ${visibleResults} 条 HF、${visibleCells} 个 N–Z 格点`;

  const warningParts = [];
  if (method.uses_l) warningParts.push(`数据库无跃迁角动量 l，本方法使用 l=${DATA.meta.l_assumed} 假设。`);
  if (method.id === "xu2022_unified") warningParts.push("全部可计算记录均在本地 Xu 实现的 C(Z,N) 定义区之外，按脚本规则使用 C=0 外推。");
  warning.textContent = warningParts.join(" ");
  warning.classList.toggle("visible", warningParts.length > 0);
  applySearch();
}

function addTipRow(dl, term, value) {
  dl.append(svgFree("dt", term));
  dl.append(svgFree("dd", value));
}

function svgFree(tag, text) {
  const el = document.createElement(tag);
  el.textContent = text;
  return el;
}

function tooltipContent(record) {
  const method = DATA.meta.methods.find(item => item.id === methodSelect.value);
  const result = record.values[method.id];
  tooltip.replaceChildren();
  tooltip.append(svgFree("div", displayIsotope(record.isotope)));
  tooltip.firstChild.className = "tip-title";
  const dl = document.createElement("dl");
  addTipRow(dl,"核素",`A=${record.A}, Z=${record.Z}, N=${record.N}`);
  addTipRow(dl,"方法",method.label);
  if (record.Ealpha_MeV !== null) addTipRow(dl,"Eα (MeV)",compactNumber(record.Ealpha_MeV,6));
  addTipRow(dl,"Qα (MeV)",record.Qalpha_MeV === null ? "—" : `${compactNumber(record.Qalpha_MeV,6)}${record.qalpha_estimated ? "  (estimated)" : ""}`);
  if (record.decay_modes_raw) addTipRow(dl,"衰变方式",record.decay_modes_raw);
  if (record.b_alpha_percent !== null) addTipRow(dl,"BRα",`${compactNumber(record.b_alpha_percent,6)}%`);
  if (record.b_SF_percent !== null) addTipRow(dl,"BRSF",`${compactNumber(record.b_SF_percent,6)}%`);
  if (record.alpha_daughter) addTipRow(dl,"α 子核",record.alpha_daughter);
  if (record.alpha_daughter_population_raw) addTipRow(dl,"子核态布居",record.alpha_daughter_population_raw);
  if (record.T_half_s !== null) addTipRow(dl,"T1/2,total (s)",preciseNumber(record.T_half_s));
  addTipRow(dl,"Tα,partial(exp) (s)",record.Talpha_partial_exp_s === null ? "—" : `${preciseNumber(record.Talpha_partial_exp_s)}${record.Talpha_qualifier ? `  [${record.Talpha_qualifier}]` : ""}`);
  addTipRow(dl,"Tα 数据字段",DATA.meta.talpha_source);
  if (result) {
    addTipRow(dl,"Tα,calc (s)",preciseNumber(result.Tcalc_s));
    addTipRow(dl,"HF",`${qualifierPrefix(record.Talpha_qualifier)}${preciseNumber(result.HF)}${record.qalpha_estimated ? " *" : ""}`);
    addTipRow(dl,"log10(HF)",compactNumber(result.log10_HF,6));
    addTipRow(dl,"模型状态",result.model_status);
  } else {
    addTipRow(dl,"HF","不可计算");
    addTipRow(dl,"原因",record.calculation_status);
  }
  if (method.uses_l) addTipRow(dl,"角动量",`l=${DATA.meta.l_assumed} (assumed)`);
  addTipRow(dl,"数据库行",String(record.source_csv_line));
  tooltip.append(dl);
}

function positionTooltip(event) {
  const gap = 14;
  const box = tooltip.getBoundingClientRect();
  let left = event.clientX + gap;
  let top = event.clientY + gap;
  if (left + box.width > window.innerWidth - 8) left = event.clientX - box.width - gap;
  if (top + box.height > window.innerHeight - 8) top = event.clientY - box.height - gap;
  tooltip.style.left = `${Math.max(8,left)}px`;
  tooltip.style.top = `${Math.max(8,top)}px`;
}

let pinned = false;
chart.addEventListener("pointerover", event => {
  if (pinned) return;
  const band = event.target.closest("[data-record-id]");
  if (!band) return;
  tooltipContent(recordMap.get(band.dataset.recordId));
  tooltip.hidden = false;
  positionTooltip(event);
});
chart.addEventListener("pointermove", event => { if (!tooltip.hidden && !pinned) positionTooltip(event); });
chart.addEventListener("pointerout", event => { if (!pinned && !event.relatedTarget?.closest?.("[data-record-id]")) tooltip.hidden = true; });
chart.addEventListener("click", event => {
  const band = event.target.closest("[data-record-id]");
  if (!band) { pinned = false; tooltip.hidden = true; return; }
  pinned = true;
  tooltipContent(recordMap.get(band.dataset.recordId));
  tooltip.hidden = false;
  positionTooltip(event);
});
chart.addEventListener("focusin", event => {
  const band = event.target.closest("[data-record-id]");
  if (!band) return;
  tooltipContent(recordMap.get(band.dataset.recordId));
  tooltip.hidden = false;
  const rect = band.getBoundingClientRect();
  positionTooltip({clientX:rect.right,clientY:rect.top});
});
chart.addEventListener("focusout", () => { if (!pinned) tooltip.hidden = true; });
document.addEventListener("keydown", event => { if (event.key === "Escape") { pinned = false; tooltip.hidden = true; } });

function applySearch() {
  const query = searchInput.value.trim().toLowerCase().replace("ᵐ","^m");
  document.querySelectorAll(".cell.highlight").forEach(node => node.classList.remove("highlight"));
  if (!query) return;
  for (const cell of DATA.cells) {
    if (cell.records.some(record => String(record.isotope).toLowerCase().includes(query))) {
      const node = chart.querySelector(`[data-cell-key="${cell.key}"]`);
      if (node) node.classList.add("highlight");
    }
  }
}

for (const method of DATA.meta.methods) {
  const option = document.createElement("option");
  option.value = method.id;
  option.textContent = method.label;
  methodSelect.append(option);
}
const legendBins = document.getElementById("legend-bins");
for (const bin of DATA.meta.color_bins) {
  const item = document.createElement("div");
  item.className = "legend-bin";
  const swatch = document.createElement("div");
  swatch.className = "legend-swatch";
  swatch.style.background = bin.color;
  const label = document.createElement("span");
  label.className = "legend-label";
  label.textContent = bin.label;
  item.append(swatch, label);
  legendBins.append(item);
}
document.getElementById("color-note").textContent = `颜色直接按 HF 固定分档，同一档内不再连续渐变；格内和悬停仍显示未截断的原始 HF。DZR、Ismail 和 Xu 中数据库没有 l，统一采用 l=${DATA.meta.l_assumed} 假设。`;
methodSelect.addEventListener("change", () => requestAnimationFrame(render));
showLabels.addEventListener("change", render);
showNA.addEventListener("change", render);
searchInput.addEventListener("input", applySearch);
document.getElementById("reset").addEventListener("click", () => {
  methodSelect.selectedIndex = 0;
  showLabels.checked = true;
  showNA.checked = true;
  searchInput.value = "";
  pinned = false;
  tooltip.hidden = true;
  render();
});

drawAxes();
render();
</script>
</body>
</html>
'''


def write_dashboard(payload: dict[str, Any], output: Path) -> None:
    json_text = json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    )
    # Prevent a source value from terminating the application/json script tag.
    json_text = json_text.replace("<", "\\u003c")
    html = HTML_TEMPLATE.replace("__HF_PAYLOAD__", json_text)
    output.write_text(html, encoding="utf-8")


def qualifier_prefix(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return {
        "approx": "≈",
        ">": ">",
        "<": "<",
        ">=": "≥",
        "<=": "≤",
        "uncertain": "?",
    }.get(str(value).strip(), "")


def compact_hf(value: float) -> str:
    if not math.isfinite(value):
        return "—"
    text = f"{value:.2g}"
    return (
        text.replace("e+0", "e")
        .replace("e-0", "e-")
        .replace("e+", "e")
    )


def hf_color(value: float) -> str:
    """Return the fixed discrete color for an HF value."""
    if value < 0.1:
        return str(HF_COLOR_BINS[0]["color"])
    if value < 4.0:
        return str(HF_COLOR_BINS[1]["color"])
    if value < 10.0:
        return str(HF_COLOR_BINS[2]["color"])
    if value <= 100.0:
        return str(HF_COLOR_BINS[3]["color"])
    return str(HF_COLOR_BINS[4]["color"])


def text_color_for_rgba(rgba: Any) -> str:
    if isinstance(rgba, str):
        value = rgba.removeprefix("#")
        rgba = tuple(int(value[index : index + 2], 16) / 255 for index in (0, 2, 4)) + (1.0,)

    def linear(channel: float) -> float:
        return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4

    luminance = 0.2126 * linear(rgba[0]) + 0.7152 * linear(rgba[1]) + 0.0722 * linear(rgba[2])
    black_contrast = (luminance + 0.05) / 0.05
    white_contrast = 1.05 / (luminance + 0.05)
    return "black" if black_contrast >= white_contrast else "white"


def write_static_overview(
    database: pd.DataFrame,
    results: pd.DataFrame,
    methods: list[MethodSpec],
    output_dir: Path,
    *,
    dpi: int,
) -> None:
    """Write a publication-friendly, non-interactive overview as PNG/PDF."""
    cache_dir = Path(tempfile.gettempdir()) / "hfdistr-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    n_min, n_max = int(database["N"].min()), int(database["N"].max())
    z_min, z_max = int(database["Z"].min()), int(database["Z"].max())

    fig, axes = plt.subplots(4, 2, figsize=(18, 21), constrained_layout=True)
    flat_axes = axes.ravel()
    for ax, method in zip(flat_axes[:7], methods):
        method_results = results.loc[
            results["method_slug"] == method.slug,
            ["record_id", "HF", "log10_HF"],
        ]
        plotted = database.merge(
            method_results, on="record_id", how="left", validate="one_to_one"
        )
        for (Z, N), group in plotted.groupby(["Z", "N"], sort=False):
            group = group.sort_values("source_csv_line")
            count = len(group)
            for band_index, (_, row) in enumerate(group.iterrows()):
                height = 0.92 / count
                y0 = Z - 0.46 + band_index * height
                valid = pd.notna(row["HF"])
                face = hf_color(float(row["HF"])) if valid else "#e5e7eb"
                linestyle = "--" if valid and bool(row["qalpha_estimated_used"]) else "-"
                rect = Rectangle(
                    (N - 0.46, y0),
                    0.92,
                    height,
                    facecolor=face,
                    edgecolor="#20242c",
                    linewidth=0.45,
                    linestyle=linestyle,
                )
                ax.add_patch(rect)
                if valid:
                    prefix = qualifier_prefix(row.get("T_alpha_half_qualifier"))
                    suffix = "*" if bool(row["qalpha_estimated_used"]) else ""
                    label = f"{prefix}{compact_hf(float(row['HF']))}{suffix}"
                else:
                    label = "—"
                if count > 1:
                    state_match = re.search(r"\^([abm])$", str(row["isotope"]))
                    if state_match:
                        tag = state_match.group(1)
                    elif group["isotope"].astype(str).str.contains(r"\^m").any():
                        tag = "g"
                    else:
                        tag = str(int(row["cell_record_index"]))
                    label = f"{tag} {label}"
                color = text_color_for_rgba(face) if valid else "#5d6470"
                ax.text(
                    N,
                    y0 + height / 2,
                    label,
                    ha="center",
                    va="center",
                    fontsize=3.7 if count == 1 else 3.1,
                    fontweight="semibold",
                    color=color,
                    clip_on=True,
                )

        ax.set_xlim(n_min - 0.7, n_max + 0.7)
        ax.set_ylim(z_min - 0.7, z_max + 0.7)
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(n_min, n_max + 1, 4))
        ax.set_yticks(np.arange(z_min, z_max + 1, 2))
        ax.set_xticks(np.arange(n_min - 0.5, n_max + 1.5, 1), minor=True)
        ax.set_yticks(np.arange(z_min - 0.5, z_max + 1.5, 1), minor=True)
        ax.grid(which="minor", color="#dfe3ea", linewidth=0.35)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.set_xlabel("Neutron number N")
        ax.set_ylabel("Proton number Z")
        ax.set_title(method.title)

    info_ax = flat_axes[7]
    info_ax.axis("off")
    info_ax.text(
        0.10,
        0.93,
        "Shared discrete HF colors",
        transform=info_ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
    )
    for index, color_bin in enumerate(HF_COLOR_BINS):
        y = 0.83 - index * 0.085
        info_ax.add_patch(
            Rectangle(
                (0.10, y),
                0.12,
                0.055,
                transform=info_ax.transAxes,
                facecolor=str(color_bin["color"]),
                edgecolor="#30343b",
                linewidth=0.6,
            )
        )
        info_ax.text(
            0.25,
            y + 0.0275,
            str(color_bin["label"]),
            transform=info_ax.transAxes,
            fontsize=9,
            va="center",
        )
    info_ax.text(
        0.10,
        0.36,
        "Color = fixed HF interval\n"
        "Text = un-clipped HF\n"
        "* / dashed edge = estimated Qα\n"
        "— = missing Qα or Tα,partial(exp)\n"
        "Split cell = multiple records at one N,Z\n\n"
        "HF numerator: database T_alpha_half_s\n"
        "l = 0 assumed where required\n"
        "Xu: C=0 outside-domain extrapolation",
        transform=info_ax.transAxes,
        fontsize=9,
        va="top",
        linespacing=1.5,
    )
    fig.suptitle(
        "Alpha-decay hindrance factors across the nuclide chart",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.002,
        "HF = Tα,partial(exp) / Tα,calc; Tα,partial(exp) is read directly from T_alpha_half_s.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(
        output_dir / "hf_nuclide_chart_overview.png",
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        output_dir / "hf_nuclide_chart_overview.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def write_metadata(
    output: Path,
    database: pd.DataFrame,
    methods: list[MethodSpec],
    *,
    database_path: Path,
    ell: float,
    summary: pd.DataFrame,
) -> None:
    metadata = {
        "database": str(database_path),
        "database_records": len(database),
        "database_unique_nz_cells": int(
            database[["Z", "N"]].drop_duplicates().shape[0]
        ),
        "eligible_records": int(database["eligible_for_hf"].sum()),
        "eligible_unique_nz_cells": int(
            database.loc[database["eligible_for_hf"], ["Z", "N"]]
            .drop_duplicates()
            .shape[0]
        ),
        "excluded_records": int((~database["eligible_for_hf"]).sum()),
        "hf_definition": HF_DEFINITION,
        "experimental_partial_half_life_source": TALPHA_SOURCE,
        "total_half_life_used_in_hf": False,
        "branching_ratio_used_to_reconstruct_talpha": False,
        "qalpha_source": Q_SOURCE,
        "l_assumed": ell,
        "color_scale": {
            "quantity": "HF",
            "type": "discrete_bins",
            "bins": list(HF_COLOR_BINS),
            "text_values_clipped": False,
        },
        "methods": [
            {
                "name": method.name,
                "slug": method.slug,
                "title": method.title,
                "note": method.note,
                "uses_assumed_l": method.uses_l,
            }
            for method in methods
        ],
        "method_summary": summary.to_dict(orient="records"),
    }
    output.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate HF from Qalpha_MeV and the experimental partial alpha "
            "half-life T_alpha_half_s, then build a self-contained HTML dashboard."
        )
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=DEFAULT_DATABASE,
        help=f"Source CSV database (default: {DEFAULT_DATABASE})",
    )
    parser.add_argument(
        "--method-dir",
        type=Path,
        default=DEFAULT_METHOD_DIR,
        help=f"Directory containing the seven method scripts (default: {DEFAULT_METHOD_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--l-assumed",
        type=float,
        default=0.0,
        help="Assumed alpha orbital angular momentum where required (default: 0)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Static overview PNG resolution (default: 300)",
    )
    parser.add_argument(
        "--experimental-q-only",
        action="store_true",
        help="Exclude rows whose Qalpha is marked nonexperimental/estimated",
    )
    parser.add_argument(
        "--no-static-overview",
        action="store_true",
        help="Skip the PNG/PDF overview and build only data files plus HTML",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    database_path = args.database.resolve()
    method_dir = args.method_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not database_path.is_file():
        raise FileNotFoundError(f"Database not found: {database_path}")
    if not method_dir.is_dir():
        raise FileNotFoundError(f"Method directory not found: {method_dir}")
    if not math.isfinite(args.l_assumed) or args.l_assumed < 0.0:
        raise ValueError("--l-assumed must be finite and non-negative")
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72")

    output_dir.mkdir(parents=True, exist_ok=True)
    methods = load_methods(method_dir)
    database = prepare_database(
        database_path,
        experimental_q_only=args.experimental_q_only,
        ell=args.l_assumed,
    )
    audit = audit_database_consistency(database)
    results = evaluate_methods(database, methods, args.l_assumed)
    wide = build_wide_results(database, results, methods)
    summary = method_summary(results, methods)

    quality_columns = [
        "record_id",
        "source_csv_line",
        "isotope",
        "state_label",
        "symbol",
        "A",
        "Z",
        "N",
        "is_parenthesized",
        "mass_number_uncertain",
        "Qalpha_MeV",
        "Ealpha_MeV",
        "Qalpha_err_MeV",
        "Qalpha_estimated",
        "Qalpha_is_experimental",
        "qalpha_estimated_used",
        "Qalpha_qualifier",
        "T_alpha_half_s",
        "T_alpha_half_err_plus_s",
        "T_alpha_half_err_minus_s",
        "T_alpha_half_qualifier",
        "T_half_s",
        "b_alpha_percent",
        "b_SF_percent",
        "T_SF_half_s",
        "decay_modes_raw",
        "alpha_daughter",
        "alpha_daughter_population_raw",
        "partial_half_life_reconstructed_s",
        "partial_half_life_consistent",
        "nonexperimental_fields",
        "eligible_for_hf",
        "calculation_status",
        "cell_record_index",
        "cell_record_count",
        "duplicate_nz_cell",
        "l_assumed",
    ]
    quality_columns = [col for col in quality_columns if col in database.columns]
    quality = database[quality_columns].copy()

    results.to_csv(
        output_dir / "hf_results_long.csv",
        index=False,
        float_format="%.12g",
    )
    wide.to_csv(
        output_dir / "hf_results_wide.csv",
        index=False,
        float_format="%.12g",
    )
    quality.to_csv(
        output_dir / "hf_data_quality.csv",
        index=False,
        float_format="%.12g",
    )
    quality.loc[~quality["eligible_for_hf"]].to_csv(
        output_dir / "hf_excluded_records.csv",
        index=False,
        float_format="%.12g",
    )
    summary.to_csv(
        output_dir / "hf_method_summary.csv",
        index=False,
        float_format="%.12g",
    )
    audit.to_csv(
        output_dir / "hf_consistency_audit.csv",
        index=False,
        float_format="%.12g",
    )

    payload = make_dashboard_payload(
        database,
        results,
        methods,
        summary,
        ell=args.l_assumed,
        database_path=database_path,
    )
    write_dashboard(payload, output_dir / "hf_nuclide_chart.html")
    write_metadata(
        output_dir / "run_metadata.json",
        database,
        methods,
        database_path=database_path,
        ell=args.l_assumed,
        summary=summary,
    )
    if not args.no_static_overview:
        write_static_overview(
            database,
            results,
            methods,
            output_dir,
            dpi=args.dpi,
        )

    print(f"Database records: {len(database)}")
    print(
        "HF-calculable records: "
        f"{int(database['eligible_for_hf'].sum())} "
        f"({database.loc[database['eligible_for_hf'], ['Z', 'N']].drop_duplicates().shape[0]} N-Z cells)"
    )
    print(f"Excluded records: {int((~database['eligible_for_hf']).sum())}")
    print(f"Methods: {len(methods)}; calculated rows: {len(results)}")
    print(
        "Consistency audit: "
        f"{int((audit['severity'] == 'error').sum())} error(s), "
        f"{int((audit['severity'] == 'warning').sum())} warning(s)"
    )
    print(f"HF definition: {HF_DEFINITION}")
    print(f"Dashboard: {output_dir / 'hf_nuclide_chart.html'}")


if __name__ == "__main__":
    main()
