#!/usr/bin/env python3
"""
scan_cutsets.py

Wrapper script to scan cut sets for the GEN final analysis.

What it does
------------
1. Reads a base cuts JSON.
2. Builds a grid of modified cut values.
3. Writes one temporary cuts JSON per trial.
4. Runs your analysis executable for each trial.
5. Parses output text files from txt/ and corrections/.
6. Copies the outputs for each trial into a unique directory.
7. Writes a summary CSV and a ranked CSV.

Important
---------
- This runs SERIALly.
- Do not run in parallel unless you first modify the C++ code so each trial writes
  unique output filenames/directories internally.
- This script assumes your executable overwrites the same standard output files
  under txt/, corrections/, rootfiles/, images/, plots/.

Example
-------
python3 scan_cutsets.py \
  --project-dir /path/to/Final_analysis \
  --exe ./final_analysis \
  --data /path/to/data.root \
  --sim-qe /path/to/sim_qe.root \
  --sim-pim /path/to/sim_pim.root \
  --sim-inel /path/to/sim_inel.root \
  --sim-n2 /path/to/sim_n2.root \
  --base-cuts cuts/GEN3.json \
  --run-quality DB/run_quality.json \
  --kin GEN3_He3 \
  --scan W2_H=1.10,1.15,1.20,1.25 \
  --scan eHCAL_L=0.08,0.10,0.12 \
  --scan coin_L=113,114,115 \
  --scan coin_H=121,122,123 \
  --scan dx_r=0.18,0.20,0.22 \
  --scan dy_r=0.28,0.30,0.32 \
  --metric total_err

Alternative scan syntax
-----------------------
You can also use:
  --scan W2_H=1.10:1.30:0.05
which means start:stop:step (inclusive if floating-point allows).

You can also scan a single value:
  --scan dx_c=0.0
"""

import argparse
import csv
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ----------------------------------------------------------------------
# Utility helpers
# ----------------------------------------------------------------------

def parse_scan_spec(spec: str) -> Tuple[str, List[float]]:
    """
    Parse a scan spec like:
      W2_H=1.1,1.2,1.3
      W2_H=1.1:1.3:0.05
      dx_c=0.0

    Returns:
      ("W2_H", [1.1, 1.2, 1.3])
    """
    if "=" not in spec:
        raise ValueError(f"Bad --scan spec '{spec}'. Expected KEY=values")

    key, rhs = spec.split("=", 1)
    key = key.strip()
    rhs = rhs.strip()

    if not key:
        raise ValueError(f"Bad --scan spec '{spec}'. Empty key")

    # comma-separated
    if "," in rhs:
        vals = [float(x.strip()) for x in rhs.split(",") if x.strip()]
        if not vals:
            raise ValueError(f"Bad --scan spec '{spec}'. No values found")
        return key, vals

    # range start:stop:step
    if rhs.count(":") == 2:
        start_s, stop_s, step_s = rhs.split(":")
        start = float(start_s)
        stop = float(stop_s)
        step = float(step_s)
        if step == 0:
            raise ValueError(f"Bad --scan spec '{spec}'. Step cannot be zero")

        vals = []
        x = start

        # robust inclusive stepping
        if step > 0:
            while x <= stop + 1e-12:
                vals.append(round(x, 12))
                x += step
        else:
            while x >= stop - 1e-12:
                vals.append(round(x, 12))
                x += step

        if not vals:
            raise ValueError(f"Bad --scan spec '{spec}'. Empty range")
        return key, vals

    # single value
    return key, [float(rhs)]


def safe_float_from_regex(text: str, pattern: str) -> Optional[float]:
    m = re.search(pattern, text, re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def read_text_if_exists(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(errors="replace")


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def copy_any(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    if src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        mkdir(dst.parent)
        shutil.copy2(src, dst)


def format_value_for_tag(v: float) -> str:
    """
    Make float safe for folder/file tags.
    Example:
      1.20 -> 1p2
      -0.05 -> m0p05
    """
    s = f"{v:.12g}"
    s = s.replace("-", "m").replace(".", "p")
    return s


def flatten_trial_tag(cuts: Dict[str, float], ordered_keys: List[str], idx: int) -> str:
    pieces = [f"{idx:04d}"]
    for k in ordered_keys:
        pieces.append(f"{k}_{format_value_for_tag(cuts[k])}")
    return "__".join(pieces)


def json_dump_pretty(obj: dict, path: Path) -> None:
    with path.open("w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


# ----------------------------------------------------------------------
# Parsing output files
# ----------------------------------------------------------------------

def parse_summary_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "Aphys":            safe_float_from_regex(text, r"^Aphys\s*=\s*([^\s]+)"),
        "err_Aphys_stat":   safe_float_from_regex(text, r"^err_Aphys_stat\s*=\s*([^\s]+)"),
        "err_Aphys_sys":    safe_float_from_regex(text, r"^err_Aphys_sys\s*=\s*([^\s]+)"),
        "err_Aacc_sys":     safe_float_from_regex(text, r"^err_Aacc_sys\s*=\s*([^\s]+)"),
        "err_Api_sys":      safe_float_from_regex(text, r"^err_Api_sys\s*=\s*([^\s]+)"),
        "err_Ain_sys":      safe_float_from_regex(text, r"^err_Ain_sys\s*=\s*([^\s]+)"),
        "err_Ap_sys":       safe_float_from_regex(text, r"^err_Ap_sys\s*=\s*([^\s]+)"),
        "err_facc_sys":     safe_float_from_regex(text, r"^err_facc_sys\s*=\s*([^\s]+)"),
        "err_fpi_sys":      safe_float_from_regex(text, r"^err_fpi_sys\s*=\s*([^\s]+)"),
        "err_fin_sys":      safe_float_from_regex(text, r"^err_fin_sys\s*=\s*([^\s]+)"),
        "err_fp_sys":       safe_float_from_regex(text, r"^err_fp_sys\s*=\s*([^\s]+)"),
        "err_Ptar_sys":     safe_float_from_regex(text, r"^err_Ptar_sys\s*=\s*([^\s]+)"),
        "err_Pn_sys":       safe_float_from_regex(text, r"^err_Pn_sys_%\s*=\s*([^\s]+)"),
        "err_Pbeam_sys":    safe_float_from_regex(text, r"^err_Pbeam_sys_%\s*=\s*([^\s]+)"),
    }


def parse_raw_asym_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "Np":        safe_float_from_regex(text, r"^Np\s*=\s*([^\s]+)"),
        "Nm":        safe_float_from_regex(text, r"^Nm\s*=\s*([^\s]+)"),
        "A_raw":     safe_float_from_regex(text, r"^A_raw\s*=\s*([^\s]+)"),
        "err_A_raw": safe_float_from_regex(text, r"^err_A_raw\s*=\s*([^\s]+)"),
    }


def parse_ehcal_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "lambda":            safe_float_from_regex(text, r"^lambda\s*=\s*([^\s]+)"),
        "sigma_lambda_stat": safe_float_from_regex(text, r"^sigma_lambda_stat\s*=\s*([^\s]+)"),
        "sigma_lambda_sys":  safe_float_from_regex(text, r"^sigma_lambda_sys\s*=\s*([^\s]+)"),
        "sigma_lambda":      safe_float_from_regex(text, r"^sigma_lambda\s*=\s*([^\s]+)"),
    }


def parse_inelastic_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "background_fraction":   safe_float_from_regex(text, r"^background_fraction\s*=\s*([^\s]+)"),
        "err_background_fraction": safe_float_from_regex(text, r"^err_inelastic_fraction\s*=\s*([^\s]+)"),
        "proton_fraction":       safe_float_from_regex(text, r"^proton_fraction\s*=\s*([^\s]+)"),
        "err_proton_fraction":   safe_float_from_regex(text, r"^err_proton_fraction\s*=\s*([^\s]+)"),
        "f_in":                  safe_float_from_regex(text, r"^f_in\s*=\s*([^\s]+)"),
        "err_f_in":              safe_float_from_regex(text, r"^err_f_in\s*=\s*([^\s]+)"),
        "f_p":                   safe_float_from_regex(text, r"^f_p\s*=\s*([^\s]+)"),
        "err_f_p":               safe_float_from_regex(text, r"^err_f_p\s*=\s*([^\s]+)"),
        "A_in":                  safe_float_from_regex(text, r"^A_in\s*=\s*([^\s]+)"),
        "err_A_in":              safe_float_from_regex(text, r"^err_A_in\s*=\s*([^\s]+)"),
        "A_p":                   safe_float_from_regex(text, r"^A_p\s*=\s*([^\s]+)"),
        "err_A_p":               safe_float_from_regex(text, r"^err_A_p\s*=\s*([^\s]+)"),
    }


def parse_pion_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "A_pi":    safe_float_from_regex(text, r"^A_pi\s*=\s*([^\s]+)"),
        "err_A_pi": safe_float_from_regex(text, r"^err_A_pi\s*=\s*([^\s]+)"),
        "f_pi":    safe_float_from_regex(text, r"^f_pi\s*=\s*([^\s]+)"),
        "err_f_pi": safe_float_from_regex(text, r"^err_f_pi\s*=\s*([^\s]+)"),
    }


def parse_accidental_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "N_plus_all":        safe_float_from_regex(text, r"^N_plus_all\s*=\s*([^\s]+)"),
        "N_minus_all":       safe_float_from_regex(text, r"^N_minus_all\s*=\s*([^\s]+)"),
        "A_acc":             safe_float_from_regex(text, r"^A_acc\s*=\s*([^\s]+)"),
        "err_A_acc":         safe_float_from_regex(text, r"^err_A_acc\s*=\s*([^\s]+)"),
        "accidental_events": safe_float_from_regex(text, r"^accidental_events\s*=\s*([^\s]+)"),
        "QE_events":         safe_float_from_regex(text, r"^QE_events\s*=\s*([^\s]+)"),
        "f_acc":             safe_float_from_regex(text, r"^f_acc\s*=\s*([^\s]+)"),
        "err_f_acc":         safe_float_from_regex(text, r"^err_f_acc\s*=\s*([^\s]+)"),
    }


def parse_nitrogen_file(path: Path) -> Dict[str, Optional[float]]:
    text = read_text_if_exists(path)
    return {
        "f_N2":     safe_float_from_regex(text, r"^f_N2\s*=\s*([^\s]+)"),
        "err_f_N2": safe_float_from_regex(text, r"^err_f_N2\s*=\s*([^\s]+)"),
    }


def parse_all_outputs(project_dir: Path, kin: str) -> Dict[str, Optional[float]]:
    txt_dir = project_dir / "txt"
    corr_dir = project_dir / "corrections"

    data = {}
    data.update(parse_summary_file(txt_dir / f"physics_neutron_asymmetry_summary_{kin}.txt"))
    data.update(parse_raw_asym_file(txt_dir / f"raw_asym_{kin}.txt"))
    data.update(parse_ehcal_file(txt_dir / f"results_eHCAL_cut_{kin}.txt"))
    data.update(parse_inelastic_file(corr_dir / f"InelasticCorrection_{kin}.txt"))
    data.update(parse_pion_file(corr_dir / f"PionCorrection_{kin}.txt"))
    data.update(parse_accidental_file(corr_dir / f"AccidentalCorrection_{kin}.txt"))
    data.update(parse_nitrogen_file(corr_dir / f"NitrogenCorrection_{kin}.txt"))

    # Derived totals
    stat = data.get("err_Aphys_stat")
    sys_ = data.get("err_Aphys_sys")
    if stat is not None and sys_ is not None:
        data["total_err"] = math.sqrt(stat * stat + sys_ * sys_)
    else:
        data["total_err"] = None

    # Approx neutron purity
    # f_n ~ 1 - (f_acc + f_pi + f_in + f_p + f_N2)
    f_acc = data.get("f_acc") or 0.0
    f_pi  = data.get("f_pi") or 0.0
    f_in  = data.get("f_in") or 0.0
    f_p   = data.get("f_p") or 0.0
    f_N2  = data.get("f_N2") or 0.0
    data["f_n_approx"] = 1.0 - (f_acc + f_pi + f_in + f_p + f_N2)

    # A simple figure of merit you may inspect later
    if data["f_n_approx"] is not None and stat not in (None, 0):
        data["fom_fn2_over_stat2"] = (data["f_n_approx"] ** 2) / (stat ** 2)
    else:
        data["fom_fn2_over_stat2"] = None

    return data


# ----------------------------------------------------------------------
# Analysis execution
# ----------------------------------------------------------------------

def run_analysis(
    project_dir: Path,
    exe: str,
    data_file: str,
    sim_qe: str,
    sim_pim: str,
    sim_inel: str,
    sim_n2: str,
    cuts_file: str,
    run_quality: Optional[str],
    kin: str,
    log_path: Path,
    dry_run: bool = False,
) -> int:
    cmd = [
        exe,
        data_file,
        sim_qe,
        sim_pim,
        sim_inel,
        sim_n2,
        cuts_file,
    ]

    if run_quality:
        cmd.append(run_quality)

    cmd.append(kin)

    with log_path.open("w") as logf:
        logf.write("COMMAND:\n")
        logf.write(" ".join(cmd) + "\n\n")

        if dry_run:
            logf.write("[dry-run] analysis not executed\n")
            return 0

        proc = subprocess.run(
            cmd,
            cwd=str(project_dir),
            stdout=logf,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return proc.returncode


# ----------------------------------------------------------------------
# Ranking
# ----------------------------------------------------------------------

def metric_value(row: Dict[str, object], metric: str) -> float:
    """
    Lower is better for:
      stat_err, sys_err, total_err
    Higher is better for:
      fom_fn2_over_stat2, f_n_approx
    """
    lower_better = {"stat_err", "sys_err", "total_err"}
    higher_better = {"fom_fn2_over_stat2", "f_n_approx"}

    if metric not in lower_better and metric not in higher_better:
        raise ValueError(f"Unknown metric '{metric}'")

    value = row.get(metric)

    if value is None:
        return float("inf") if metric in lower_better else float("-inf")

    try:
        value = float(value)
    except Exception:
        return float("inf") if metric in lower_better else float("-inf")

    if metric in lower_better:
        return value

    return -value


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Scan cut sets for Final_analysis")
    ap.add_argument("--project-dir", required=True, help="Path to Final_analysis directory")
    ap.add_argument("--exe", default="./final_analysis", help="Executable path relative to project-dir, or absolute path")

    ap.add_argument("--data", required=True, help="Data ROOT file")
    ap.add_argument("--sim-qe", required=True, help="QE sim ROOT file")
    ap.add_argument("--sim-pim", required=True, help="Pion sim ROOT file")
    ap.add_argument("--sim-inel", required=True, help="Inelastic sim ROOT file")
    ap.add_argument("--sim-n2", required=True, help="N2 sim ROOT file")

    ap.add_argument("--base-cuts", required=True, help="Base cuts JSON, relative to project-dir or absolute")
    ap.add_argument("--run-quality", default=None, help="Run quality JSON/CSV")
    ap.add_argument("--kin", required=True, help="Kinematic label, e.g. GEN3_He3")

    ap.add_argument(
        "--scan",
        action="append",
        default=[],
        help="Scan specification, e.g. W2_H=1.1,1.2,1.3 or W2_H=1.1:1.3:0.05",
    )

    ap.add_argument(
        "--metric",
        default="total_err",
        choices=["stat_err", "sys_err", "total_err", "fom_fn2_over_stat2", "f_n_approx"],
        help="Ranking metric",
    )

    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory for scan results. Default: <project-dir>/scan_results/<timestamp>_<kin>",
    )

    ap.add_argument(
        "--copy",
        nargs="*",
        default=["txt", "corrections"],
        help="Directories/files under project-dir to snapshot per trial. Example: txt corrections rootfiles plots",
    )

    ap.add_argument("--resume", action="store_true", help="Skip trials whose done.json already exists")
    ap.add_argument("--dry-run", action="store_true", help="Generate cut files and commands, but do not run analysis")

    args = ap.parse_args()

    project_dir = Path(args.project_dir).resolve()
    if not project_dir.exists():
        print(f"ERROR: project dir does not exist: {project_dir}", file=sys.stderr)
        return 2

    exe = args.exe
    if not os.path.isabs(exe):
        exe_path = (project_dir / exe).resolve()
    else:
        exe_path = Path(exe).resolve()

    if not args.dry_run and not exe_path.exists():
        print(f"ERROR: executable not found: {exe_path}", file=sys.stderr)
        return 2

    base_cuts_path = Path(args.base_cuts)
    if not base_cuts_path.is_absolute():
        base_cuts_path = (project_dir / base_cuts_path).resolve()
    if not base_cuts_path.exists():
        print(f"ERROR: base cuts file not found: {base_cuts_path}", file=sys.stderr)
        return 2

    run_quality_path = None
    if args.run_quality:
        run_quality_path = Path(args.run_quality)
        if not run_quality_path.is_absolute():
            run_quality_path = (project_dir / run_quality_path).resolve()
        if not run_quality_path.exists():
            print(f"ERROR: run quality file not found: {run_quality_path}", file=sys.stderr)
            return 2

    with base_cuts_path.open() as f:
        base_cuts = json.load(f)

    scan_map: Dict[str, List[float]] = {}
    for spec in args.scan:
        key, vals = parse_scan_spec(spec)
        scan_map[key] = vals

    if not scan_map:
        print("ERROR: no scans defined. Use at least one --scan KEY=...", file=sys.stderr)
        return 2

    ordered_keys = list(scan_map.keys())
    value_lists = [scan_map[k] for k in ordered_keys]
    grid = list(itertools.product(*value_lists))

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    if args.outdir is None:
        outdir = project_dir / "scan_results" / f"{timestamp}_{args.kin}"
    else:
        outdir = Path(args.outdir).resolve()

    cuts_dir = outdir / "cuts"
    trials_dir = outdir / "trials"
    mkdir(cuts_dir)
    mkdir(trials_dir)

    print(f"Project dir : {project_dir}")
    print(f"Executable  : {exe_path}")
    print(f"Base cuts   : {base_cuts_path}")
    print(f"Kinematic   : {args.kin}")
    print(f"Trials      : {len(grid)}")
    print(f"Output dir  : {outdir}")
    print()

    rows: List[Dict[str, object]] = []

    for idx, combo in enumerate(grid, start=1):
        cut_updates = {k: v for k, v in zip(ordered_keys, combo)}
        tag = flatten_trial_tag(cut_updates, ordered_keys, idx)

        trial_dir = trials_dir / tag
        mkdir(trial_dir)

        done_json = trial_dir / "done.json"
        if args.resume and done_json.exists():
            print(f"[{idx:4d}/{len(grid)}] resume-skip {tag}")
            with done_json.open() as f:
                rows.append(json.load(f))
            continue

        # Build cuts JSON for this trial
        cuts_this = dict(base_cuts)
        cuts_this.update(cut_updates)

        cuts_json_path = cuts_dir / f"{tag}.json"
        json_dump_pretty(cuts_this, cuts_json_path)

        # Run analysis
        print(f"[{idx:4d}/{len(grid)}] running {tag}")
        log_path = trial_dir / "run.log"

        rc = run_analysis(
            project_dir=project_dir,
            exe=str(exe_path),
            data_file=args.data,
            sim_qe=args.sim_qe,
            sim_pim=args.sim_pim,
            sim_inel=args.sim_inel,
            sim_n2=args.sim_n2,
            cuts_file=str(cuts_json_path),
            run_quality=str(run_quality_path) if run_quality_path else None,
            kin=args.kin,
            log_path=log_path,
            dry_run=args.dry_run,
        )

        # Parse outputs
        parsed = parse_all_outputs(project_dir, args.kin)

        # Rename metric column names to simpler names for ranking
        row: Dict[str, object] = {
            "trial_index": idx,
            "tag": tag,
            "return_code": rc,
            "cuts_json": str(cuts_json_path),
            "trial_dir": str(trial_dir),
        }

        # save scanned variables
        for k in ordered_keys:
            row[k] = cut_updates[k]

        # save parsed outputs
        row["Aphys"] = parsed.get("Aphys")
        row["stat_err"] = parsed.get("err_Aphys_stat")
        row["sys_err"] = parsed.get("err_Aphys_sys")
        row["total_err"] = parsed.get("total_err")

        row["A_raw"] = parsed.get("A_raw")
        row["err_A_raw"] = parsed.get("err_A_raw")
        row["Np"] = parsed.get("Np")
        row["Nm"] = parsed.get("Nm")

        row["lambda"] = parsed.get("lambda")
        row["sigma_lambda_stat"] = parsed.get("sigma_lambda_stat")
        row["sigma_lambda_sys"] = parsed.get("sigma_lambda_sys")
        row["sigma_lambda"] = parsed.get("sigma_lambda")

        row["f_acc"] = parsed.get("f_acc")
        row["err_f_acc"] = parsed.get("err_f_acc")
        row["A_acc"] = parsed.get("A_acc")
        row["err_A_acc"] = parsed.get("err_A_acc")

        row["f_pi"] = parsed.get("f_pi")
        row["err_f_pi"] = parsed.get("err_f_pi")
        row["A_pi"] = parsed.get("A_pi")
        row["err_A_pi"] = parsed.get("err_A_pi")

        row["f_in"] = parsed.get("f_in")
        row["err_f_in"] = parsed.get("err_f_in")
        row["A_in"] = parsed.get("A_in")
        row["err_A_in"] = parsed.get("err_A_in")

        row["f_p"] = parsed.get("f_p")
        row["err_f_p"] = parsed.get("err_f_p")
        row["A_p"] = parsed.get("A_p")
        row["err_A_p"] = parsed.get("err_A_p")

        row["f_N2"] = parsed.get("f_N2")
        row["err_f_N2"] = parsed.get("err_f_N2")

        row["background_fraction"] = parsed.get("background_fraction")
        row["proton_fraction"] = parsed.get("proton_fraction")

        row["f_n_approx"] = parsed.get("f_n_approx")
        row["fom_fn2_over_stat2"] = parsed.get("fom_fn2_over_stat2")

        # snapshot outputs after run
        if rc == 0 or args.dry_run:
            for item in args.copy:
                src = project_dir / item
                dst = trial_dir / item
                copy_any(src, dst)

        # save row json for resume
        with done_json.open("w") as f:
            json.dump(row, f, indent=2, sort_keys=True)

        rows.append(row)

        if rc != 0:
            print(f"  WARNING: return code = {rc}")

    if not rows:
        print("No rows produced.", file=sys.stderr)
        return 1

    # Write full CSV
    all_keys = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                all_keys.append(k)

    csv_all = outdir / "scan_results_all.csv"
    with csv_all.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(rows)

    # Ranked CSV
    ranked_rows = sorted(rows, key=lambda r: metric_value(r, args.metric))
    csv_ranked = outdir / f"scan_results_ranked_by_{args.metric}.csv"
    with csv_ranked.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(ranked_rows)

    # Print top results
    print()
    print(f"Wrote: {csv_all}")
    print(f"Wrote: {csv_ranked}")
    print()
    print(f"Top 10 by {args.metric}:")
    for i, row in enumerate(ranked_rows[:10], start=1):
        print(
            f"{i:2d}. tag={row.get('tag')}  "
            f"{args.metric}={row.get(args.metric)}  "
            f"Aphys={row.get('Aphys')}  "
            f"stat={row.get('stat_err')}  "
            f"sys={row.get('sys_err')}  "
            f"f_n~={row.get('f_n_approx')}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
