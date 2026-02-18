#!/usr/bin/env python3
"""
compute_proton_Araw_from_worldparam.py

UPDATED: Now reads Pbeam, PHe3, Pp and their uncertainties from the kin CSV.

Required kin CSV columns (per-row):
  Q2, eps_bar, tau_bar, Px_bar, Pz_bar,
  Pbeam, dPbeam, PHe3, dPHe3, Pp, dPp

(You can still pass --Pp, etc., but they only act as fallbacks if a row is missing a value.)

See previous header for physics and error-propagation details.
"""

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from typing import List, Dict, Tuple

MU_P = 2.79284734462  # proton magnetic moment


# -----------------------------
# Utilities
# -----------------------------
def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False


def read_proton_lookup(path: str) -> List[Tuple[float, float, float, float, float, float, float]]:
    """
    Read lookup table with 7 columns:
      Q2  GE/GD  dGE/GD  dGEpar/GD  GM/(mu GD)  dGM/(mu GD)  dGMpar/(mu GD)
    Lines beginning with # are ignored.
    """
    rows = []
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 7:
                continue
            if not all(_is_float(x) for x in parts[:7]):
                continue
            vals = tuple(float(x) for x in parts[:7])
            rows.append(vals)

    if not rows:
        raise RuntimeError(f"No valid rows read from lookup file: {path}")

    rows.sort(key=lambda r: r[0])  # Ensure increasing Q2
    return rows


def lin_interp(x: float, xp: List[float], fp: List[float]) -> float:
    """Simple linear interpolation with clamping at ends."""
    if x <= xp[0]:
        return fp[0]
    if x >= xp[-1]:
        return fp[-1]

    lo, hi = 0, len(xp) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xp[mid] <= x:
            lo = mid
        else:
            hi = mid

    x0, x1 = xp[lo], xp[hi]
    y0, y1 = fp[lo], fp[hi]
    if x1 == x0:
        return y0
    t = (x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


@dataclass
class LambdaResult:
    Lam: float
    dLam_stat: float
    dLam_par: float
    dLam_tot: float


def lambda_from_lookup(Q2: float, lookup_rows) -> LambdaResult:
    q2 = [r[0] for r in lookup_rows]
    ge = [r[1] for r in lookup_rows]   # GE/GD
    dge = [r[2] for r in lookup_rows]
    dgepar = [r[3] for r in lookup_rows]

    gm = [r[4] for r in lookup_rows]   # GM/(mu GD)
    dgm = [r[5] for r in lookup_rows]
    dgmpar = [r[6] for r in lookup_rows]

    GEo = lin_interp(Q2, q2, ge)
    dGEo = lin_interp(Q2, q2, dge)
    dGEo_par = lin_interp(Q2, q2, dgepar)

    GMo = lin_interp(Q2, q2, gm)
    dGMo = lin_interp(Q2, q2, dgm)
    dGMo_par = lin_interp(Q2, q2, dgmpar)

    if GMo == 0.0:
        raise ZeroDivisionError(f"Interpolated GM/(mu GD) == 0 at Q2={Q2}")

    Lam = GEo / (MU_P * GMo)

    def prop(dGEo_i: float, dGMo_i: float) -> float:
        rel2 = 0.0
        if GEo != 0.0:
            rel2 += (dGEo_i / GEo) ** 2
        if GMo != 0.0:
            rel2 += (dGMo_i / GMo) ** 2
        return abs(Lam) * math.sqrt(rel2)

    dLam_stat = prop(dGEo, dGMo)
    dLam_par = prop(dGEo_par, dGMo_par)
    dLam_tot = math.sqrt(dLam_stat**2 + dLam_par**2)

    return LambdaResult(Lam=Lam, dLam_stat=dLam_stat, dLam_par=dLam_par, dLam_tot=dLam_tot)


@dataclass
class ApResult:
    Ap_phys: float
    dAp_phys: float
    Ap: float
    dAp: float


def compute_Ap(Q2: float,
               eps_bar: float, tau_bar: float, Px_bar: float, Pz_bar: float,
               Lam: float, dLam: float,
               Pbeam: float, dPbeam: float,
               PHe3: float, dPHe3: float,
               Pp: float, dPp: float) -> ApResult:
    """
    Implements Eq. 5.35, 5.36 and simplified error propagation (holding eps,tau,Px,Pz constants)
    """
    if tau_bar <= 0.0:
        raise ValueError(f"tau_bar must be >0 (got {tau_bar}) at Q2={Q2}")
    if eps_bar == 0.0:
        raise ValueError(f"eps_bar must be != 0 (got {eps_bar}) at Q2={Q2}")

    a = math.sqrt(2.0 * eps_bar * (1.0 - eps_bar) / tau_bar)
    b = math.sqrt(max(0.0, 1.0 - eps_bar**2))
    D = 1.0 + (tau_bar / eps_bar) * (Lam**2)

    Ap_phys = -(1.0 / D) * (Lam * a * Px_bar + b * Pz_bar)

    # d(Ap_phys)/dLambda in your Eq. 5.39 form
    dApphys_dLam = -(1.0 / D) * (a * Px_bar + 2.0 * Ap_phys * (tau_bar / eps_bar) * Lam)

    dAp_phys = abs(dApphys_dLam) * dLam

    Ap = Pbeam * PHe3 * Pp * Ap_phys

    # relative uncertainty combination (Eq. 5.40 style)
    rel2 = 0.0
    if Ap_phys != 0.0:
        rel2 += (dAp_phys / Ap_phys) ** 2
    if Pbeam != 0.0:
        rel2 += (dPbeam / Pbeam) ** 2
    if PHe3 != 0.0:
        rel2 += (dPHe3 / PHe3) ** 2
    if Pp != 0.0:
        rel2 += (dPp / Pp) ** 2

    dAp = abs(Ap) * math.sqrt(rel2) if rel2 > 0 else 0.0

    return ApResult(Ap_phys=Ap_phys, dAp_phys=dAp_phys, Ap=Ap, dAp=dAp)


def read_kin_csv(path: str) -> List[Dict[str, float]]:
    """
    Read kin CSV with REQUIRED columns:
      Q2, eps_bar, tau_bar, Px_bar, Pz_bar,
      Pbeam, dPbeam, PHe3, dPHe3, Pp, dPp
    """
    required = [
        "Q2", "eps_bar", "tau_bar", "Px_bar", "Pz_bar",
        "Pbeam", "dPbeam", "PHe3", "dPHe3", "Pp", "dPp"
    ]
    rows_out: List[Dict[str, float]] = []

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in {path}")
        fields = [x.strip() for x in reader.fieldnames]

        for r in required:
            if r not in fields:
                raise RuntimeError(
                    f"Missing required column '{r}' in {path}.\n"
                    f"Found columns: {fields}\n"
                    f"Expected at least: {required}"
                )

        for i, row in enumerate(reader, start=2):  # start=2 => header is line 1
            d: Dict[str, float] = {}
            for k in required:
                v = row.get(k, "")
                v = "" if v is None else v.strip()
                if v == "":
                    raise RuntimeError(f"Empty value for required column '{k}' at line {i} in {path}")
                d[k] = float(v)
            rows_out.append(d)

    if not rows_out:
        raise RuntimeError(f"No data rows found in {path}")
    return rows_out


def write_out_csv(path: str, rows: List[Dict[str, float]], field_order: List[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=field_order)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in field_order})


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compute proton Ap (raw asymmetry for correction) from world-data lookup parameterization."
    )
    ap.add_argument("--lookup", required=True, help="Path to proton_lookup.dat")
    ap.add_argument("--kin", required=True, help="CSV with required kinematics + polarization columns")
    ap.add_argument("--out", default="Ap_proton_out.csv", help="Output CSV")

    # These are now true fallbacks only (kept for convenience / backward compat)
    ap.add_argument("--Pbeam", type=float, default=float("nan"))
    ap.add_argument("--dPbeam", type=float, default=float("nan"))
    ap.add_argument("--PHe3", type=float, default=float("nan"))
    ap.add_argument("--dPHe3", type=float, default=float("nan"))
    ap.add_argument("--Pp", type=float, default=float("nan"))
    ap.add_argument("--dPp", type=float, default=float("nan"))

    ap.add_argument("--use_total_lambda_unc", action="store_true",
                    help="Use total dLambda = sqrt(stat^2 + par^2). Default: use stat-only.")
    ap.add_argument("--print", action="store_true", help="Also print a human-readable summary to stdout.")

    args = ap.parse_args()

    lookup = read_proton_lookup(args.lookup)
    kin_rows = read_kin_csv(args.kin)

    out_rows: List[Dict[str, float]] = []

    for row in kin_rows:
        Q2 = row["Q2"]
        eps_bar = row["eps_bar"]
        tau_bar = row["tau_bar"]
        Px_bar = row["Px_bar"]
        Pz_bar = row["Pz_bar"]

        lam = lambda_from_lookup(Q2, lookup)
        dLam_used = lam.dLam_tot if args.use_total_lambda_unc else lam.dLam_stat

        # now read per-row pols from CSV (fallback to args only if NaN provided)
        def pick(val_row, val_arg):
            if val_row is not None:
                return val_row
            return val_arg

        Pbeam = row["Pbeam"] if "Pbeam" in row else args.Pbeam
        dPbeam = row["dPbeam"] if "dPbeam" in row else args.dPbeam
        PHe3 = row["PHe3"] if "PHe3" in row else args.PHe3
        dPHe3 = row["dPHe3"] if "dPHe3" in row else args.dPHe3
        Pp = row["Pp"] if "Pp" in row else args.Pp
        dPp = row["dPp"] if "dPp" in row else args.dPp

        ap_res = compute_Ap(
            Q2=Q2,
            eps_bar=eps_bar, tau_bar=tau_bar, Px_bar=Px_bar, Pz_bar=Pz_bar,
            Lam=lam.Lam, dLam=dLam_used,
            Pbeam=Pbeam, dPbeam=dPbeam,
            PHe3=PHe3, dPHe3=dPHe3,
            Pp=Pp, dPp=dPp
        )

        out_rows.append({
            "Q2": Q2,
            "eps_bar": eps_bar,
            "tau_bar": tau_bar,
            "Px_bar": Px_bar,
            "Pz_bar": Pz_bar,

            "Lambda": lam.Lam,
            "dLambda_stat": lam.dLam_stat,
            "dLambda_par": lam.dLam_par,
            "dLambda_tot": lam.dLam_tot,
            "dLambda_used": dLam_used,

            "Ap_phys": ap_res.Ap_phys,
            "dAp_phys": ap_res.dAp_phys,

            "Pbeam": Pbeam,
            "dPbeam": dPbeam,
            "PHe3": PHe3,
            "dPHe3": dPHe3,
            "Pp": Pp,
            "dPp": dPp,

            "Ap": ap_res.Ap,
            "dAp": ap_res.dAp,
        })

    field_order = [
        "Q2", "eps_bar", "tau_bar", "Px_bar", "Pz_bar",
        "Lambda", "dLambda_stat", "dLambda_par", "dLambda_tot", "dLambda_used",
        "Ap_phys", "dAp_phys",
        "Pbeam", "dPbeam", "PHe3", "dPHe3", "Pp", "dPp",
        "Ap", "dAp"
    ]

    write_out_csv(args.out, out_rows, field_order)

    if args.print:
        print("\n=== Proton asymmetry results ===")
        for r in out_rows:
            print(
                f"Q2={r['Q2']:.4g}  "
                f"Lambda={r['Lambda']:.6g} ± {r['dLambda_used']:.3g}  "
                f"Ap_phys={r['Ap_phys']:.6g} ± {r['dAp_phys']:.3g}  "
                f"Ap={r['Ap']:.6g} ± {r['dAp']:.3g}"
            )
        print(f"\nWrote: {args.out}\n")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)

