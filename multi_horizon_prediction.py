#!/usr/bin/env python3
"""
multi_horizon_prediction.py — how far ahead can co-resistance exposure predict?

Two questions this script answers, both extensions of stage_prediction() in
reproduce_paper.py:

  (1) HORIZON.  Train ONCE on all data up to a single cutoff year t0, then from
      that same fixed model predict N years into the future (target year
      t0 + N) for N = 1..5 and compare against the real resistance observed in
      year t0 + N.  The training history — and therefore the co-resistance
      network and the rho features — is IDENTICAL for every horizon; the only
      thing that changes is how far out the target year sits.  This isolates
      the pure effect of forecast distance and shows how predictability decays.

  (2) TESTED-IN-THE-FUTURE PAIRS.  A candidate is scored ONLY if the
      pathogen-antibiotic pair is actually tested (n_tested >= 30) in the
      future target year t0 + N.  Untested future pairs are meaningless in
      practice (there is no measurement to compare against), so they are never
      counted as right or wrong.  The script reports, per horizon, how many
      predicted pairs are dropped for lack of a future measurement.

INPUT  : same species x year x antibiotic summary CSV as the paper.
OUTPUT : <datadir>/horizon_prediction_auc.csv        one row per horizon
         <datadir>/horizon_prediction_candidates.csv  every scored pair
         <figdir>/fig6_horizon.{pdf,png}              decay curve

Usage:
    python multi_horizon_prediction.py
    python multi_horizon_prediction.py --cutoff 2019 --horizons 1 2 3 4 5
"""
import argparse, os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------- config / IO
HERE = os.path.dirname(os.path.abspath(__file__))
DEF_INPUT = f"{HERE}/results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv"

ap = argparse.ArgumentParser(description="Single-cutoff multi-horizon resistance prediction, evaluated only on future-tested pairs.")
ap.add_argument("--input", default=DEF_INPUT, help="species x year x antibiotic summary CSV")
ap.add_argument("--outdir", default=HERE, help="base output dir (default: project root)")
ap.add_argument("--figdir", default=None, help="override figure dir (default: <outdir>/paper/figures)")
ap.add_argument("--datadir", default=None, help="override data dir (default: <outdir>/results/analysis)")
ap.add_argument("--cutoff", type=int, default=2019, help="training cutoff year t0 (history = 2004..t0); default 2019 so N=5 -> 2024")
ap.add_argument("--horizons", type=int, nargs="+", default=[1, 2, 3, 4, 5], help="forecast horizons N (years ahead of t0)")
args = ap.parse_args()

FIG = args.figdir or f"{args.outdir}/paper/figures"
DATA = args.datadir or f"{args.outdir}/results/analysis"
os.makedirs(FIG, exist_ok=True); os.makedirs(DATA, exist_ok=True)

MIN = 30            # min tested isolates per cell (matches the summary file's own filter)
HIST_START = 2004   # expanding history always begins here
T0 = args.cutoff

print(f"input   : {args.input}")
print(f"cutoff  : train on {HIST_START}..{T0}  (same history for every horizon)")
print(f"horizons: {args.horizons}  -> target years {[T0 + N for N in sorted(args.horizons)]}")
print(f"data    : {DATA}")
print(f"figures : {FIG}\n")
df = pd.read_csv(args.input)
YMAX = int(df.Year.max())

# ---------------------------------------------------------------- core helpers
# (identical construction to reproduce_paper.py stage_prediction, so the N=1
#  horizon reproduces the paper's rolling-origin numbers on this cutoff.)
def win(y0, y1):
    """Pooled resistance-fraction matrix (species x antibiotic) over [y0,y1]."""
    d = df[(df.Year >= y0) & (df.Year <= y1)]
    g = d.groupby(["Species", "Antibiotic"])[["n_R", "n_tested"]].sum().reset_index()
    g = g[g.n_tested >= MIN]; g["prop_R"] = g.n_R / g.n_tested
    return g.pivot(index="Species", columns="Antibiotic", values="prop_R")

def rca_val(p):
    tot = np.nansum(p.values); rp = p.sum(1); ra = p.sum(0); exp = np.outer(rp, ra) / tot
    return pd.DataFrame(p.values / np.where(exp > 0, exp, np.nan), index=p.index, columns=p.columns)

def rca_bin(p, cut=1.0):
    return (rca_val(p) >= cut)

def cooc_cols(R):
    """Co-occurrence projection onto antibiotics: W_ab = # species resistant to both a and b."""
    ab = list(R.columns); A = R.fillna(0).astype(float).to_numpy().T
    W = A @ A.T; np.fill_diagonal(W, 0.0); return pd.DataFrame(W, index=ab, columns=ab)

def cooc_rows(R):
    """Co-occurrence projection onto pathogens: W_pq = # antibiotics resisted by both p and q."""
    sp = list(R.index); A = R.fillna(0).astype(float).to_numpy()
    W = A @ A.T; np.fill_diagonal(W, 0.0); return pd.DataFrame(W, index=sp, columns=sp)

def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float); n1 = y.sum(); n0 = len(y) - n1
    if n1 < 3 or n0 < 3: return np.nan
    r = pd.Series(s).rank().values; return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

# ---------------------------------------------------------------- train ONCE at t0
# History, binary matrix, co-resistance networks and rho features are computed a
# single time; every horizon reuses them and only swaps the future target year.
pr_h = win(HIST_START, T0); Rh = rca_bin(pr_h)
Waa = cooc_cols(Rh); Was = Waa.sum(1)
Wpp = cooc_rows(Rh); Wps = Wpp.sum(1)
pop = Rh.mean(0); AB = list(Rh.columns); SP = list(Rh.index); Rh0 = Rh.astype(float)

# candidate pairs (currently susceptible AND currently tested) — fixed at t0,
# with their history-based features. Future testing/labels are applied per horizon.
base = []
for p in SP:
    Mrow = Rh0.loc[p]
    for a in AB:
        if bool(Rh.loc[p, a]): continue           # must be susceptible now
        if pd.isna(pr_h.loc[p, a]): continue       # must be tested now
        rho_abx = (Waa[a] * Mrow).sum() / Was[a] if Was[a] > 0 else 0.0
        colv = Rh0[a].reindex(SP).fillna(0)
        rho_path = (Wpp.loc[p] * colv).sum() / Wps[p] if Wps[p] > 0 else 0.0
        base.append((p, a, rho_abx, rho_path, pop[a]))
BASE = pd.DataFrame(base, columns=["Species", "Antibiotic", "rho_abx", "rho_path", "prev"])
print(f"[1/2] trained at t0={T0}: {len(BASE)} currently-susceptible tested candidate pairs\n")

# ---------------------------------------------------------------- evaluate each horizon
def evaluate_horizon(N):
    """Score the fixed candidate set against the future target year t0+N, keeping
    ONLY pairs that are actually tested in that future year (question 2)."""
    ty = T0 + N
    pr_n = win(ty, ty); Rn = rca_bin(pr_n)
    kept, dropped = [], 0
    for r in BASE.itertuples(index=False):
        if r.Species not in pr_n.index or r.Antibiotic not in pr_n.columns or pd.isna(pr_n.loc[r.Species, r.Antibiotic]):
            dropped += 1                            # not tested in the future -> not evaluable
            continue
        kept.append(dict(N=N, target_year=ty, Species=r.Species, Antibiotic=r.Antibiotic,
                         rho_abx=r.rho_abx, rho_path=r.rho_path, prev=r.prev,
                         label=int(bool(Rn.loc[r.Species, r.Antibiotic]))))
    return pd.DataFrame(kept), dropped

print("[1/2] evaluating horizons ...")
all_cand = []; summary = []
for N in sorted(args.horizons):
    if T0 + N > YMAX:
        print(f"      N={N}: target {T0 + N} beyond data ({YMAX}) — skipped"); continue
    C, dropped = evaluate_horizon(N)
    if len(C) == 0:
        print(f"      N={N}: no evaluable pairs — skipped"); continue
    all_cand.append(C)
    kept = len(C); pos = int(C.label.sum())
    row = dict(horizon=N, target_year=T0 + N,
               n_candidates=kept, n_positive=pos, base_rate=round(C.label.mean(), 4),
               n_dropped_future_untested=dropped,
               frac_kept=round(kept / (kept + dropped), 4) if (kept + dropped) else np.nan,
               AUC_rho_abx=auc(C.label, C.rho_abx),
               AUC_rho_path=auc(C.label, C.rho_path),
               AUC_prevalence=auc(C.label, C.prev))
    summary.append(row)
    print(f"      N={N} (->{T0 + N}): {kept:5d} tested pairs kept, {dropped:5d} dropped (untested in future) | "
          f"pos {pos} (base {row['base_rate']:.3f}) | "
          f"AUC abx={row['AUC_rho_abx']:.3f} path={row['AUC_rho_path']:.3f} prev={row['AUC_prevalence']:.3f}")

S = pd.DataFrame(summary)
S.to_csv(f"{DATA}/horizon_prediction_auc.csv", index=False)
pd.concat(all_cand, ignore_index=True).to_csv(f"{DATA}/horizon_prediction_candidates.csv", index=False)
print(f"      wrote horizon_prediction_auc.csv and horizon_prediction_candidates.csv")

# ---------------------------------------------------------------- figure
print("[2/2] rendering fig6_horizon ...")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 12.5,
    'font.family': 'sans-serif', 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
# (a) AUC vs horizon
for col, c, lab in [("AUC_rho_abx", "#c9184a", "antibiotic-net  " + r"$\rho_{abx}$"),
                    ("AUC_rho_path", "#0077b6", "pathogen-net  " + r"$\rho_{path}$"),
                    ("AUC_prevalence", "#adb5bd", "prevalence baseline")]:
    ax[0].plot(S.horizon, S[col], "o-", color=c, lw=2, ms=7, label=lab)
ax[0].axhline(0.5, ls="--", color="k", alpha=0.5, label="random (AUC=0.5)")
ax[0].set_xticks(S.horizon); ax[0].set_ylim(0.45, 1.0)
ax[0].set_xlabel("forecast horizon N (years ahead of %d)" % T0); ax[0].set_ylabel("out-of-sample AUC")
ax[0].set_title(f"(a) Predictability vs horizon\n(train ≤{T0}, tested-in-future pairs only)")
ax[0].legend(frameon=False, fontsize=9)
# (b) candidate count / positives per horizon
w = 0.38; x = np.arange(len(S))
ax[1].bar(x - w / 2, S.n_candidates, w, color="#0077b6", alpha=0.85, label="scored pairs (tested in future)")
ax[1].bar(x + w / 2, S.n_positive, w, color="#c9184a", alpha=0.85, label="positives (acquired R)")
ax[1].set_xticks(x); ax[1].set_xticklabels([f"{n}\n(→{T0 + n})" for n in S.horizon])
ax[1].set_xlabel("forecast horizon N (years ahead)"); ax[1].set_ylabel("count")
ax[1].set_title("(b) Evaluable pairs per horizon")
ax[1].legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(f"{FIG}/fig6_horizon.pdf"); fig.savefig(f"{FIG}/fig6_horizon.png", dpi=200); plt.close(fig)
print(f"      wrote {FIG}/fig6_horizon.png\n")

print("DONE.")
print(S.to_string(index=False))
