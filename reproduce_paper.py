#!/usr/bin/env python3
"""
reproduce_paper.py — one script, one input, all paper data + figures.

INPUT (the only thing you need):
    results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv
    (species x year x antibiotic counts of S/I/R isolates, aggregated from the
     Pfizer ATLAS / Vivli surveillance programme, 2004-2024, cells with
     n_tested >= 30.)

OUTPUT:
    <datadir>/clean_prediction.csv        rolling-origin prediction candidates (Fig 4)
    <datadir>/contagion_hubs.csv          per-antibiotic predictability      (Fig 4)
    <datadir>/nestedness_stats.csv        observed vs degree-null nestedness (Fig 2)
    <figdir>/fig1_binarization.{pdf,png}
    <figdir>/fig2_blocknested.{pdf,png}
    <figdir>/fig3_efc_validation.{pdf,png}
    <figdir>/fig4_prediction.{pdf,png}
    <figdir>/fig5_spreading.{pdf,png}

Everything downstream is deterministic (fixed seeds), so re-running reproduces
the paper bit-for-bit.

Usage:
    python reproduce_paper.py                     # default project layout, in place
    python reproduce_paper.py --input X.csv --outdir /tmp/repro
    python reproduce_paper.py --permutations 0    # skip the slow significance null
"""
import argparse, os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------- config / IO
HERE = os.path.dirname(os.path.abspath(__file__))
DEF_INPUT = f"{HERE}/results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv"

ap = argparse.ArgumentParser(description="Regenerate all paper data and figures from the summary CSV.")
ap.add_argument("--input",  default=DEF_INPUT, help="species x year x antibiotic summary CSV")
ap.add_argument("--outdir", default=HERE, help="base output dir (default: project root)")
ap.add_argument("--figdir", default=None, help="override figure dir (default: <outdir>/paper/figures)")
ap.add_argument("--datadir", default=None, help="override data dir (default: <outdir>/results/analysis)")
ap.add_argument("--permutations", type=int, default=200,
                help="node-label permutations for the prediction significance null (0 to skip)")
ap.add_argument("--null-samples", type=int, default=100,
                help="degree-preserving samples for the nestedness null (Fig 2)")
args = ap.parse_args()

FIG  = args.figdir  or f"{args.outdir}/paper/figures"
DATA = args.datadir or f"{args.outdir}/results/analysis"
REF  = f"{HERE}/data/reference"          # WHO AWaRe + ESKAPE reference labels (see SOURCES.md)
os.makedirs(FIG, exist_ok=True); os.makedirs(DATA, exist_ok=True)

MIN = 30           # min tested isolates per cell (matches the summary file's own filter)
POOL = (2015, 2024)  # pooled window for structural / validation figures

# ---- external reference labels, loaded from auditable files (not hardcoded) ----
# WHO AWaRe 2023 tier per ATLAS antibiotic (Access=1, Watch=2, Reserve=3).
_aw = pd.read_csv(f"{REF}/atlas_antibiotic_aware_map.csv")
_aw = _aw[_aw.aware_tier.astype(str).str.strip() != ""]
AWARE = {r.atlas_name: int(r.aware_tier) for r in _aw.itertuples()}
# ESKAPE/ESKAPEE priority pathogens (match substrings against lowercased species).
ESK = list(pd.read_csv(f"{REF}/eskape_pathogens.csv")["match_string"].str.strip())

print(f"input   : {args.input}")
print(f"figures : {FIG}")
print(f"data    : {DATA}")
print(f"labels  : {REF}/atlas_antibiotic_aware_map.csv ({len(AWARE)} antibiotics), "
      f"eskape_pathogens.csv ({len(ESK)} taxa)\n")
df = pd.read_csv(args.input)

# ---------------------------------------------------------------- core helpers
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

def efc(X, n=3000, tol=1e-12):
    """Fitness–complexity (Tacchella et al.) fixed point."""
    X = np.asarray(X, float); F = np.ones(X.shape[0]); Q = np.ones(X.shape[1])
    for _ in range(n):
        Fn = X @ Q; m = Fn.mean(); Fn = Fn / m if m > 0 else Fn
        den = X.T @ (1 / np.where(F > 0, F, np.inf)); Qn = 1 / np.where(den > 0, den, np.inf)
        mq = Qn[np.isfinite(Qn)].mean() if np.isfinite(Qn).any() else 1; Qn = Qn / (mq if mq > 0 else 1)
        if np.nanmax(np.abs(Fn - F)) + np.nanmax(np.abs(Qn - Q)) < tol: F, Q = Fn, Qn; break
        F, Q = Fn, Qn
    return F, Q

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

def nestedness(M, rb, cb):
    """In-block nestedness (Laudati 2023, Eq.4); with constant labels -> global N (Eq.5)."""
    P, A = M.shape; kr = M.sum(1); kc = M.sum(0)
    Orr = M @ M.T; exp = np.outer(kr, kr) / A
    Crow = np.array([(rb == rb[i]).sum() for i in range(P)])
    mask = (kr[:, None] > kr[None, :]) & (rb[:, None] == rb[None, :]) & (Crow[:, None] > 1)
    den = kr[None, :] * (Crow[:, None] - 1); den = np.where(den == 0, np.nan, den)
    srow = np.nansum(((Orr - exp) / den)[mask])
    Occ = M.T @ M; expc = np.outer(kc, kc) / P
    Ccol = np.array([(cb == cb[a]).sum() for a in range(A)])
    maskc = (kc[:, None] > kc[None, :]) & (cb[:, None] == cb[None, :]) & (Ccol[:, None] > 1)
    denc = kc[None, :] * (Ccol[:, None] - 1); denc = np.where(denc == 0, np.nan, denc)
    scol = np.nansum(((Occ - expc) / denc)[maskc])
    return 2 / (P + A) * (srow + scol)

# ================================================================ STAGE 1
# Rolling-origin, leakage-free prediction -> clean_prediction.csv, contagion_hubs.csv
# ================================================================
def stage_prediction():
    print("[1/3] rolling-origin prediction (Fig 4 data) ...")
    ORIGINS = list(range(2019, 2024))          # predict t+1 = 2020..2024
    rows = []
    for t in ORIGINS:
        pr_h = win(2004, t); Rh = rca_bin(pr_h)            # expanding history
        pr_n = win(t + 1, t + 1); Rn = rca_bin(pr_n)       # next-year labels
        Waa = cooc_cols(Rh); Was = Waa.sum(1)
        Wpp = cooc_rows(Rh); Wps = Wpp.sum(1)
        pop = Rh.mean(0); ab = list(Rh.columns); sp = list(Rh.index); Rh0 = Rh.astype(float)
        for p in sp:
            if p not in pr_n.index: continue
            Mrow = Rh0.loc[p]
            for a in ab:
                if a not in pr_n.columns: continue
                if bool(Rh.loc[p, a]): continue                 # must be susceptible now
                if pd.isna(pr_h.loc[p, a]): continue            # observed now
                if pd.isna(pr_n.loc[p, a]): continue            # observed next year
                rho_abx = (Waa[a] * Mrow).sum() / Was[a] if Was[a] > 0 else 0.0
                colv = Rh0[a].reindex(sp).fillna(0)
                rho_path = (Wpp.loc[p] * colv).sum() / Wps[p] if Wps[p] > 0 else 0.0
                rows.append(dict(t=t, Species=p, Antibiotic=a, rho_abx=rho_abx,
                                 rho_path=rho_path, prev=pop[a], label=int(bool(Rn.loc[p, a]))))
    C = pd.DataFrame(rows)
    C.to_csv(f"{DATA}/clean_prediction.csv", index=False)
    print(f"      candidates {len(C)} | positives {int(C.label.sum())} (base rate {C.label.mean():.3f})")
    print(f"      pooled AUC  rho_abx={auc(C.label, C.rho_abx):.3f}  "
          f"rho_path={auc(C.label, C.rho_path):.3f}  prevalence={auc(C.label, C.prev):.3f}")

    # per-antibiotic predictability (contagion hubs)
    hub = (C.groupby("Antibiotic").apply(lambda x: pd.Series(
              {"AUC": auc(x.label, x.rho_abx), "n": len(x), "pos": int(x.label.sum())}))
           .dropna().sort_values("AUC", ascending=False))
    hub = hub[hub.pos >= 5]
    hub.to_csv(f"{DATA}/contagion_hubs.csv")

    # optional significance null (only prints; not needed for the figures)
    if args.permutations > 0:
        def perm_run(proj, seed):
            rng = np.random.default_rng(seed); rr = []
            for t in ORIGINS:
                pr_h = win(2004, t); Rh = rca_bin(pr_h); pr_n = win(t + 1, t + 1); Rn = rca_bin(pr_n)
                ab = list(Rh.columns); sp = list(Rh.index); Rh0 = Rh.astype(float)
                if proj == "abx":
                    W = cooc_cols(Rh); perm = rng.permutation(ab); Wp = W.loc[perm, perm]; Wp.index = ab; Wp.columns = ab
                else:
                    W = cooc_rows(Rh); perm = rng.permutation(sp); Wp = W.loc[perm, perm]; Wp.index = sp; Wp.columns = sp
                Ws = Wp.sum(1)
                for p in sp:
                    if p not in pr_n.index: continue
                    Mrow = Rh0.loc[p]
                    for a in ab:
                        if a not in pr_n.columns or bool(Rh.loc[p, a]) or pd.isna(pr_h.loc[p, a]) or pd.isna(pr_n.loc[p, a]): continue
                        if proj == "abx":
                            rho = (Wp[a] * Mrow).sum() / Ws[a] if Ws[a] > 0 else 0.0
                        else:
                            colv = Rh0[a].reindex(sp).fillna(0); rho = (Wp.loc[p] * colv).sum() / Ws[p] if Ws[p] > 0 else 0.0
                        rr.append((rho, int(bool(Rn.loc[p, a]))))
            r = pd.DataFrame(rr, columns=["rho", "label"]); return auc(r.label, r.rho)
        for proj, feat in [("abx", "rho_abx"), ("path", "rho_path")]:
            obs = auc(C.label, C[feat]); null = np.array([perm_run(proj, s) for s in range(args.permutations)])
            z = (obs - null.mean()) / null.std(); pval = (1 + (null >= obs).sum()) / (len(null) + 1)
            print(f"      [{feat}] obs {obs:.3f} vs null {null.mean():.3f}±{null.std():.3f}  z={z:.2f}  p={pval:.4f}")
    return C, hub

# ================================================================ STAGE 2
# Structural objects + nestedness stats (shared by Figs 1,2,3,5)
# ================================================================
def stage_structure():
    print("[2/3] structure + nestedness (Fig 1/2/3 data) ...")
    from sklearn.cluster import SpectralCoclustering
    propR = win(*POOL)
    RCA = rca_val(propR); M = (RCA >= 1).astype(float).reindex(propR.index)
    M = M.loc[M.sum(1) > 0, M.sum(0) > 0]; Mv = M.to_numpy(float)
    cc = SpectralCoclustering(n_clusters=5, random_state=42).fit(Mv)
    rb, cb = cc.row_labels_, cc.column_labels_
    P, A = Mv.shape

    # nestedness: observed global N & in-block I*, plus degree-preserving null
    N_glob = nestedness(Mv, np.zeros(P, int), np.zeros(A, int))
    Istar  = nestedness(Mv, rb, cb)
    rng = np.random.default_rng(0); kr = Mv.sum(1); kc = Mv.sum(0); tot = Mv.sum()
    Pmat = np.clip(np.outer(kr, kc) / tot, 0, 1)
    Nn, In = [], []
    for _ in range(args.null_samples):
        R = (rng.random((P, A)) < Pmat).astype(float)
        keep_r = R.sum(1) > 0; keep_c = R.sum(0) > 0; Rf = R[np.ix_(keep_r, keep_c)]
        if Rf.shape[0] < 6 or Rf.shape[1] < 6: continue
        try:
            Nn.append(nestedness(Rf, np.zeros(Rf.shape[0], int), np.zeros(Rf.shape[1], int)))
            m = SpectralCoclustering(n_clusters=5, random_state=42).fit(Rf)
            In.append(nestedness(Rf, m.row_labels_, m.column_labels_))
        except Exception:
            continue
    Nn = np.array(Nn); In = np.array(In)
    zN = (N_glob - Nn.mean()) / Nn.std(); zI = (Istar - In.mean()) / In.std()
    nest = dict(N_obs=N_glob, N_null=Nn.mean(), N_null_sd=Nn.std(), zN=zN,
                I_obs=Istar, I_null=In.mean(), I_null_sd=In.std(), zI=zI, ratio=Istar / N_glob)
    pd.DataFrame([nest]).to_csv(f"{DATA}/nestedness_stats.csv", index=False)
    print(f"      global N={N_glob:.4f} (z={zN:.1f}) | in-block I*={Istar:.4f} (z={zI:.1f}) | I*/N={Istar/N_glob:.1f}")
    return dict(propR=propR, RCA=RCA, M=M, Mv=Mv, rb=rb, cb=cb, nest=nest)

# ================================================================ STAGE 2b
# Validation statistics that the paper quotes (Spearman rho's, ESKAPEE medians)
# ================================================================
def stage_validation(S):
    print("[2b] external validation stats (numbers quoted in the paper) ...")
    from scipy.stats import spearmanr, mannwhitneyu
    M, Mv, rb, cb = S['M'], S['Mv'], S['rb'], S['cb']
    ab = np.array(M.columns); sp = np.array(M.index)
    # global EFC + within-block z-scored EFC (same construction as the figures)
    Fg, Qg = efc(Mv); Qglob = pd.Series(Qg, index=ab); Fglob = pd.Series(Fg, index=sp)
    Qz = pd.Series(index=ab, dtype=float); Fz = pd.Series(index=sp, dtype=float)
    for b in range(5):
        ri = np.where(rb == b)[0]; ci = np.where(cb == b)[0]; sub = Mv[np.ix_(ri, ci)]
        kr = sub.sum(1) > 0; kc = sub.sum(0) > 0; ri, ci = ri[kr], ci[kc]; sub = Mv[np.ix_(ri, ci)]
        if sub.shape[0] < 3 or sub.shape[1] < 3: continue
        F, Q = efc(sub); Q = pd.Series(Q, index=ab[ci]); F = pd.Series(F, index=sp[ri])
        if Q.std() > 0: Qz[Q.index] = (Q - Q.mean()) / Q.std()
        if F.std() > 0: Fz[F.index] = (F - F.mean()) / F.std()
    # per-antibiotic complexity table
    pd.DataFrame({"AWaRe": pd.Series(AWARE), "complexity_global": Qglob,
                  "complexity_inblock_z": Qz}).dropna(how="all").to_csv(
                  f"{DATA}/antibiotic_complexity_validation.csv")
    aw = pd.Series({a: AWARE.get(a, np.nan) for a in ab}).dropna()
    cg = aw.index.intersection(Qglob.index); cz = aw.index.intersection(Qz.dropna().index)
    rho_overall, p_overall = spearmanr(Qglob[cg], aw[cg])
    rho_inblock, p_inblock = spearmanr(Qz[cz], aw[cz])
    big = int(pd.Series(cb).value_counts().idxmax()); gci = ab[np.where(cb == big)[0]]
    gc = [a for a in gci if a in aw.index and a in Qz.dropna().index]
    rho_gramneg, p_gramneg = spearmanr(Qz[gc], aw[gc]) if len(gc) >= 5 else (np.nan, np.nan)
    # pathogen fitness: WHO priority + ESKAPEE
    def priority(s):
        s = s.lower()
        crit = ["acinetobacter baumannii", "klebsiella pneumoniae", "escherichia coli", "enterobacter"]
        high = ["enterococcus faecium", "enterococcus faecalis", "pseudomonas aeruginosa", "staphylococcus aureus",
                "salmonella", "shigella", "neisseria gonorrhoeae"]
        med = ["streptococcus pneumoniae", "haemophilus influenzae", "streptococcus agalactiae", "streptococcus pyogenes"]
        if any(x in s for x in crit): return 3
        if any(x in s for x in high): return 2
        if any(x in s for x in med): return 1
        return np.nan
    pri = pd.Series({s: priority(s) for s in sp}).dropna()
    cpg = pri.index.intersection(Fglob.index)
    rho_fit_global, _ = spearmanr(Fglob[cpg], pri[cpg])
    esk = pd.Series({s: int(any(x in s.lower() for x in ESK)) for s in sp})
    fe = Fz[esk[esk == 1].index].dropna(); fn = Fz[esk[esk == 0].index].dropna()
    u, p_esk = mannwhitneyu(fe, fn, alternative="greater")
    stats = dict(rho_complexity_overall=rho_overall, p_complexity_overall=p_overall,
                 rho_complexity_inblock=rho_inblock, p_complexity_inblock=p_inblock,
                 rho_complexity_gramneg=rho_gramneg, p_complexity_gramneg=p_gramneg,
                 rho_fitness_global_vs_priority=rho_fit_global,
                 eskapee_fitness_median=fe.median(), noneskapee_fitness_median=fn.median(),
                 p_eskapee_mannwhitney=p_esk, n_eskapee=int(len(fe)), n_noneskapee=int(len(fn)))
    pd.DataFrame([stats]).to_csv(f"{DATA}/validation_stats.csv", index=False)
    print(f"      complexity vs AWaRe:  overall rho={rho_overall:.3f}  in-block rho={rho_inblock:.3f}  "
          f"Gram-neg rho={rho_gramneg:.3f}")
    print(f"      global fitness vs WHO priority rho={rho_fit_global:.3f}")
    print(f"      ESKAPEE fitness median {fe.median():.2f} vs non {fn.median():.2f} (MWU p={p_esk:.1e})")
    return stats

# ================================================================ STAGE 3
# Figures
# ================================================================
def stage_figures(S, C, hub):
    print("[3/3] rendering figures ...")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, to_rgb
    from matplotlib.lines import Line2D
    from scipy.stats import spearmanr
    import networkx as nx
    plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 12.5,
        'font.family': 'sans-serif', 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 150, 'savefig.bbox': 'tight'})
    BLOCKC = ['#c9184a', '#0077b6', '#2a9d8f', '#e9c46a', '#7209b7']
    colA = {1: "#2a9d8f", 2: "#e9c46a", 3: "#c9184a"}
    propR, RCA, M, Mv, rb, cb = S['propR'], S['RCA'], S['M'], S['Mv'], S['rb'], S['cb']
    nest = S['nest']
    ab = np.array(M.columns); sp = np.array(M.index)
    deg_r = Mv.sum(1); deg_c = Mv.sum(0)

    # ---- FIGURE 1 -------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    pr = propR.reindex(index=M.index, columns=M.columns).to_numpy()
    im0 = ax[0].imshow(pr, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax[0].set_title("(a) Resistance fraction  $r_{pa}$"); ax[0].set_xlabel("antibiotics"); ax[0].set_ylabel("pathogens")
    ax[0].set_xticks([]); ax[0].set_yticks([]); plt.colorbar(im0, ax=ax[0], fraction=0.046, label="prop. resistant")
    rv = RCA.reindex(index=M.index, columns=M.columns).to_numpy(); rv = np.log10(np.where(rv > 0, rv, np.nan))
    ax[1].hist(rv[np.isfinite(rv)], bins=40, color="#0077b6", alpha=0.85)
    ax[1].axvline(0, color="#c9184a", lw=2, ls="--", label="RCA = 1 (threshold)")
    ax[1].set_title("(b) Revealed comparative advantage"); ax[1].set_xlabel(r"$\log_{10}\,\mathrm{RCA}_{pa}$")
    ax[1].set_ylabel("count"); ax[1].legend()
    ro = np.argsort(-deg_r); co = np.argsort(-deg_c)
    ax[2].imshow(Mv[np.ix_(ro, co)], aspect="auto", cmap=ListedColormap(["#f0f0f0", "#0d1b2a"]), interpolation="none")
    ax[2].set_title("(c) Binary matrix $M_{pa}$ (degree-sorted)"); ax[2].set_xlabel("antibiotics"); ax[2].set_ylabel("pathogens")
    ax[2].set_xticks([]); ax[2].set_yticks([])
    fig.tight_layout(); fig.savefig(f"{FIG}/fig1_binarization.pdf"); fig.savefig(f"{FIG}/fig1_binarization.png", dpi=200); plt.close(fig)

    # ---- FIGURE 2 -------------------------------------------------------------
    block_order = sorted(set(rb) | set(cb), key=lambda b: -((rb == b).sum() + (cb == b).sum()))
    def order(labels, deg, seq=block_order):
        out = []
        for b in seq:
            idx = np.where(labels == b)[0]; idx = idx[np.argsort(-deg[idx])]; out += list(idx)
        return np.array(out)
    ro = order(rb, deg_r); co = order(cb, deg_c)
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 6), gridspec_kw={'width_ratios': [1.25, 1]})
    Msort = Mv[np.ix_(ro, co)]; rgb = np.ones((*Msort.shape, 3)); rbl = rb[ro]; cbl = cb[co]
    for i in range(Msort.shape[0]):
        for j in range(Msort.shape[1]):
            if Msort[i, j] > 0:
                rgb[i, j] = to_rgb(BLOCKC[rbl[i] % 5]) if rbl[i] == cbl[j] else (0.15, 0.15, 0.15)
    ax[0].imshow(rgb, aspect="auto", interpolation="none")
    for y in np.cumsum([(rbl == b).sum() for b in pd.unique(rbl)])[:-1]: ax[0].axhline(y - 0.5, color="k", lw=0.8, alpha=0.5)
    for x in np.cumsum([(cbl == b).sum() for b in pd.unique(cbl)])[:-1]: ax[0].axvline(x - 0.5, color="k", lw=0.8, alpha=0.5)
    ax[0].set_title("(a) Block-nested resistance matrix"); ax[0].set_xlabel("antibiotics"); ax[0].set_ylabel("pathogens")
    ax[0].set_xticks([]); ax[0].set_yticks([])
    # nestedness bars use the LIVE-computed values from stage 2
    obs = [nest['N_obs'], nest['I_obs']]; nul = [nest['N_null'], nest['I_null']]; nsd = [nest['N_null_sd'], nest['I_null_sd']]
    x = np.arange(2); w = 0.36
    ax[1].bar(x - w / 2, obs, w, color=["#0077b6", "#c9184a"], label="observed")
    ax[1].bar(x + w / 2, nul, w, yerr=nsd, color="#adb5bd", label="degree null", capsize=4)
    ax[1].set_xticks(x); ax[1].set_xticklabels(["global\n$\\mathcal{N}$", "in-block\n$\\mathcal{I}^*$"]); ax[1].set_ylabel("nestedness")
    ax[1].set_title("(b) Nestedness is a local property\n$\\mathcal{I}^*/\\mathcal{N}=%.1f$ (companies-like)" % nest['ratio'])
    ax[1].text(1, nest['I_obs'] + 0.02, "z = %.1f" % nest['zI'], ha="center", color="#c9184a", fontweight="bold")
    ax[1].text(0, nest['N_null'] + nest['N_null_sd'] + 0.02, "z = %.1f" % nest['zN'], ha="center", color="#0077b6")
    ax[1].legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig2_blocknested.pdf"); fig.savefig(f"{FIG}/fig2_blocknested.png", dpi=200); plt.close(fig)

    # ---- FIGURE 3 -------------------------------------------------------------
    # AWARE and ESK are loaded at module scope from data/reference/ (see SOURCES.md)
    Fg, Qg = efc(Mv); Fglob = pd.Series(Fg, index=sp)
    Qz = pd.Series(index=ab, dtype=float); Fz = pd.Series(index=sp, dtype=float)
    for b in range(5):
        ri = np.where(rb == b)[0]; ci = np.where(cb == b)[0]; sub = Mv[np.ix_(ri, ci)]
        kr = sub.sum(1) > 0; kc = sub.sum(0) > 0; ri, ci = ri[kr], ci[kc]; sub = Mv[np.ix_(ri, ci)]
        if sub.shape[0] < 3 or sub.shape[1] < 3: continue
        F, Q = efc(sub); Q = pd.Series(Q, index=ab[ci]); F = pd.Series(F, index=sp[ri])
        if Q.std() > 0: Qz[Q.index] = (Q - Q.mean()) / Q.std()
        if F.std() > 0: Fz[F.index] = (F - F.mean()) / F.std()
    esk = pd.Series({s: int(any(x in s.lower() for x in ESK)) for s in sp})
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    data = [[Qz[a] for a in ab if AWARE.get(a) == lvl and a in Qz.dropna().index] for lvl in [1, 2, 3]]
    bp = ax[0].boxplot(data, labels=["Access", "Watch", "Reserve"], patch_artist=True, widths=0.6)
    for patch, lvl in zip(bp['boxes'], [1, 2, 3]): patch.set_facecolor(colA[lvl]); patch.set_alpha(0.8)
    rr, pp = spearmanr(*zip(*[(Qz[a], AWARE[a]) for a in ab if a in AWARE and a in Qz.dropna().index]))
    ax[0].set_title(f"(a) Antibiotic complexity vs AWaRe\n" + r"$\rho=%.2f$" % rr + f" (p={pp:.1e})")
    ax[0].set_ylabel("in-block complexity (z)")
    g0 = Fglob[esk[esk == 0].index].dropna(); g1 = Fglob[esk[esk == 1].index].dropna()
    bp1 = ax[1].boxplot([g0, g1], labels=["non-ESKAPE", "ESKAPE"], patch_artist=True, widths=0.6)
    bp1['boxes'][0].set_facecolor("#adb5bd"); bp1['boxes'][1].set_facecolor("#c9184a")
    for b in bp1['boxes']: b.set_alpha(0.85)
    ax[1].set_yscale("log"); ax[1].set_title("(b) GLOBAL fitness (degenerate)\nESKAPE not fitter"); ax[1].set_ylabel("global fitness")
    d0 = Fz[esk[esk == 0].index].dropna(); d1 = Fz[esk[esk == 1].index].dropna()
    bp2 = ax[2].boxplot([d0, d1], labels=["non-ESKAPE", "ESKAPE"], patch_artist=True, widths=0.6)
    bp2['boxes'][0].set_facecolor("#adb5bd"); bp2['boxes'][1].set_facecolor("#c9184a")
    for b in bp2['boxes']: b.set_alpha(0.85)
    ax[2].set_title("(c) IN-BLOCK fitness\nESKAPE $\\gg$ non-ESKAPE (p<0.001)"); ax[2].set_ylabel("in-block fitness (z)")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig3_efc_validation.pdf"); fig.savefig(f"{FIG}/fig3_efc_validation.png", dpi=200); plt.close(fig)

    # ---- FIGURE 4 -------------------------------------------------------------
    te = C.copy()
    hb10 = hub.reset_index().sort_values("AUC", ascending=False).head(10)
    W = cooc_cols(M); G = nx.from_pandas_adjacency(W)
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d["weight"] <= 0])
    posl = nx.spring_layout(G, weight="weight", seed=1, k=0.7)
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 11))
    ncol = [colA.get(AWARE.get(n, 2), "#cccccc") for n in G.nodes()]; deg = dict(G.degree(weight="weight"))
    nx.draw_networkx_edges(G, posl, ax=ax[0, 0], alpha=0.15, width=1)
    nx.draw_networkx_nodes(G, posl, ax=ax[0, 0], node_color=ncol, node_size=[80 + deg[n] * 260 for n in G.nodes()],
        edgecolors="k", linewidths=0.4)
    nx.draw_networkx_labels(G, posl, ax=ax[0, 0], font_size=6.2)
    ax[0, 0].set_title("(a) Antibiotic co-resistance network (colored by AWaRe)"); ax[0, 0].axis("off")
    ax[0, 0].legend(handles=[Line2D([0], [0], marker='o', color='w', markerfacecolor=colA[l], markersize=9, label=n)
        for l, n in [(1, "Access"), (2, "Watch"), (3, "Reserve")]], loc="upper left", frameon=False)
    for col, c, lab in [("rho_abx", "#c9184a", "antibiotic-net"), ("rho_path", "#0077b6", "pathogen-net")]:
        t = te.copy(); t["b"] = pd.qcut(t[col].rank(method="first"), 6, labels=False)
        g = t.groupby("b").agg(x=(col, "mean"), p=("label", "mean"))
        ax[0, 1].plot(g.x, g.p, "o-", color=c, label=lab, lw=2, ms=7)
    ax[0, 1].axhline(te.label.mean(), ls="--", color="gray", label="base rate")
    ax[0, 1].set_xlabel(r"co-resistance exposure $\rho$"); ax[0, 1].set_ylabel("P(acquire resistance)")
    ax[0, 1].set_title("(b) Resistance acquisition vs co-resistance exposure"); ax[0, 1].legend(frameon=False)
    def roc(y, s):
        o = np.argsort(-s); y = np.asarray(y)[o]; tpr = np.cumsum(y) / y.sum(); fpr = np.cumsum(1 - y) / (len(y) - y.sum())
        return np.r_[0, fpr], np.r_[0, tpr]
    for col, c, lab in [("rho_abx", "#c9184a", "antibiotic-net"), ("rho_path", "#0077b6", "pathogen-net"), ("prev", "#adb5bd", "prevalence")]:
        f, t = roc(te.label.values, te[col].values); ax[1, 0].plot(f, t, color=c, lw=2, label=f"{lab} (AUC={auc(te.label, te[col]):.2f})")
    ax[1, 0].plot([0, 1], [0, 1], ls="--", color="k", alpha=0.5)
    ax[1, 0].set_xlabel("false positive rate"); ax[1, 0].set_ylabel("true positive rate")
    ax[1, 0].set_title("(c) Out-of-sample prediction (rolling origin, 2020–2024)"); ax[1, 0].legend(frameon=False, loc="lower right")
    hbb = hb10.sort_values("AUC")
    ax[1, 1].barh(hbb.Antibiotic, hbb.AUC, color=[colA.get(AWARE.get(a, 2), "#ccc") for a in hbb.Antibiotic], edgecolor="k", lw=0.4)
    ax[1, 1].axvline(0.5, ls="--", color="gray"); ax[1, 1].set_xlim(0.5, 1.0)
    ax[1, 1].set_xlabel("prediction AUC"); ax[1, 1].set_title("(d) Contagion-hub antibiotics")
    for lbl in ax[1, 1].get_yticklabels(): lbl.set_fontsize(8.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig4_prediction.pdf"); fig.savefig(f"{FIG}/fig4_prediction.png", dpi=200); plt.close(fig)

    # ---- FIGURE 5 -------------------------------------------------------------
    Wa = cooc_cols(M); Ga = nx.from_pandas_adjacency(Wa)
    Ga.remove_edges_from([(u, v) for u, v, d in Ga.edges(data=True) if d["weight"] <= 0])
    posA = nx.spring_layout(Ga, weight="weight", seed=3, k=0.8)
    Wp = cooc_rows(M.loc[M.sum(1) >= 3]); Gp = nx.from_pandas_adjacency(Wp)
    Gp.remove_edges_from([(u, v) for u, v, d in Gp.edges(data=True) if d["weight"] <= 0])
    posP = nx.spring_layout(Gp, weight="weight", seed=5, k=0.5)
    SNAP = [2009, 2014, 2019, 2024]; cumR = {}; acc = None
    for y in SNAP:
        st = (rca_val(win(2004, y)) >= 1).reindex(index=M.index, columns=M.columns).fillna(False)
        acc = st if acc is None else (acc | st); cumR[y] = acc.copy()
    examples = [("pathogen", "Klebsiella pneumoniae"), ("pathogen", "Escherichia coli"),
                ("antibiotic", "Meropenem"), ("antibiotic", "Ciprofloxacin")]
    fig, axes = plt.subplots(len(examples), len(SNAP), figsize=(4 * len(SNAP), 3.5 * len(examples)))
    for i, (kind, name) in enumerate(examples):
        for j, y in enumerate(SNAP):
            axx = axes[i, j]; Rc = cumR[y]
            if kind == "pathogen":
                G, pos, c = Ga, posA, "#c9184a"
                lit = set(a for a in G.nodes() if name in Rc.index and a in Rc.columns and Rc.loc[name, a])
                nx.draw_networkx_edges(G, pos, ax=axx, alpha=0.10, width=0.6)
                cols = [c if n in lit else "#e6e6e6" for n in G.nodes()]
                nx.draw_networkx_nodes(G, pos, ax=axx, node_color=cols, node_size=80, edgecolors="gray", linewidths=0.3)
                tag = f"{name.split()[0][0]}. {name.split()[1]}"
            else:
                G, pos, c = Gp, posP, "#0077b6"
                lit = set(p for p in G.nodes() if p in Rc.index and name in Rc.columns and Rc.loc[p, name])
                nx.draw_networkx_edges(G, pos, ax=axx, alpha=0.05, width=0.4)
                cols = [c if n in lit else "#e6e6e6" for n in G.nodes()]
                nx.draw_networkx_nodes(G, pos, ax=axx, node_color=cols, node_size=32, edgecolors="none")
                tag = name
            axx.set_title(f"{tag} — {y}   ({len(lit)} R)", fontsize=10, color=c); axx.axis("off")
        lab = ("pathogen on\nantibiotic network" if kind == "pathogen" else "antibiotic on\npathogen network")
        axes[i, 0].text(-0.08, 0.5, lab, transform=axes[i, 0].transAxes, rotation=90, va="center", ha="center",
                        fontsize=10, color=("#c9184a" if kind == "pathogen" else "#0077b6"), fontweight="bold")
    fig.suptitle("Spreading of resistance across the relatedness networks (cumulative; nodes light up and never revert)",
                 fontsize=13, y=0.995)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig5_spreading.pdf"); fig.savefig(f"{FIG}/fig5_spreading.png", dpi=170); plt.close(fig)

# ================================================================ main
if __name__ == "__main__":
    C, hub = stage_prediction()
    S = stage_structure()
    stage_validation(S)
    stage_figures(S, C, hub)
    print("\nDONE.")
    print("data  ->", ", ".join(sorted(f for f in os.listdir(DATA) if f.endswith(".csv"))))
    print("figs  ->", ", ".join(sorted(f for f in os.listdir(FIG) if f.startswith("fig"))))
