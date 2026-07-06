#!/usr/bin/env python3
"""
build_dashboard_data.py — precompute compact JSON for the static Resistome Atlas
dashboard, using the SAME core functions as reproduce_paper.py so every number
shown in the browser matches the manuscript.

Outputs (into dashboard/data/):
    meta.json             headline stats, year range, methodology numbers
    antibiotics.json      43 antibiotic nodes: AWaRe tier, complexity, layout, AUC
    network.json          pruned co-resistance edges for the antibiotic network
    pathogens.json        394-pathogen index (name, ESKAPE, fitness, counts)
    pathogen_detail.json  per-pathogen timelines, network-lighting, risk forecast

Run:  python dashboard/build_dashboard_data.py
"""
import os, json, warnings, math
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INPUT = f"{ROOT}/results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv"
ANA   = f"{ROOT}/results/analysis"
REF   = f"{ROOT}/data/reference"
OUT   = f"{HERE}/data"
os.makedirs(OUT, exist_ok=True)

MIN = 30
POOL = (2015, 2024)
YEARS = list(range(2004, 2025))

# ---- reference labels (same source files as the paper) ----------------------
_aw = pd.read_csv(f"{REF}/atlas_antibiotic_aware_map.csv")
_aw = _aw[_aw.aware_tier.astype(str).str.strip() != ""]
AWARE = {r.atlas_name: int(r.aware_tier) for r in _aw.itertuples()}
AWARE_CAT = {r.atlas_name: str(r.aware_category) for r in _aw.itertuples()}
ESK = list(pd.read_csv(f"{REF}/eskape_pathogens.csv")["match_string"].str.strip())
ESK_LABEL = {r.match_string.strip(): str(r.label) for r in
             pd.read_csv(f"{REF}/eskape_pathogens.csv").itertuples()}

df = pd.read_csv(INPUT)

# ---- core helpers (verbatim maths from reproduce_paper.py) ------------------
def win(y0, y1):
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

def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float); n1 = y.sum(); n0 = len(y) - n1
    if n1 < 3 or n0 < 3: return np.nan
    r = pd.Series(s).rank().values; return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def r2(x):  # round helper -> keep JSON small
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), 4)

def priority(s):
    s = s.lower()
    crit = ["acinetobacter baumannii", "klebsiella pneumoniae", "escherichia coli", "enterobacter"]
    high = ["enterococcus faecium", "enterococcus faecalis", "pseudomonas aeruginosa", "staphylococcus aureus",
            "salmonella", "shigella", "neisseria gonorrhoeae"]
    med = ["streptococcus pneumoniae", "haemophilus influenzae", "streptococcus agalactiae", "streptococcus pyogenes"]
    if any(x in s for x in crit): return 3
    if any(x in s for x in high): return 2
    if any(x in s for x in med): return 1
    return 0

print("[1] building canonical antibiotic node set + co-resistance network ...")
ALL_ABX = sorted(df.Antibiotic.unique())
ab_index = {a: i for i, a in enumerate(ALL_ABX)}
propR = win(*POOL)
RCA = rca_val(propR)
# structural matrix (paper's Fig 2/3 construction: drop empty rows AND columns)
Mfull = (RCA >= 1).astype(float).reindex(columns=ALL_ABX).fillna(0.0)
Mfull = Mfull.loc[Mfull.sum(1) > 0, Mfull.sum(0) > 0]
Mv = Mfull.to_numpy(float)
sp_net = np.array(Mfull.index); ab_net = np.array(Mfull.columns)

# Raw co-occurrence counts drive the FORECAST (see below), but their range (5..500+)
# collapses a force layout into a hub-dominated clump. For *display only* we position
# and connect nodes with the degree-normalised co-occurrence (cosine similarity), which
# spreads the drugs evenly and exposes the co-resistance clusters. The edges shown are
# therefore the strongest co-occurrence relationships, controlled for how common each drug is.
Craw = cooc_cols(Mfull).to_numpy()
kdeg = Mfull.sum(0).reindex(ab_net).to_numpy()            # # species resistant to each drug
den = np.sqrt(np.outer(kdeg, kdeg)); den[den == 0] = 1.0
Wviz = Craw / den; np.fill_diagonal(Wviz, 0.0)
Waa = pd.DataFrame(Wviz, index=ab_net, columns=ab_net)    # cosine-normalised, for layout+edges
Fg, Qg = efc(Mv)
Fglob = pd.Series(Fg, index=sp_net)
Qglob = pd.Series(Qg, index=ab_net)

# spectral blocks (same as paper) for in-block z-scored fitness/complexity
from sklearn.cluster import SpectralCoclustering
cc = SpectralCoclustering(n_clusters=5, random_state=42).fit(Mv)
rb, cb = cc.row_labels_, cc.column_labels_
Qz = pd.Series(index=ab_net, dtype=float); Fz = pd.Series(index=sp_net, dtype=float)
for b in range(5):
    ri = np.where(rb == b)[0]; ci = np.where(cb == b)[0]; sub = Mv[np.ix_(ri, ci)]
    kr = sub.sum(1) > 0; kc = sub.sum(0) > 0; ri, ci = ri[kr], ci[kc]; sub = Mv[np.ix_(ri, ci)]
    if sub.shape[0] < 3 or sub.shape[1] < 3: continue
    F, Q = efc(sub); Q = pd.Series(Q, index=ab_net[ci]); F = pd.Series(F, index=sp_net[ri])
    if Q.std() > 0: Qz[Q.index] = (Q - Q.mean()) / Q.std()
    if F.std() > 0: Fz[F.index] = (F - F.mean()) / F.std()

# ---- deterministic 2-D layout via networkx spring (matches paper's Fig4 look)
import networkx as nx
G = nx.from_pandas_adjacency(Waa)
G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d["weight"] <= 0])
pos = nx.spring_layout(G, weight="weight", seed=1, k=1.7, iterations=500)
deg = dict(G.degree(weight="weight"))
xs = np.array([pos[n][0] for n in ab_net]); ys = np.array([pos[n][1] for n in ab_net])
xs = (xs - xs.min()) / (xs.ptp() or 1); ys = (ys - ys.min()) / (ys.ptp() or 1)

# collision relaxation in the same pixel frame the front-end uses (W=1000,H=700,pad=40,
# radius = 6 + sqrt(deg/max)*12): the spring keeps strongly-linked drugs on top of one
# another, so nudge overlapping disks apart until each has a little breathing room.
Wpx, Hpx, pad = 1000.0, 700.0, 40.0
mx = max(deg.values()) or 1.0
rad = np.array([6.0 + np.sqrt(deg.get(n, 0.0) / mx) * 12.0 for n in ab_net])
px = pad + xs * (Wpx - 2 * pad); py = pad + (1 - ys) * (Hpx - 2 * pad)
for _ in range(250):
    moved = False
    for i in range(len(ab_net)):
        for j in range(i + 1, len(ab_net)):
            dx = px[j] - px[i]; dy = py[j] - py[i]; dist = float(np.hypot(dx, dy))
            if dist < 1e-6: dx, dy, dist = 0.9, 0.3, 0.95
            need = rad[i] + rad[j] + 15.0
            if dist < need:
                sh = (need - dist) / 2.0; ux, uy = dx / dist, dy / dist
                px[i] -= ux * sh; py[i] -= uy * sh; px[j] += ux * sh; py[j] += uy * sh
                moved = True
    if not moved: break
px = np.clip(px, pad, Wpx - pad); py = np.clip(py, pad, Hpx - pad)
xs = (px - pad) / (Wpx - 2 * pad); ys = 1 - (py - pad) / (Hpx - 2 * pad)
LAYOUT = {n: (float(xs[i]), float(ys[i])) for i, n in enumerate(ab_net)}

# ---- per-antibiotic predictive AUC (contagion predictability) --------------
clean = pd.read_csv(f"{ANA}/clean_prediction.csv")
abx_auc = {}
for a, grp in clean.groupby("Antibiotic"):
    v = auc(grp.label, grp.rho_abx)
    if v is not None and np.isfinite(v) and grp.label.sum() >= 5:
        abx_auc[a] = float(v)

deg = dict(G.degree(weight="weight"))
antibiotics = []
for a in ALL_ABX:
    innet = a in LAYOUT
    antibiotics.append(dict(
        name=a,
        tier=AWARE.get(a, 0),
        tier_name=AWARE_CAT.get(a, "Unclassified"),
        x=r2(LAYOUT[a][0]) if innet else None,
        y=r2(LAYOUT[a][1]) if innet else None,
        complexity=r2(Qz.get(a)) if a in Qz.index and np.isfinite(Qz.get(a, np.nan)) else None,
        degree=r2(deg.get(a, 0.0)) if innet else None,
        auc=r2(abx_auc.get(a)),
    ))

# ---- prune edges for a legible network (top-k per node, unioned) -----------
print("[2] pruning network edges for display ...")
edges = []
seen = set()
W = Waa.to_numpy()
K = 4
for i, a in enumerate(ab_net):
    order = np.argsort(-W[i])
    kept = 0
    for j in order:
        if j == i or W[i, j] <= 0: continue
        key = tuple(sorted((ab_index[a], ab_index[ab_net[j]])))
        if key not in seen:
            seen.add(key)
            edges.append([key[0], key[1], r2(W[i, j])])
        kept += 1
        if kept >= K: break
wmax = max((e[2] for e in edges), default=1.0)
for e in edges:
    e[2] = r2(e[2] / wmax)  # normalise weight 0..1 for line opacity/width

# ================================================================ PATHOGENS
print("[3] per-pathogen timelines + cumulative network lighting ...")
# cumulative RCA>=1 per snapshot year (matches Fig 5 "light up and never revert")
first_lit = {}   # species -> {abx_index: first_year}
cum = None
for y in YEARS:
    st = (rca_val(win(2004, y)) >= 1).reindex(columns=ALL_ABX)
    st = st.fillna(False)
    cum = st if cum is None else (cum | st)
    for s in cum.index:
        row = cum.loc[s]
        d = first_lit.setdefault(s, {})
        for a in ALL_ABX:
            if bool(row.get(a, False)) and a not in d:
                d[a] = y

# per (species, antibiotic) real resistance-fraction timeline (from raw rows)
timelines = {}
for r in df.itertuples():
    timelines.setdefault(r.Species, {}).setdefault(r.Antibiotic, []).append(
        [int(r.Year), round(float(r.prop_R), 3), int(r.n_tested), int(r.n_R)])
for s in timelines:
    for a in timelines[s]:
        timelines[s][a].sort()

# ---- forward risk forecast: history through 2024 -> P(resistance next) ------
print("[4] forward resistance-risk forecast (history -> 2025) ...")
propR_all = win(2004, 2024)
Rh = (rca_val(propR_all) >= 1)
Rh0 = Rh.astype(float)
Wf = cooc_cols(Rh); Wf_s = Wf.sum(1)
# calibrate rho_abx -> probability with a monotone logistic fit on the back-test
from sklearn.linear_model import LogisticRegression
cal = LogisticRegression()
cal.fit(clean[["rho_abx"]].values, clean["label"].values)
def prob_from_rho(rho):
    return float(cal.predict_proba([[rho]])[0, 1])
base_rate = float(clean.label.mean())

forecasts = {}
ab_hist = list(Rh.columns)
for p in Rh.index:
    if Wf_s is None: break
    Mrow = Rh0.loc[p]
    tested = set(timelines.get(p, {}).keys())
    cand = []
    for a in ab_hist:
        if a not in tested: continue
        if bool(Rh.loc[p, a]): continue                 # already resistant
        if pd.isna(propR_all.loc[p, a]): continue       # must be observed
        if Wf_s[a] <= 0: continue
        cur = timelines[p][a][-1] if timelines.get(p, {}).get(a) else None
        # only forecast drugs the pathogen can still largely be treated with
        if cur is not None and cur[1] >= 0.5: continue
        rho = float((Wf[a] * Mrow).sum() / Wf_s[a])
        cand.append(dict(abx=a, rho=r2(rho), prob=r2(prob_from_rho(rho)),
                         cur_propR=(cur[1] if cur else None), cur_year=(cur[0] if cur else None)))
    cand.sort(key=lambda d: -(d["rho"] or 0))
    if cand:
        forecasts[p] = cand[:8]

# ---- pathogen index --------------------------------------------------------
n_resist_now = {s: sum(1 for a in first_lit.get(s, {})) for s in first_lit}
fit_pct = Fz.rank(pct=True) * 100  # in-block fitness percentile (the meaningful axis)
pathogens = []
detail = {}
for s in sorted(df.Species.unique()):
    sl = df[df.Species == s]
    esk_hit = [ESK_LABEL[x] for x in ESK if x in s.lower()]
    esk = "ESKAPEE" if any("ESKAPEE" in e for e in esk_hit) else ("ESKAPE" if esk_hit else None)
    n_rec = int(len(sl))
    yrs = sorted(sl.Year.unique())
    tl = timelines.get(s, {})
    lit = first_lit.get(s, {})
    fc = forecasts.get(s, [])
    pathogens.append(dict(
        name=s,
        esk=esk,
        priority=priority(s),
        fitness_z=r2(Fz.get(s)) if s in Fz.index and np.isfinite(Fz.get(s, np.nan)) else None,
        fitness_pct=r2(fit_pct.get(s)) if s in fit_pct.index and np.isfinite(fit_pct.get(s, np.nan)) else None,
        n_abx=len(tl),
        n_resist=n_resist_now.get(s, 0),
        y0=int(min(yrs)) if (yrs := yrs) and yrs else None,
        y1=int(max(yrs)) if yrs else None,
        n_records=n_rec,
        top_risk=(fc[0]["prob"] if fc else None),
    ))
    detail[s] = dict(
        timeline={a: v for a, v in tl.items()},
        lit={ab_index[a]: y for a, y in lit.items()},
        forecast=fc,
    )

# ================================================================ META
print("[5] meta / headline stats ...")
nest = pd.read_csv(f"{ANA}/nestedness_stats.csv").iloc[0].to_dict()
val = pd.read_csv(f"{ANA}/validation_stats.csv").iloc[0].to_dict()
pooled_auc_abx = auc(clean.label, clean.rho_abx)
pooled_auc_prev = auc(clean.label, clean.prev)
meta = dict(
    years=YEARS,
    n_species=int(df.Species.nunique()),
    n_antibiotics=int(df.Antibiotic.nunique()),
    n_records=int(len(df)),
    n_countries=int(df.n_countries.max()) if "n_countries" in df else None,
    year_min=int(df.Year.min()), year_max=int(df.Year.max()),
    nest_ratio=r2(nest["ratio"]), nest_zI=r2(nest["zI"]), nest_zN=r2(nest["zN"]),
    auc_abx=r2(pooled_auc_abx), auc_prev=r2(pooled_auc_prev),
    rho_complexity_inblock=r2(val["rho_complexity_inblock"]),
    p_complexity_inblock=val["p_complexity_inblock"],
    eskapee_fit_med=r2(val["eskapee_fitness_median"]),
    noneskapee_fit_med=r2(val["noneskapee_fitness_median"]),
    p_eskapee=val["p_eskapee_mannwhitney"],
    base_rate=r2(base_rate),
    n_forecast_pathogens=len(forecasts),
)

# ================================================================ WRITE
def dump(name, obj):
    p = f"{OUT}/{name}"
    with open(p, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    print(f"    wrote {name}  ({os.path.getsize(p)/1024:.0f} KB)")

dump("meta.json", meta)
dump("antibiotics.json", antibiotics)
dump("network.json", dict(edges=edges))
dump("pathogens.json", pathogens)
dump("pathogen_detail.json", detail)
print("DONE.")
