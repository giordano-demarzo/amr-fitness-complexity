# The Resistome Atlas — interactive dashboard

A **zero-dependency static website** that lets anyone explore the paper's results:
search a pathogen, watch its resistance ignite across the antibiotic co-resistance
network year by year, read its resistance timeline, and see a forward-looking
forecast of the drugs it is most likely to defeat next.

Everything is precomputed with the paper's own functions, so every number on screen
matches `reproduce_paper.py`.

## Run it locally

The page loads JSON with `fetch()`, so it must be served over HTTP (opening
`index.html` from disk will be blocked by the browser):

```bash
cd dashboard
python3 -m http.server 8000
# then open http://localhost:8000
```

## Rebuild the data

If the input surveillance table or the paper outputs change, regenerate the JSON:

```bash
python3 dashboard/build_dashboard_data.py
```

This reads:
- `results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv` (the input table)
- `results/nbi_projection/*.csv` (the paper's committed outputs)
- `data/reference/*` (WHO AWaRe + ESKAPEE labels)

and writes compact JSON into `dashboard/data/`:

| file | contents |
|------|----------|
| `meta.json` | headline stats, year range, AUCs, nestedness ratio |
| `antibiotics.json` | 43 antibiotic nodes: AWaRe tier, complexity, network layout, predictive AUC |
| `network.json` | pruned co-resistance edges for display |
| `pathogens.json` | 394-pathogen index (ESKAPE status, fitness, counts) |
| `pathogen_detail.json` | per-pathogen resistance timelines, network-ignition years, risk forecast |

## Host it online

It is a plain static folder — no build step, no server code. Any static host works:

- **GitHub Pages** — push the repo and enable Pages, then serve from `/dashboard`
  (or copy the folder to the Pages root).
- **Netlify / Vercel / Cloudflare Pages** — drag-and-drop the `dashboard/` folder,
  or point the project at it with no build command and publish directory `dashboard`.
- **Any bucket** (S3 / GCS / Azure) — upload the folder and enable static hosting.

### Shareable links

Selecting a pathogen updates the URL, so links are shareable and deep-linkable:

```
?p=Klebsiella%20pneumoniae        # open straight to a pathogen
?p=Escherichia%20coli&y=2012      # …at a specific year in the timeline
```

## Files

```
index.html                 layout + copy
style.css                  dark scientific theme (AWaRe palette)
app.js                     rendering + interaction (hand-written SVG, no libraries)
build_dashboard_data.py    data pipeline -> data/*.json
data/                      generated JSON (committed for convenience)
```

> Research prototype for the Vivli AMR Surveillance Open Data Challenge. Forecasts
> are statistical signals from historical co-resistance patterns, **not clinical guidance**.
