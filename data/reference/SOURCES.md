# Data sources & provenance

Every dataset the analysis consumes, where it came from, and how it maps into
the pipeline (`reproduce_paper.py`). Access date for all downloads: 2026-07-04.

---

## 1. Primary dataset — antimicrobial susceptibility counts

**File used by the pipeline:**
`results/temporal_diagnostics/pathogen_year_antibiotic_summary_minN_30.csv`

- **What it is:** species × year × antibiotic counts of susceptible / intermediate /
  resistant isolates (columns: Species, Year, Antibiotic, n_tested, n_R, n_S, n_I,
  prop_R, …). 394 species, 43 antibiotics, 2004–2024; cells with `n_tested ≥ 30`.
- **Source:** the **Pfizer ATLAS** (Antimicrobial Testing Leadership and Surveillance)
  programme, distributed through the **Vivli** data-sharing platform
  (https://vivli.org/) as part of the Vivli AMR Data Challenge.
- **Raw origin:** isolate-level download `atlas_vivli_2004_2024.csv` (not stored in
  this repo copy), aggregated to the species × year × antibiotic summary above by the
  upstream `src/` pipeline.
- **Citation:** Pfizer ATLAS / Vivli AMR surveillance programme (2004–2024).

## 2. Antibiotic stewardship classification — WHO AWaRe 2023

**Files:**
- `data/reference/who_aware_classification_2023.xlsx` — verbatim WHO download.
- `data/reference/who_aware_classification_2023.csv` — the 255 classified
  antibiotics extracted from that file (Antibiotic, Class, ATC, Category, EML).
- `data/reference/atlas_antibiotic_aware_map.csv` — the 43 ATLAS antibiotics mapped
  to their AWaRe tier (Access=1, Watch=2, Reserve=3), with the matched WHO name and
  a note for any name/formulation caveat. **This is the file the pipeline reads.**

- **Source:** World Health Organization, *AWaRe classification of antibiotics for
  evaluation and monitoring of use, 2023*, publication **WHO-MHP-HPS-EML-2023.04**.
  Landing page: https://www.who.int/publications/i/item/WHO-MHP-HPS-EML-2023.04
  File: https://iris.who.int/server/api/core/bitstreams/abba5c2a-8457-431c-9695-16cb2317dd0e/content
- **Used for:** Fig 3a (antibiotic complexity vs AWaRe tier) and Fig 4 node colours.
- **Corrections vs the earlier hand-typed labels (verified against this file):**
  - **Cefoxitin**: was Access → **Watch**.
  - **Aztreonam**: was Watch → **Reserve** (moved to Reserve in the 2023 revision).
  - Newly mapped (previously unclassified): Minocycline → Reserve (WHO lists only
    `Minocycline_IV`), Oxacillin → Access, Tetracycline → Access,
    Quinupristin-dalfopristin → Reserve (WHO name `Dalfopristin/quinupristin`).
  - IV/oral note: for metronidazole, vancomycin, colistin, minocycline the WHO tier
    used is the parenteral (IV) entry; ATLAS MIC testing does not distinguish route.

## 3. Priority-pathogen labels — ESKAPE / ESKAPEE

**File:** `data/reference/eskape_pathogens.csv` (match_string, label, note, source).

- **Source:** the **ESKAPE** acronym — Rice LB, *J Infect Dis* 2008;197:1079
  (Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae,
  Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacter spp.).
- **Extension:** *Escherichia coli* is added, making the set **ESKAPEE** (the common
  7-organism extension). Kept by explicit user decision ("more data is better"); the
  paper therefore refers to this group as ESKAPEE.
- **Used for:** Fig 3b/c (pathogen fitness, ESKAPEE vs non-ESKAPEE).
