# Reproducing this analysis in Phoenix WinNonlin

This is the manual half of the bridge. Run it once, export the result, and the
R side picks it up automatically.

## Why manual

Phoenix's documented command-line mode drives the **NLME engine**, not the NCA
object. There is no supported `Phoenix.exe --run-nca` switch. Scripted NCA in
Phoenix goes through one of:

- **AutoPilot Toolkit** — separately licensed module built for exactly this
- **R Script object inside a Phoenix workflow** — R runs *within* Phoenix, so
  Phoenix stays the orchestrator
- **Phoenix WebServices / PKS** — enterprise deployments

For a portfolio repo, export/import is the honest and reproducible pattern:
it works with a plain Phoenix WinNonlin license, produces a checked-in artifact,
and does not depend on a module the reader may not have.

If your site has AutoPilot, the same settings apply; you would automate steps
3–7 instead of clicking them.

---

## 1. Import the data

File → Import → `phoenix/input/theoph_for_phoenix.csv`

Columns: `Subject`, `Time`, `Conc`, `Dose`, `Weight`

## 2. Create the NCA object

Right-click the imported worksheet → Send To → NCA

## 3. Map columns

| Phoenix context | Column   |
|-----------------|----------|
| Sort            | `Subject`|
| Time            | `Time`   |
| Concentration   | `Conc`   |

## 4. Dosing

Dosing tab → **Constant** → `Dose` column
(Or type 320 if using the flat-dose convention directly.)

Model type: **Plasma (200-202)**
Dose type: **Extravascular**

## 5. Calculation settings — these must match `R/00_config.R`

| `NCA_RULES` field   | Value             | Phoenix control                                        |
|---------------------|-------------------|--------------------------------------------------------|
| `auc_method`        | `lin up/log down` | Calculation Method → **Linear Log Trapezoidal**        |
| `min_hl_points`     | `3`               | Lambda Z → minimum 3 points                            |
| `allow_tmax_in_hl`  | `FALSE`           | Lambda Z → exclude Cmax from the regression            |
| `adj_r2_factor`     | `0.0001`          | Lambda Z → **Best Fit** (uses adjusted R² with 0.0001) |
| `max_aucinf_pext`   | `20`              | Acceptance criteria → AUC %extrapolated ≤ 20           |

> Swap to **Linear Trapezoidal / Linear Interpolation** if you set
> `auc_method = "linear"` and `noncompart_down = "Linear"` in the config.
> The two must move together — `assert_rules_consistent()` enforces this on the
> R side.

Weighting for lambda-z: **Uniform** (Phoenix default; matches both R engines).

## 6. Run

Click **Execute**.

## 7. Export

Object Browser → **Final Parameters Pivoted** → right-click → Export → CSV

Save as:

```
phoenix/reference/winnonlin_theoph_mine.csv
```

## 8. Pick it up in R

```r
source("R/04_engine_winnonlin.R")
wnl <- read_winnonlin_csv("phoenix/reference/winnonlin_theoph_mine.csv")
```

`run_all.R` finds any CSV in `phoenix/reference/` automatically.

---

## Recording provenance

For the comparison to mean anything six months later, record alongside the
export:

- Phoenix version (Help → About) and build number
- Licence tier and modules enabled
- The settings above, as actually set (screenshot the Options tab)
- Input file SHA-256

`outputs/provenance.txt` captures the R half; add the Phoenix half by hand or
via `phoenix/reference/README.md`.

---

## Expected result

Against the published WinNonlin 6.3/7.0 reference for Theophylline, both R
engines should agree to within 0.1% on every parameter. If your own Phoenix run
disagrees with the published reference, the cause is almost always one of:

1. **AUC method mismatch** — linear vs linear-log is the biggest single source
2. **Dose convention** — flat 320 mg vs subject-specific `Dose × Weight`
3. **Lambda-z point selection** — Best Fit vs a manually fixed range
4. **Cmax inclusion** — whether the peak is eligible for the terminal fit
5. **Phoenix version** — output column names and some defaults shifted across
   6.x → 8.x

Work down that list before suspecting the arithmetic.
