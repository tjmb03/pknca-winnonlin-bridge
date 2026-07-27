"""
Independent verification of the NCA comparison.

Implements NCA from first principles in pure Python -- no PKNCA, no NonCompart,
no R. Compares against the Phoenix WinNonlin reference CSV.

If this agrees with WinNonlin, then the R engines agreeing with WinNonlin is
corroborated by a third, structurally unrelated implementation.
"""
import csv
import math
from collections import defaultdict

DOSE = 320.0

# ---------------------------------------------------------------- load ----
rows = defaultdict(list)
with open("/tmp/theoph_raw.csv") as f:
    for r in csv.DictReader(f):
        rows[int(r["Subject"])].append((float(r["Time"]), float(r["conc"])))
for s in rows:
    rows[s].sort()


# ------------------------------------------------------------ integrate ----
def auc_aumc_linup_logdown(t, c):
    """Linear-up / log-down trapezoidal AUC and AUMC.

    Descending segments with both concentrations > 0 use the log rule;
    everything else (ascending, or a zero endpoint) falls back to linear.
    This is Phoenix's 'Linear Log Trapezoidal'.
    """
    auc = 0.0
    aumc = 0.0
    for i in range(len(t) - 1):
        t1, t2 = t[i], t[i + 1]
        c1, c2 = c[i], c[i + 1]
        dt = t2 - t1
        if c2 < c1 and c1 > 0 and c2 > 0:
            lr = math.log(c1 / c2)              # positive
            auc += dt * (c1 - c2) / lr
            # AUMC log rule. Derived from integrating t*c(t) for
            # c(t) = c1*exp(-k(t-t1)):
            #   dt*(t2*c2 - t1*c1)/L  -  dt^2*(c2 - c1)/L^2,  L = ln(c2/c1)
            # The dt and dt^2 factors are easy to drop -- doing so leaves AUC
            # correct and AUMC wrong, which is exactly the signature seen on
            # the first run of this script.
            L = math.log(c2 / c1)               # negative
            aumc += dt * (t2 * c2 - t1 * c1) / L - dt**2 * (c2 - c1) / (L ** 2)
        else:
            auc += dt * (c1 + c2) / 2.0
            aumc += dt * (t1 * c1 + t2 * c2) / 2.0
    return auc, aumc


# ------------------------------------------------------------- lambda-z ----
def loglin_fit(t, c):
    """OLS of ln(c) on t. Returns (slope, r2, adj_r2, n)."""
    n = len(t)
    y = [math.log(v) for v in c]
    mt = sum(t) / n
    my = sum(y) / n
    sxy = sum((a - mt) * (b - my) for a, b in zip(t, y))
    sxx = sum((a - mt) ** 2 for a in t)
    slope = sxy / sxx
    inter = my - slope * mt
    ss_res = sum((b - (inter + slope * a)) ** 2 for a, b in zip(t, y))
    ss_tot = sum((b - my) ** 2 for b in y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    adj = 1 - (1 - r2) * (n - 1) / (n - 2) if n > 2 else float("nan")
    return slope, r2, adj, n


def best_fit_lambda_z(t, c, tmax, min_pts=3, tol=1e-4):
    """Phoenix 'Best Fit' terminal slope.

    Candidate windows are the last k points for k = 3, 4, ...; points at or
    before Tmax are excluded. The window with the highest adjusted R-squared
    wins; a longer window whose adjusted R-squared is within 1e-4 of the best
    is preferred (WinNonlin's tie-break toward more points).
    """
    idx = [i for i in range(len(t)) if t[i] > tmax and c[i] > 0]
    if len(idx) < min_pts:
        return None

    best = None
    for k in range(min_pts, len(idx) + 1):
        w = idx[-k:]
        tt = [t[i] for i in w]
        cc = [c[i] for i in w]
        slope, r2, adj, n = loglin_fit(tt, cc)
        if slope >= 0:
            continue
        cand = dict(lamz=-slope, r2=r2, adj=adj, n=n,
                    lower=tt[0], upper=tt[-1])
        if best is None:
            best = cand
        elif cand["adj"] > best["adj"] + tol:
            best = cand
        elif cand["adj"] > best["adj"] - tol and cand["n"] > best["n"]:
            best = cand
    return best


# ------------------------------------------------------------- per-subj ----
results = {}
for s in sorted(rows):
    t = [p[0] for p in rows[s]]
    c = [p[1] for p in rows[s]]

    cmax = max(c)
    tmax = t[c.index(cmax)]

    meas = [i for i in range(len(t)) if c[i] > 0]
    tlast, clast = t[meas[-1]], c[meas[-1]]

    auclast, aumclast = auc_aumc_linup_logdown(t, c)

    lz = best_fit_lambda_z(t, c, tmax)
    if lz is None:
        continue
    lamz = lz["lamz"]
    hl = math.log(2) / lamz

    aucinf = auclast + clast / lamz
    aumcinf = aumclast + tlast * clast / lamz + clast / (lamz ** 2)
    pext = 100 * (aucinf - auclast) / aucinf
    clf = DOSE / aucinf
    vzf = DOSE / (aucinf * lamz)
    mrtinf = aumcinf / aucinf

    results[s] = dict(
        CMAX=cmax, TMAX=tmax, TLST=tlast, CLST=clast,
        AUCLST=auclast, AUMCLST=aumclast, LAMZ=lamz, LAMZHL=hl,
        LAMZNPT=lz["n"], R2ADJ=lz["adj"], LAMZLL=lz["lower"], LAMZUL=lz["upper"],
        AUCIFO=aucinf, AUMCIFO=aumcinf, AUCPEO=pext, CLFO=clf, VZFO=vzf,
        MRTEVIFO=mrtinf,
    )

# ------------------------------------------------------------- compare ----
WNL_COL = {
    "CMAX": "Cmax", "TMAX": "Tmax", "TLST": "Tlast", "CLST": "Clast",
    "AUCLST": "AUClast", "AUMCLST": "AUMClast", "LAMZ": "Lambda_z",
    "LAMZHL": "HL_Lambda_z", "LAMZNPT": "No_points_lambda_z",
    "R2ADJ": "Rsq_adjusted", "LAMZLL": "Lambda_z_lower",
    "LAMZUL": "Lambda_z_upper", "AUCIFO": "AUCINF_obs",
    "AUMCIFO": "AUMCINF_obs", "AUCPEO": "AUC_%Extrap_obs",
    "CLFO": "Cl_F_obs", "VZFO": "Vz_F_obs", "MRTINF": "MRTINF_obs",
    "MRTEVIFO": "MRTINF_obs",
}

wnl = {}
with open("phoenix/reference/winnonlin_theoph_log.csv") as f:
    for r in csv.DictReader(f):
        wnl[int(r["Subject"])] = r

print("Independent Python NCA  vs  Phoenix WinNonlin reference")
print("=" * 74)

worst = 0.0
worst_where = ""
n_comp = 0
n_fail = 0

for p in ["CMAX", "TMAX", "TLST", "CLST", "LAMZNPT", "LAMZLL", "LAMZUL",
          "LAMZ", "LAMZHL", "R2ADJ", "AUCLST", "AUCIFO", "AUCPEO",
          "AUMCLST", "AUMCIFO", "CLFO", "VZFO", "MRTEVIFO"]:
    diffs = []
    for s in sorted(results):
        mine = results[s][p]
        ref = float(wnl[s][WNL_COL[p]])
        n_comp += 1
        if abs(ref) < 1e-12:
            d = 0.0 if abs(mine) < 1e-12 else float("inf")
        else:
            d = 100 * (mine - ref) / ref
        diffs.append(abs(d))
        if abs(d) > 0.1:
            n_fail += 1
        if abs(d) > worst:
            worst, worst_where = abs(d), f"{p} subj {s}"
    mx = max(diffs)
    flag = "OK" if mx < 0.1 else "**MISMATCH**"
    print(f"  {p:<10s} n={len(diffs):2d}  max |diff| = {mx:.3e} %   {flag}")

print("=" * 74)
print(f"Comparisons: {n_comp}   Exceeding 0.1%: {n_fail}")
print(f"Worst: {worst:.3e} %  ({worst_where})")
print()
print("Subject 1, computed here vs WinNonlin:")
for p in ["AUCLST", "AUCIFO", "LAMZHL", "CLFO", "VZFO", "AUMCIFO", "MRTEVIFO"]:
    print(f"  {p:<9s} {results[1][p]:>16.8f}   {float(wnl[1][WNL_COL[p]]):>16.8f}")
