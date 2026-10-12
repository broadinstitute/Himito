#!/usr/bin/env python3
"""Offline re-implementation of Himito `call`'s co-occurrence permutation test.

Reads `<prefix>.raw_matrix.csv` (variant x read, written before strand-bias and
the permutation test) and re-runs the test so a change to the test can be
checked without rebuilding Himito. Statistic, null and BH are the same as
src/call.rs (sum of Jaccard against every other row; shuffle a row among its
covered reads; empirical p = (1 + #null >= obs) / (R + 1); BH over tested).

Modes
  --mode fixed     R = --rounds for every tested variant (today's behaviour)
  --mode adaptive  the Plan 2026-10-11 T1 algorithm (the ORACLE for T1):
                   stage 1 draws --rounds per variant, stopping a variant early
                   once it has --stop-after exceedances; stage 2 extends only the
                   variants still at zero exceedances to
                   R_target = ceil(headroom * m / (alpha * k_floor)), capped at
                   --max-rounds, repeating until R_target stops growing.
  --mode normal    z-score of the observation against its own null (mean, sd),
                   one-sided normal tail. Diagnostic only.
  --mode max-partner  ONT-recall T8 prototype: for each tested variant, over
                   partners with >= --min-shared co-alt reads, take the
                   hypergeometric upper-tail p on co-covered reads, multiply by
                   the partner count (Bonferroni), then BH as usual.

Validated 2026-10-11 on HapMap ONT r10 (HG00097.raw_matrix.csv): fixed R=1000
passes 0 truth (min q 0.0107 > 0.01); adaptive passes 54 truth / 0 non-truth,
matching the real `--permutation-rounds 10000` run (88 TP overall).

The column filter mirrors call.rs approximately: a row is tested when its
matrix alt frequency <= --test-hf (call.rs uses allele_count / graph depth,
which differs slightly). Rows above --null-hf get no null and are untested.

usage: perm_retest.py --matrix X.raw_matrix.csv --truth T.vcf.gz [--mode adaptive]
needs: numpy, pysam
"""
import argparse, csv, math, sys, time
import numpy as np
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--matrix", required=True)
ap.add_argument("--truth", required=True)
ap.add_argument("--mode", choices=["fixed", "adaptive", "normal", "max-partner"], default="adaptive")
ap.add_argument("--rounds", type=int, default=1000, help="fixed R, or stage-1 minimum")
ap.add_argument("--max-rounds", type=int, default=100000)
ap.add_argument("--stop-after", type=int, default=20, help="sequential stop: exceedances")
ap.add_argument("--headroom", type=float, default=10.0)
ap.add_argument("--alpha", type=float, default=0.01)
ap.add_argument("--test-hf", type=float, default=0.2)
ap.add_argument("--null-hf", type=float, default=0.7)
ap.add_argument("--min-shared", type=int, default=2, help="max-partner: min co-alt reads")
ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()

names, rows = [], []
with open(a.matrix) as fh:
    rd = csv.reader(fh); next(rd)
    for r in rd:
        names.append(r[0]); rows.append([float(x) if x != "" else np.nan for x in r[1:]])
M = np.array(rows)
cov = ~np.isnan(M)
alt = (np.nan_to_num(M) > 0.5) & cov
V, N = M.shape
af = alt.sum(1) / np.maximum(cov.sum(1), 1)
Y = alt.astype(np.float32); C = cov.astype(np.float32)
inter = Y @ Y.T; A = Y @ C.T
union = A + A.T - inter
J = np.where(union > 0, inter / np.maximum(union, 1), 0.0); np.fill_diagonal(J, 0)
obs = J.sum(1)
truth = {f"m.{r.pos}{r.ref}>{r.alts[0]}": float(r.samples[0]["HL"][0]) for r in pysam.VariantFile(a.truth)}
rng = np.random.default_rng(a.seed)

def draws(i, n):
    """n null draws of the summed Jaccard for row i."""
    idx = np.where(cov[i])[0]; vals = alt[i, idx]
    X = np.zeros((n, N), dtype=np.float32)
    for r in range(n):
        X[r, idx] = rng.permutation(vals)
    it = X @ Y.T; b = X @ C.T; c = (Y @ C[i])[None, :]
    u = b + c - it
    Jr = np.where(u > 0, it / np.maximum(u, 1), 0.0); Jr[:, i] = 0
    return Jr.sum(1)

tested = [i for i in range(V) if af[i] <= a.test_hf and af[i] <= a.null_hf]
m = len(tested)
state = {i: [0, 0, None] for i in tested}  # exceedances, draws, null list (normal mode)
t0 = time.time()

def run(i, target, chunk=200):
    e, d, _ = state[i]
    while d < target and e < a.stop_after:
        n = min(chunk, target - d)
        nd = draws(i, n)
        e += int((nd >= obs[i] - 1e-9).sum()); d += n
        if a.mode == "fixed" and d >= target:
            break
    state[i][0], state[i][1] = e, d

def hypergeom_sf(a, Npop, K, n):
    """P(X >= a) for X ~ Hypergeom(N=Npop, K=K, n=n)."""
    if a <= 0:
        return 1.0
    if Npop <= 0 or n <= 0 or K <= 0:
        return 1.0
    # survival: 1 - cdf(a-1)
    from math import lgamma
    def logC(n_, k_):
        if k_ < 0 or k_ > n_:
            return float("-inf")
        return lgamma(n_ + 1) - lgamma(k_ + 1) - lgamma(n_ - k_ + 1)
    log_den = logC(Npop, n)
    if log_den == float("-inf"):
        return 1.0
    s = 0.0
    lo = max(0, n + K - Npop)
    hi = min(a - 1, K, n)
    # cdf via summing pmf 0..a-1, then sf = 1 - cdf; for large a use upper sum
    if a - 1 - lo < (min(K, n) - (a - 1)):
        for x in range(lo, hi + 1):
            s += math.exp(logC(K, x) + logC(Npop - K, n - x) - log_den)
        return max(0.0, min(1.0, 1.0 - s))
    s = 0.0
    for x in range(max(a, lo), min(K, n) + 1):
        s += math.exp(logC(K, x) + logC(Npop - K, n - x) - log_den)
    return max(0.0, min(1.0, s))

def max_partner_p(i):
    """Bonferroni-adjusted min hypergeometric p over partners with >= min-shared co-alts."""
    partners = []
    for j in range(V):
        if j == i:
            continue
        both = int((alt[i] & alt[j]).sum())
        if both < a.min_shared:
            continue
        cocov = int((cov[i] & cov[j]).sum())
        n_i = int((alt[i] & cov[j]).sum())   # query alts among co-covered
        K_j = int((alt[j] & cov[i]).sum())   # partner alts among co-covered
        # both already counted on full; restrict to co-covered (same set)
        a_ij = both
        p_ij = hypergeom_sf(a_ij, cocov, K_j, n_i)
        partners.append((p_ij, j, a_ij, cocov))
    if not partners:
        return 1.0, 0, None
    partners.sort(key=lambda t: t[0])
    p_min, j_best, a_best, _ = partners[0]
    p_adj = min(1.0, p_min * len(partners))
    return p_adj, len(partners), (names[j_best], a_best, p_min)

if a.mode == "normal":
    p = {}
    for i in tested:
        nd = draws(i, a.rounds)
        z = (obs[i] - nd.mean()) / max(nd.std(), 1e-12)
        p[i] = 0.5 * math.erfc(z / math.sqrt(2))
elif a.mode == "max-partner":
    p = {}
    best = {}
    for i in tested:
        p[i], npart, info = max_partner_p(i)
        best[i] = (npart, info)
    print(f"max-partner: tested={m} min_shared={a.min_shared}", file=sys.stderr)
else:
    for i in tested:
        if a.mode == "fixed":
            nd = draws(i, a.rounds)
            state[i] = [int((nd >= obs[i] - 1e-9).sum()), a.rounds, None]
        else:
            run(i, a.rounds)
    if a.mode == "adaptive":
        target = a.rounds
        for _ in range(5):
            floor = [i for i in tested if state[i][0] == 0]
            if not floor:
                break
            new_target = min(a.max_rounds, math.ceil(a.headroom * m / (a.alpha * len(floor))))
            if new_target <= target:
                break
            target = new_target
            for i in floor:
                run(i, target)
        print(f"adaptive: final R_target {target}; total draws {sum(s[1] for s in state.values())}", file=sys.stderr)
    p = {i: (state[i][0] + 1) / (state[i][1] + 1) for i in tested}

ids = list(p); pv = np.array([p[i] for i in ids])
o = np.argsort(pv); q = np.empty(m); prev = 1.0
for rank in range(m, 0, -1):
    j = o[rank - 1]; prev = min(prev, pv[j] * m / rank); q[j] = prev
ok = q <= a.alpha
ist = np.array([names[i] in truth for i in ids])
print(f"{a.matrix}: mode={a.mode} rows={V} reads={N} tested={m} ({ist.sum()} truth) "
      f"time={time.time()-t0:.0f}s min_q={q.min():.4g}")
print(f"  pass at alpha {a.alpha}: truth {int((ok & ist).sum())}  non-truth {int((ok & ~ist).sum())}")
print("  non-truth passing: " + ", ".join(f"{names[i]}(hf {af[i]:.3f})" for k, i in enumerate(ids) if ok[k] and not ist[k]))
print("  truth failing: " + ", ".join(f"{names[i]}(HL {truth[names[i]]:.3g}, q {q[k]:.2g})" for k, i in enumerate(ids) if not ok[k] and ist[k]))
if a.mode == "max-partner":
    # Newly rescued vs adaptive baseline: truth that pass here
    print("  truth passing (name, q, best_partner, co-alts, raw_p, n_partners):")
    for k, i in enumerate(ids):
        if ok[k] and ist[k]:
            npart, info = best[i]
            if info is None:
                print(f"    {names[i]} q={q[k]:.3g} (no partner)")
            else:
                pj, a_ij, p_raw = info
                print(f"    {names[i]} q={q[k]:.3g} partner={pj} co={a_ij} p_raw={p_raw:.2g} npart={npart}")
