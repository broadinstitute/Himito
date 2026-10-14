#!/usr/bin/env python3
"""Himito co-occurrence permutation statistic, re-implemented with options for calibration.
stat_i = sum_{j in O, j!=i} |A_i&A_j| / (|A_i&C_j| + |A_j&C_i| - |A_i&A_j|)
null: alt labels of row i shuffled among its covered reads (optionally within read-quality strata).
"""
import sys, csv, re, os, argparse, numpy as np
ap=argparse.ArgumentParser()
ap.add_argument('matrix'); ap.add_argument('out')
ap.add_argument('--select', default='snp', help='snp | all | file:<path>')
ap.add_argument('--draws', type=int, default=1000)
ap.add_argument('--extend-to', type=int, default=0)
ap.add_argument('--stop-after', type=int, default=10)
ap.add_argument('--others-min-ac', type=int, default=0)
ap.add_argument('--others-min-hf', type=float, default=0.0)
ap.add_argument('--floor-m', type=int, default=0)
ap.add_argument('--floor-h', type=float, default=0.0)
ap.add_argument('--partner-m', type=int, default=-1)
ap.add_argument('--partner-h', type=float, default=-1.0)
ap.add_argument('--ift', type=float, default=0.0, help='only test indels with HF>=ift')
ap.add_argument('--pois-eps', type=float, default=0.0, help='partners = SNPs with AC >= T(cov), T = Poisson(eps*cov) upper quantile')
ap.add_argument('--pois-tail', type=float, default=0.002)
ap.add_argument('--min-dist', type=int, default=0)
ap.add_argument('--others-snp-only', type=int, default=0)
ap.add_argument('--min-rq', type=float, default=0.0)
ap.add_argument('--strata', type=int, default=0, help='rq-quantile strata for the shuffle (0=off)')
ap.add_argument('--strata-by', default='rq', help='rq | burden')
ap.add_argument('--rq', default='read_rq.tsv')
ap.add_argument('--seed', type=int, default=42)
a=ap.parse_args()
rng=np.random.default_rng(a.seed)
cache=a.matrix+'.npz'
if os.path.exists(cache):
    z=np.load(cache, allow_pickle=True); A=z['A']; C=z['C']; names=list(z['names']); reads=list(z['reads'])
else:
    csv.field_size_limit(10**9); names=[]; Ar=[]; Cr=[]
    with open(a.matrix) as f:
        r=csv.reader(f); hdr=next(r); reads=hdr[1:]; R=len(reads)
        for row in r:
            v=row[1:]
            Cr.append(np.fromiter((x!='' for x in v),dtype=bool,count=R))
            Ar.append(np.fromiter((x!='' and float(x)>0.5 for x in v),dtype=bool,count=R))
            names.append(row[0])
    A=np.array(Ar); C=np.array(Cr); np.savez(cache, A=A, C=C, names=np.array(names,dtype=object), reads=np.array(reads,dtype=object))
rq={}
for l in open(a.rq):
    n,q,_=l.split('\t'); rq[n]=float(q)
rqv=np.array([rq.get(n,1.0) for n in reads])
if a.min_rq>0:
    keep=rqv>=a.min_rq; A=A[:,keep]; C=C[:,keep]; rqv=rqv[keep]; reads=[r for r,k in zip(reads,keep) if k]
nrow,R=A.shape
ac=A.sum(1); cov=C.sum(1); hf=np.where(cov>0, ac/np.maximum(cov,1), 0)
def is_snp(n):
    m=re.match(r'm\.(\d+)([ACGTN]+)>([ACGTN]+)',n); return len(m.group(2))==1 and len(m.group(3))==1
posv=np.array([int(re.match(r'm\.(\d+)',n).group(1)) for n in names])
snpmask=np.array([is_snp(n) for n in names])
floor_ok = ac >= np.maximum(a.floor_m, np.ceil(a.floor_h*cov))
import math
def pois_T(dp, eps, tail, m0=3):
    lam=eps*dp; k=0; p=math.exp(-lam); cdf=p
    while 1-cdf >= tail:
        k+=1; p*=lam/k; cdf+=p
    return max(m0, k+1)   # smallest k with P(X>=k) < tail
if a.pois_eps>0:
    T=np.array([pois_T(int(c), a.pois_eps, a.pois_tail) for c in cov])
    O=np.flatnonzero((ac>=T) & snpmask)
elif a.partner_m>=0:
    pfloor = ac >= np.maximum(a.partner_m, np.ceil(max(a.partner_h,0)*cov))
    O=np.flatnonzero(pfloor & snpmask)
else:
    O=np.flatnonzero((ac>=a.others_min_ac)&(hf>=a.others_min_hf)&(snpmask | (a.others_snp_only==0)))
AT=np.ascontiguousarray(A[O].T.astype(np.float32)); CT=np.ascontiguousarray(C[O].T.astype(np.float32))  # R x |O|
pos_in_O={j:k for k,j in enumerate(O)}
if a.select=='family': idx=[i for i,n in enumerate(names) if floor_ok[i] and (snpmask[i] or hf[i]>=a.ift)]
elif a.select=='snp': idx=[i for i,n in enumerate(names) if is_snp(n)]
elif a.select=='all': idx=list(range(nrow))
else: want=set(x.strip() for x in open(a.select[5:])); idx=[i for i,n in enumerate(names) if n in want]
# strata labels
if a.strata>0:
    if a.strata_by=='rq': key=rqv
    else: key=-A.sum(0).astype(float)  # per-read alt burden (higher burden -> lower key)
    qs=np.quantile(key, np.linspace(0,1,a.strata+1)[1:-1]); lab=np.searchsorted(qs,key)
else: lab=np.zeros(R,dtype=int)
print(f'rows={nrow} reads={R} others={len(O)} selected={len(idx)}', file=sys.stderr)
out=open(a.out,'w'); out.write('name\tac\tcov\thf\tobs\tnull_mean\tnull_sd\texceed\tdraws\tp\n')
for i in idx:
    ci=np.flatnonzero(C[i]); ai=np.flatnonzero(A[i]); na=len(ai)
    if na==0 or na==len(ci): continue
    h=na/len(ci)
    if h>0.2: continue                       # above heteroplasmic_frequency_threshold: not tested
    if is_snp(names[i]) and h>=0.10 and na>=20: continue   # COOC exemption
    K=(A[O].astype(np.float32) @ C[i].astype(np.float32))   # |A_j & C_i|
    self_k=pos_in_O.get(i,-1)
    d=np.abs(posv[O]-posv[i]); d=np.minimum(d,16569-d); near=np.flatnonzero(d<a.min_dist)
    # per-stratum pools
    pools=[]; 
    for s in np.unique(lab[ci]):
        pool=ci[lab[ci]==s]; k=int((lab[ai]==s).sum())
        if k>0: pools.append((pool,k))
    def stat_from_idx(SI):   # SI: B x na int indices
        X=AT[SI].sum(1); Y=CT[SI].sum(1)
        U=Y+K[None,:]-X; J=np.where(U>0, X/np.maximum(U,1), 0.0)
        if self_k>=0: J[:,self_k]=0
        if a.min_dist>0: J[:,near]=0
        return J.sum(1)
    obs=stat_from_idx(ai[None,:])[0]
    B=int(max(1,min(500, 4e7/(na*len(O)))))
    ex=0; done=0; ns=[]
    while done<a.draws and ex<a.stop_after:
        b=min(B,a.draws-done)
        SI=np.empty((b,na),dtype=np.int64); c0=0
        for pool,k in pools:
            # draw k without replacement per row: argsort of random keys
            rk=rng.random((b,len(pool))); sel=np.argpartition(rk, k-1, axis=1)[:,:k] if k<len(pool) else np.tile(np.arange(len(pool)),(b,1))
            SI[:,c0:c0+k]=pool[sel]; c0+=k
        z=stat_from_idx(SI); ns.append(z); ex+=int((z>=obs-1e-9).sum()); done+=b
    if a.extend_to>done and ex==0:
        while done<a.extend_to and ex<a.stop_after:
            b=min(B,a.extend_to-done)
            SI=np.empty((b,na),dtype=np.int64); c0=0
            for pool,k in pools:
                rk=rng.random((b,len(pool))); sel=np.argpartition(rk, k-1, axis=1)[:,:k] if k<len(pool) else np.tile(np.arange(len(pool)),(b,1))
                SI[:,c0:c0+k]=pool[sel]; c0+=k
            z=stat_from_idx(SI); ns.append(z); ex+=int((z>=obs-1e-9).sum()); done+=b
    nz=np.concatenate(ns)
    out.write(f"{names[i]}\t{na}\t{len(ci)}\t{na/len(ci):.6f}\t{obs:.4f}\t{nz.mean():.4f}\t{nz.std():.4f}\t{ex}\t{done}\t{(ex+1)/(done+1):.4g}\n"); out.flush()
