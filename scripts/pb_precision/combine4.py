#!/usr/bin/env python3
"""combine4.py table cand.vcf out.vcf alpha mode ift [floor_h] [eps] [ift_in_tracts]
(set HIMITO_REF=path/to/rCRS.fasta; cand.vcf = `call -p 1 --permutation-rounds 1 -v 0` output; table = perm3.py output)
mode=all : every SNP with AC>=max(3,ceil(floor_h*DP)) is linkage-tested
mode=gray: SNPs with AC>=T(DP) (Poisson eps=0.002, tail=0.002) are exempt; 3<=AC<T are tested
Indels: tested if HF>=ift; PASS iff HF>=ift (len<5) and not failing."""
import sys, re, math
tab, VCF, out, alpha, mode, ift = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5], float(sys.argv[6])
fh = float(sys.argv[7]) if len(sys.argv)>7 else 0.0
import os
ref='N'+''.join(l.strip() for l in open(os.environ.get('HIMITO_REF','rCRS.fasta')) if l[0]!='>').upper()
EPS=float(sys.argv[8]) if len(sys.argv)>8 else 0.002
IFT_TRACT=float(sys.argv[9]) if len(sys.argv)>9 else None
TRACTS=[(301,316),(512,526),(16181,16194)]
def ift_at(pos):
    if IFT_TRACT is not None and any(a<=pos<=b for a,b in TRACTS): return IFT_TRACT
    return ift
def pois_T(dp, eps=None, tail=0.002, m0=3):
    eps = EPS if eps is None else eps
    lam=eps*dp; k=0; p=math.exp(-lam); cdf=p
    while 1-cdf >= tail:
        k+=1; p*=lam/k; cdf+=p
    return max(m0, k+1)
def norm(pos,r,a):
    while len(r)>0 and len(a)>0 and r[-1]==a[-1]:
        r=r[:-1]; a=a[:-1]
        if len(r)==0 or len(a)==0: pos-=1; r=ref[pos]+r; a=ref[pos]+a
    while len(r)>1 and len(a)>1 and r[0]==a[0]: r=r[1:]; a=a[1:]; pos+=1
    return (pos,r,a)
def parse(n):
    m=re.match(r'm\.(\d+)([ACGTN]+)>([ACGTN]+)',n); return int(m.group(1)),m.group(2),m.group(3)
def floor(dp): return max(3, math.ceil(fh*dp))
rows=[]
for l in list(open(tab))[1:]:
    r=l.rstrip('\n').split('\t'); n=r[0]; p,ra,aa=parse(n); snp=len(ra)==1 and len(aa)==1
    ac_,cov_,hf_=int(r[1]),int(r[2]),float(r[3])
    if snp:
        if ac_<floor(cov_): continue
        if mode=='gray' and ac_>=pois_T(cov_): continue
    else:
        if hf_<ift_at(p) or ac_<3: continue
    rows.append((n,float(r[9])))
m=len(rows); order=sorted(range(m),key=lambda i:rows[i][1]); q=[1]*m; prev=1.0
for rank in range(m-1,-1,-1):
    i=order[rank]; prev=min(prev,rows[i][1]*m/(rank+1)); q[i]=prev
fail={norm(*parse(r[0])) for r,qq in zip(rows,q) if qq>alpha}
with open(out,'w') as o:
    for l in open(VCF):
        if l[0]=='#': o.write(l); continue
        f=l.rstrip('\n').split('\t'); ad=int(f[9].split(':')[1]); hf=float(f[9].split(':')[2])
        dp=int(f[7].split('DP=')[1].split(';')[0])
        if ad<3: continue
        snp=len(f[3])==1 and len(f[4])==1
        if snp and ad<floor(dp): continue
        if norm(int(f[1]),f[3],f[4]) in fail: continue
        if not snp:
            flt=[x for x in f[6].split(';') if x not in ('Potential_Artifact','PASS')]
            if abs(len(f[3])-len(f[4]))<5 and hf<ift_at(int(f[1])): flt.append('Potential_Artifact')
            f[6]=';'.join(flt) if flt else 'PASS'
        o.write('\t'.join(f)+'\n')
print(f'tests={m} failed={len(fail)}',file=sys.stderr)
