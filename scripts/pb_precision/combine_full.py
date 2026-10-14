#!/usr/bin/env python3
"""Emulate the proposed PacBio call path from a perm2 table over all AC>=min_ac rows.
usage: combine_full.py table out.vcf alpha min_ac indel_false_thr [himito.vcf] [rCRS.fasta]
himito.vcf = the -v 0 baseline VCF whose records/flags are re-filtered (default pb_benchmark/HG00097_pacbio.vcf)."""
import sys, re
tab, out, alpha, min_ac, ift = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
VCF = sys.argv[6] if len(sys.argv)>6 else 'pb_benchmark/HG00097_pacbio.vcf'
REF = sys.argv[7] if len(sys.argv)>7 else 'rCRS.fasta'
ref='N'+''.join(l.strip() for l in open(REF) if l[0]!='>').upper()
def norm(pos,r,a):
    while len(r)>0 and len(a)>0 and r[-1]==a[-1]:
        r=r[:-1]; a=a[:-1]
        if len(r)==0 or len(a)==0: pos-=1; r=ref[pos]+r; a=ref[pos]+a
    while len(r)>1 and len(a)>1 and r[0]==a[0]: r=r[1:]; a=a[1:]; pos+=1
    return (pos,r,a)
def parse(n):
    m=re.match(r'm\.(\d+)([ACGTN]+)>([ACGTN]+)',n); return int(m.group(1)),m.group(2),m.group(3)
rows=[l.rstrip('\n').split('\t') for l in open(tab)][1:]
rows=[r for r in rows if int(r[1])>=min_ac]
m=len(rows); ps=[float(r[9]) for r in rows]; order=sorted(range(m),key=lambda i:ps[i]); q=[1]*m; prev=1.0
for rank in range(m-1,-1,-1):
    i=order[rank]; prev=min(prev,ps[i]*m/(rank+1)); q[i]=prev
fail={norm(*parse(r[0])) for r,qq in zip(rows,q) if qq>alpha}
npass=0
with open(out,'w') as o:
    for l in open(VCF):
        if l[0]=='#': o.write(l); continue
        f=l.rstrip('\n').split('\t'); ad=int(f[9].split(':')[1]); hf=float(f[9].split(':')[2])
        if ad<min_ac: continue
        k=norm(int(f[1]),f[3],f[4])
        if k in fail: continue
        snp=len(f[3])==1 and len(f[4])==1
        if not snp:
            flt=[x for x in f[6].split(';') if x not in ('Potential_Artifact','PASS')]
            if abs(len(f[3])-len(f[4]))<5 and hf<ift: flt.append('Potential_Artifact')
            f[6]=';'.join(flt) if flt else 'PASS'
        npass+= f[6]=='PASS'
        o.write('\t'.join(f)+'\n')
print(f'tests={m} failed={len(fail)} pass={npass}',file=sys.stderr)
