#!/bin/bash
# strat2.sh evaldir -> "P Padj S F | classes"
python3 - "$1" <<'PY'
import sys,gzip,collections
D=sys.argv[1]
corr={(9342,'A'),(15622,'C'),(10092,'T'),(12731,'C'),(16162,'G'),(15150,'A'),(5532,'A')}
def cls(h): return 'maj' if h>=0.5 else '10%' if h>=0.08 else '2%' if h>=0.015 else '0.5%'
tp=collections.Counter(); fn=collections.Counter()
for fname,ctr in [('tp-baseline.vcf.gz',tp),('fn.vcf.gz',fn)]:
    for l in gzip.open(f'{D}/{fname}','rt'):
        if l[0]=='#': continue
        f=l.split('\t'); ctr[cls(float(f[9].split(':')[1]))]+=1
tpc=sum(1 for l in gzip.open(f'{D}/tp.vcf.gz','rt') if l[0]!='#')
fps=[l.split('\t') for l in gzip.open(f'{D}/fp.vcf.gz','rt') if l[0]!='#']
fp=len(fps); fpadj=sum(1 for f in fps if (int(f[1]),f[4]) not in corr)
fpi=sum(1 for f in fps if not(len(f[3])==1 and len(f[4])==1))
T=sum(tp.values()); N=T+sum(fn.values())
P=tpc/(tpc+fp) if tpc+fp else float('nan'); Pa=tpc/(tpc+fpadj) if tpc+fpadj else float('nan'); S=T/N; F=2*P*S/(P+S) if P+S else 0
print(f"P={P:.3f} Padj={Pa:.3f} S={S:.3f} F={F:.3f} FP={fp}(indel {fpi}) | "+' '.join(f"{c}:{tp[c]}/{tp[c]+fn[c]}" for c in ['maj','10%','2%','0.5%']))
PY
