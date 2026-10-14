"""Background substitution spectrum: usage noise.py mt.bam rCRS.fasta truth.vcf.gz [min_rq]"""
import pysam, numpy as np, gzip, sys
BAM,REF,TRUTH=sys.argv[1],sys.argv[2],sys.argv[3]
min_rq=float(sys.argv[4]) if len(sys.argv)>4 else 0
ref=''.join(l.strip() for l in open(REF) if l[0]!='>').upper()
bam=pysam.AlignmentFile(BAM)
def rf(read):
    return (not read.is_secondary) and read.get_tag('rq')>=min_rq
cnt=np.array(bam.count_coverage('chrM',0,16569,quality_threshold=0,read_callback=rf))  # 4 x L  (A,C,G,T)
truth=set()
for l in gzip.open(TRUTH,'rt'):
    if l[0]=='#':continue
    f=l.split('\t')
    if len(f[3])==1 and len(f[4])==1: truth.add((int(f[1]),f[4]))
extra={(9342,'A'),(15622,'C'),(10092,'T'),(12731,'C'),(16162,'G'),(15150,'A'),(4357,'C'),(5532,'A')}
B='ACGT'; dp=cnt.sum(0); noise=[]
for p in range(16569):
    major=cnt[:,p].argmax()
    for b in range(4):
        if b==major or B[b]==ref[p]: continue
        if (p+1,B[b]) in truth or (p+1,B[b]) in extra: continue
        noise.append((cnt[b,p], dp[p]))
noise=np.array(noise); ac=noise[:,0]; 
print(f'min_rq={min_rq} alleles={len(ac)} median DP={np.median(noise[:,1]):.0f} mean noise AC={ac.mean():.3f} (rate {ac.sum()/noise[:,1].sum():.2e})')
for k in [3,5,10,15,20,25,30,40,50]:
    print(f'  noise alleles with AC>={k}: {(ac>=k).sum()}')
o=np.argsort(-ac.astype(np.int64))[:12]; print('  top noise AC:', [int(ac[i]) for i in o])
