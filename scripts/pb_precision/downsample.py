#!/usr/bin/env python3
"""downsample.py in.bam out.bam target_depth seed : keep whole reads (all records of a qname) by hashed qname."""
import sys, hashlib, numpy as np, pysam
src,out,d,seed=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
mean=np.array(pysam.AlignmentFile(src).count_coverage('chrM',0,16569,quality_threshold=0)).sum(0).mean()
f=d/mean
h=lambda q: int(hashlib.md5(f'{seed}:{q}'.encode()).hexdigest()[:8],16)/0xffffffff
b=pysam.AlignmentFile(src); o=pysam.AlignmentFile(out,'wb',template=b)
for a in b.fetch('chrM'):
    if h(a.query_name)<f: o.write(a)
o.close(); pysam.index(out)
