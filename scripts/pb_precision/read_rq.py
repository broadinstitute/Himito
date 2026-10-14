#!/usr/bin/env python3
"""read_rq.py mt.bam > read_rq.tsv : PacBio rq/np per primary read (needs pysam)."""
import sys, pysam
for a in pysam.AlignmentFile(sys.argv[1]).fetch(until_eof=True):
    if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
    print(f"{a.query_name}\t{a.get_tag('rq') if a.has_tag('rq') else 1.0}\t{a.get_tag('np') if a.has_tag('np') else -1}")
