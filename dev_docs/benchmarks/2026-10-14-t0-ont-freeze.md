# T0 ONT golden freeze (plan 2026-10-13-pacbio-precision)
commit: 4d07968a43e3f45e75b4ba5212894df7d85f54fd
binary: target/release/Himito (downsample_list fix included)
frozen: 2026-10-14T18:30:50Z

## quick-start HapMap_ont_merge.bam
### ont-r10
quick-start: (rescored from frozen VCF)
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None                115            114          5          6     0.9580       0.9504     0.9542
== _t0_r10.vcf  (truth n=121)
   TP 114  FP 5  FN 7  recall 0.9421  precision 0.9580
   by truth HL (TP, FN): 1.0:(6, 1), 0.98:(1, 0), 0.96:(5, 0), 0.935:(1, 0), 0.86:(1, 0), 0.855:(2, 0), 0.84:(6, 0), 0.835:(15, 1), 0.105:(3, 0), 0.1:(17, 0), 0.045:(1, 0), 0.025:(2, 0), 0.02:(42, 5), 0.005:(12, 0)
   FP: 302A>AC, 9342G>A, 10092C>T, 12731T>C, 15622T>C
   truth called but not PASS: 

### ont-denoised
quick-start: 678s total; strand-bias + permutation step: 12s
During ALT comparison no ROC data will be produced, as a sample is required by the selected ROC score field: GQ (FORMAT)
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None                110            110          4         11     0.9649       0.9091     0.9362
== regress_ont-denoised.vcf  (truth n=121)
   TP 110  FP 4  FN 11  recall 0.9091  precision 0.9649
   by truth HL (TP, FN): 1.0:(6, 1), 0.98:(0, 1), 0.96:(5, 0), 0.935:(1, 0), 0.86:(1, 0), 0.855:(2, 0), 0.84:(6, 0), 0.835:(15, 1), 0.105:(2, 1), 0.1:(17, 0), 0.045:(0, 1), 0.025:(2, 0), 0.02:(41, 6), 0.005:(12, 0)
   FP: 9342G>A, 10092C>T, 12731T>C, 15622T>C
   truth called but not PASS: 310T>TC(HL 0.98, Strand_bias)
