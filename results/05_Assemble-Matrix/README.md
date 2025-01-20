## final_matrix.tsv.gz

This matrix is the final product of the workflow and is used as input in one form or another for all subsequent analyses. The matrix is organized by single nucleotide positions in the reference genome (rows) with corresponding annotations as columns. The positions here correspond to all hg38 positions contained within all CBS identified in the study. 

Note on the frequency of NaNs: Not all annotations exist for all positions contained within the matrix. For example, if there is no observed variant from gnomAD at a given position, then all variant-based annotations (deltaPWM, CADD) would be NaN.

### Column labels

| field              | description |
| :---------------- | :------: | 
| pid | "position ID" of a given nucleotide in hg38 |
| mid | "motif ID" used in the study |
| vid | Variant ID of an intersecting SNV in gnomAD v3. NaN if no high-quality variant exists at the given position. |
| af | Allele frequency of an intersecting SNV in gnomAD v3. NaN if no high-quality variant exists at the given position. |
| ac | Allele count an intersecting SNV in gnomAD v3. NaN if no high-quality variant exists at the given position. |
| an | Allele number of an intersecting SNV in gnomAD v3. NaN if no high-quality variant exists at the given position. |
| caddPhred | PHRED-Scaled CADD score of a given variant. NaN if no high-quality variant exists at the given position in gnomAD v3. See https://cadd.bihealth.org/info for the differences between scaled and raw CADD scores. Note: PHRED-Scaled CADD scores were used to assess constraint. |
| singleton | Flags if an intersecting variant has an allele count of 1. Used in the calculation of MAPS scores. NaN if no variant exists. |
| isvar | Flags if a position has an instersecting variant from gnomAD. |
| context | The trinucleotide context of the position based on hg38. |
| score | The PWM score of the CTCF motif relative to the maximum score possible under the PWM. |
| pval | Statistical significance of the observed score under rel_pwm. |
| rDHS | The accession code of the position's intersecting rDHS.  |
| activity | The meta-analyzed activity score of the rDHS, after masking scores of -10. |
| quantile_rdhswide | Masked activity quantile computed across all known cCREs. |
| gerp++ | GERP++ score for the position in hg19. Scores were obtained by lifting back the position's coordinates to hg19. NaN corresponds to positions that couldn't be scored. |
| phastcons100 | Phastcons100 score for the position. NaN corresponds to positions that couldn't be scored. |
| linsight | LINSIGHT score for the position in hg19. Scores were obtained by lifting back the position's coordinates to hg19. NaN corresponds to positions that couldn't be scored. |
| phylop100 | Phastcons100 score for the position. NaN corresponds to positions that couldn't be scored. |
| dpwm | Change in PWM score calculated by atSNP (REF-ALT) |
| dpwm_log10p | -log10P-value of the observed dPWM score calculated through atSNP's importance sampling algorithm. |
| dpwm_class | L or loss of binding; G or gain of binding. Based on the sign of the dpwm statistic. |