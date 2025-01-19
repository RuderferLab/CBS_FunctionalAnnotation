## activity-scores.all-CBS.tsv

This file contains activity scores for CBS identified in the study. 

### Column labels

| field              | description |
| :---------------- | :------: | 
| rDHS | ENCODE accession of the rDHS |
| mid | Motif ID used in the study. We refer to these motifs as CBSs. |
| cid | CBS ID used in study |
| activity | CTCF binding activity score from meta-analyzing CTCF's Z-score at a given rDHS across 214 available biosamples. 
| quantile | Masked activity scores after binning into quantiles. Quantiles were used in the study to evaluate the functional properties of activity scores. |


## variant-scores.all-CBS.tsv

This file contains all CBS variants from gnomAD v3 assessed in the study. Each variant has been anntotated with its activity score and dPWM statistics. 


### Column labels

| field              | description |
| :---------------- | :------: | 
| vid | variant identity from gnomAD v3 (chrom-pos1-ref-alt) |
| rDHS | ENCODE accession of the rDHS |
| mid | motif ID used in study |
| cid | CBS ID used in study |
| dpwm | Observed dPWM score by substracting REF-ALT PWM values.
| dpwm_class | L or loss of binding; G or gain of binding. Based on the sign of the dpwm statistic.
| dpwm_sig | Flag whether the observed dPWM score was accompanied by a p-value < 0.05. Assessed using atSNP's importance sampling algorithm. 
| activity | Masked activity score.
| quantile | Masked activity quantile. |