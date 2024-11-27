## activity-scores.all-CBS.tsv

This file contains activity scores for CBS identified in the study. 

### Column labels

| field              | description |
| :---------------- | :------: | 
| rDHS | ENCODE accession of the rDHS |
| motif_ID | ID and information pertaining to the CTCF motif overlapped by the given rDHS. Some rDHS have more than 1 overlapping motif. In the study, we refer to these data as CBSs. |
| activity | CTCF binding activity score from meta-analyzing CTCF's Z-score at a given rDHS across 214 available biosamples. 
| activity_masked | Recalculation of activity with -10 scores masked. These scores correspond to rDHS-biosample pairs that held a raw-signal value of 0. 
| activity_masked_rDHS_quantile | Masked activity scores after binning into quantiles. Quantiles were used in the study to evaluate the functional properties of activity scores. |


## variant-scores.all-CBS.tsv

This file contains all CBS variants from gnomAD v3 assessed in the study. Each variant has been anntotated with its activity score and dPWM statistics. 


### Column labels

| field              | description |
| :---------------- | :------: | 
| VariantID | variant identity from gnomAD v3 |
| rDHS | ENCODE accession of the rDHS |
| dPWM | Observed dPWM score by substracting REF-ALT PWM values.
| dPWM-Class | LoB or loss of binding; GoB or gain of binding. Based on the sign of the dPWM statistic.
| dPWM-Sig | Flag whether the observed dPWM score was accompanied by a p-value < 0.05. Assessed using atSNP's importance sampling algorithm. 
| Activity | Masked Activity score.
| Activity_Quantile | Masked Activity quantile. |