Created: 2024-09-04

**Last edited: 2025-04-16**

<h2>Authors</h2></summary>

Colby Tubbs, Mary Lauren Benton, Evonne McArthur, John A. Capra and Douglas M. Ruderfer @ the [Ruderfer Lab](https://ruderferlab.org)

<h2>Paper</h2></summary>

Tubbs, C., Benton, M. L., McArthur, E., Capra, J. A. & Ruderfer, D. M. Identifying deleterious noncoding variation through gain and loss of CTCF binding activity. [*The American Journal of Human Genetics*](doi:10.1016/j.ajhg.2025.02.009).


<h2>Publicly available data</h2></summary>
All datasets used in this study are publicly available. 

- Genetic data from gnomAD v3.1.2 was downloaded for all chromosomes as a Hail matrix table. |  [Link](https://gnomad.broadinstitute.org/data)


- CTCF’s binding activity Z-scores assessed in the cCRE framework were downloaded from the ENCODE-SCREEN web portal. | [Link](https://downloads.wenglab.org/cCREs/matrices/all/GRCh38.CTCF-zscore.rDHS-V3.txt.gz)

- CTCF’s canonical PWM (MA0139.2) was downloaded from the JASPAR 2024 web portal. | [Link](https://jaspar.elixir.no/api/v1/matrix/MA0139.2.jaspar)

- Conservation scores were downloaded from the UCSC genome browser. GERP++ scores in build hg19, Phastcons100 and PhyloP100 scores in build hg38 are found at:

    - GERP++ in build hg19 ([Link](https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw))

    - Phastcons100 in build hg38. | [Link](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way)

    - Phlop100 in build hg38. | [Link](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phylop100way)

    - LINSIGHT scores  in build hg19 from Huang et al. 2017. | [Link](http://compgen.cshl.edu/LINSIGHT/LINSIGHT.bw)

<h2>Repository code</h2></summary>

This respository is structured as a Snakemake workflow. To demonstrate how the study dataset was constructed, we've included all workflow input and outputs. All code used to generate these files are contained in the `workflow` directory and all intermediate outputs for this sample dataset are in the `results` directory.


For details on the generation of all figures included in the manuscript, see `workflow/notebooks/plots.ipynb`. For code to generate maps scores see `workflow/notebooks/maps.ipynb`.  We provide two summary files generated from the study; for all CBS annotated by their intersecting rDHS, CTCF motif and activity score see `results/Manuscript/files/activity-scores.all-CBS.tsv.gz`, for all CBS variants annotated by their activity and dPWM statistics see `results/Manuscript/files/variants-scores.all-CBS.tsv.gz`. 

<h2>Executing the pipeline</h2></summary>


To execute the pipeline without HPC resources, execute the command:

```snakemake --rerun-triggers mtime --use-conda --conda-prefix ~/snakemake_condas/ --use-conda --conda-prefix ~/snakemake_condas/ -c12 --rerun-incomplete```


To execute the pipeline with HPC resources, execute the command:
```snakemake --rerun-triggers mtime --use-conda --conda-prefix ~/snakemake_condas/ --jobs 25 --cluster-config "config/cluster.yaml" --cluster "sbatch --mem={resources.mem_mb} --time={resources.runtime} --cpus-per-task {threads} --output={cluster.output} --error={cluster.error} --account={cluster.account}"```