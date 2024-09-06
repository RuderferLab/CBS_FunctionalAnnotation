# Run pipeline, locally:
snakemake --rerun-triggers mtime --use-conda --conda-prefix ~/snakemake_condas/ --use-conda --conda-prefix ~/snakemake_condas/ -c12 --rerun-incomplete -np


# Run pipeline, clutser
snakemake --rerun-triggers mtime --use-conda --conda-prefix ~/snakemake_condas/ --jobs 25 --cluster-config "config/cluster.yaml" --cluster "sbatch --mem={resources.mem_mb} --time={resources.runtime} --cpus-per-task {threads} --output={cluster.output} --error={cluster.error} --account={cluster.account}" -np
