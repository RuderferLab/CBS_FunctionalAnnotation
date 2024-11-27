import pandas as pd
import numpy as np

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
ATSNP = snakemake.input[1:]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_atsnp_results(filepath: str) -> pd.DataFrame:
    """Reads atsnp results"""
    # Note; float_precision='round_trip'
    return pd.read_csv(
        filepath, sep="\t", engine="c", compression="gzip", float_precision="round_trip"
    )


def main():
    """d"""
    # Read results
    results = []
    for chrom_result in ATSNP:
        print(chrom_result)
        data = read_atsnp_results(chrom_result)
        results.append(data)

    # Concat
    results = pd.concat(results)

    # Drop allele and motif cols
    results.drop(["motif"], axis=1, inplace=True)

    # Rename fields
    results.rename(
        columns={
            "snpid": "vid",
            "pval_ref": "ref_pwm_pval",
            "pval_snp": "alt_pwm_pval",
            "pval_diff": "atsnp_pval_diff",
            "pval_rank": "atsnp_pval_rank",
            "snpbase": "alt",
            "log_lik_ref": "ref_pwm",
            "log_lik_snp": "alt_pwm",
            "log_lik_ratio": "dpwm",
        },
        inplace=True,
    )

    # Add log10 fields- add pseudocount if 0
    results["ref_pwm_pval_log10"] = -np.log10(results["ref_pwm_pval"].replace(0, 1e-9))
    results["alt_pwm_pval_log10"] = -np.log10(results["alt_pwm_pval"].replace(0, 1e-9))
    results["dpwm_rank_pval_log10"] = -np.log10(
        results["atsnp_pval_rank"].replace(0, 1e-9)
    )

    # Tidy up round
    results = results.round(5)

    # Write out
    fields = [
        "vid",
        "ref_pwm",
        "alt_pwm",
        "dpwm",
        "ref_pwm_pval_log10",
        "alt_pwm_pval_log10",
        "dpwm_rank_pval_log10",
    ]
    results[fields].to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
