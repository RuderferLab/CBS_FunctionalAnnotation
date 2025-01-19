import pandas as pd
import numpy as np

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
ATSNP = snakemake.input[1:]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    return pd.read_csv(
        filepath, sep="\t", engine="c", usecols=["vid"], dtype={"vid": str}
    )


def read_atsnp_results(filepath: str) -> pd.DataFrame:
    """Reads atsnp results"""
    # Note; float_precision='round_trip'
    fields = ["snpid", "log_lik_ratio", "pval_rank"]
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        compression="gzip",
        float_precision="round_trip",
        usecols=fields,
    )


def main():
    # Read track
    track = read_track(TRACK)

    # Select variants
    track = track[track["vid"] != "."]

    # Read results
    results = []
    for chrom_result in ATSNP:
        print(chrom_result)
        data = read_atsnp_results(chrom_result)
        data.columns = ["vid"] + ["pval", "dpwm"]
        results.append(data)

    # Concat
    results = pd.concat(results)

    # Add log10 fields- add pseudocount if 0
    results["dpwm_log10p"] = -np.log10(results["pval"].replace(0, 10e-7))
    results["dpwm_class"] = results["dpwm"].apply(lambda x: "L" if x > 0 else "G")

    # Subset to these fields
    fields = ["vid", "dpwm", "dpwm_log10p", "dpwm_class"]
    results = results[fields]

    # Merge with track of all variants
    track = track.merge(results, on="vid", how="left")

    # Not all variants have atsnp results, address missingness
    missingness_scheme = {"dpwm": np.nan, "dpwm_log10p": np.nan, "dpwm_class": "."}
    track = track.fillna(missingness_scheme)

    # Ensure dtypes
    dtypes = {"vid": str, "dpwm": float, "dpwm_log10p": float, "dpwm_class": str}
    track = track.astype(dtypes)

    # Write out
    track.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
