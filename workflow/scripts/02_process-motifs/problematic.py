import pandas as pd
import pyranges as pr

# Snakemake
TRACK = snakemake.input[0]  # type:ignore
OUTPUT = snakemake.output[0]  # type:ignore

# Globals
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    dtype = {
        "Chromosome": str,
        "Start": int,
        "End": int,
        "name": str,
        "score": float,
        "pval": float,
        "strand": str,
    }
    return pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=list(dtype.keys()),
        dtype=dtype,
        engine="c",
    )


def main():
    # Create pyranges object
    track = read_track(TRACK)

    # Subset to primary chromosomes
    track = track[track["Chromosome"].isin(CHROMOSOMES)]

    # Create pyranges object
    track = pr.PyRanges(track)

    # Cluster overlapping motifs
    track = track.cluster(slack=0).as_df()

    # Select highese scoring motif per cluster
    track.sort_values(["Cluster", "score"], ascending=[True, False], inplace=True)
    track.drop_duplicates("Cluster", keep="first", inplace=True)

    # Add motif id
    track.reset_index(drop=True, inplace=True)
    track["mid"] = "M" + (track.index + 1).astype(str).str.zfill(7)

    # Clean up
    track.drop(columns=["Cluster", "name"], inplace=True)
    dtypes = {
        "Chromosome": str,
        "Start": int,
        "End": int,
        "score": int,
        "pval": int,
        "strand": str,
        "mid": str,
    }
    track = track.astype(dtypes)
    
    # Save
    track.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
