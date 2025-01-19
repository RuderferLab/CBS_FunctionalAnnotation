"""jaspar-motifs.py"""

import pandas as pd
import pyranges as pr

# I/O
EXONS = "resources/data/gencode/gencode.hg38.exons.pc.bed.gz"
ENCODE_BLACKLIST = "resources/data/encode/hg38.blacklist.bed.gz"
HG38_GAPS = "resources/data/genome/hg38.gaps.bed.gz"
MOTIF_TRACK = "resources/data/jaspar/MA0139.2.tsv.gz"
OUTPUT = "resources/files/CTCF-motifs.MA01392.hg38.tsv.processed.gz"

# Parameters
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ------------- #
# Functions     #
# ------------- #


def read_pyranges(filepath: str) -> pr.PyRanges:
    """Returns bed on file as PyRanges object"""
    return pr.read_bed(filepath)

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
    """Main"""
    # Read in datasets as PyRanges
    exons = read_pyranges(EXONS)
    encode_blacklist = read_pyranges(ENCODE_BLACKLIST)
    hg38_gaps = read_pyranges(HG38_GAPS)

    # Read in motif track as DataFrame
    track = read_track(MOTIF_TRACK)

    # Subset to primary chromosomes
    track = track[track["Chromosome"].isin(CHROMOSOMES)]

    # Create pyranges object
    track = pr.PyRanges(track)

    # Exclude exons, encode blacklist regions, and hg38 gaps
    track = track.overlap(exons, invert=True).overlap(encode_blacklist, invert=True).overlap(hg38_gaps, invert=True)

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
    track.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
