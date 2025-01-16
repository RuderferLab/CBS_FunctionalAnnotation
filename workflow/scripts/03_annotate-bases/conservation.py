"""d"""

import pyBigWig
from pyliftover import LiftOver
import pandas as pd
import numpy as np

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
LINSIGHT = snakemake.input[1]  # type: ignore
PHASTCONS = snakemake.input[2]  # type: ignore
GERP = snakemake.input[3]  # type: ignore
PHYLOP = snakemake.input[4]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    """Returns track as DataFrame"""
    dtypes = {"pid": str, "mid": str, "mid_idx": int}
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        names=dtypes.keys(),
        dtype=dtypes,
    )


def convert_coord(row, lo_chain) -> int:
    """D"""
    try:
        return lo_chain.convert_coordinate(row["chrm"], row["pos1"])[0][1]
    except:
        return "NaN"


def main():
    """d"""
    print("Processing conservation scores...")

    # Setup lift chain
    lo = LiftOver("hg38", "hg19")

    # Connect to bigwigs
    linsight_bw = pyBigWig.open(LINSIGHT)
    phast_bw = pyBigWig.open(PHASTCONS)
    gerp_bw = pyBigWig.open(GERP)
    phylop_bw = pyBigWig.open(PHYLOP)

    # Read in track
    track = read_track(TRACK)
    
    # Add hg38 coords from pid
    track[['chrm', 'pos1']] = track['pid'].str.split("-", expand=True)
    track['pos1'] = track['pos1'].astype(int)
    track['pos0'] = track['pos1'] - 1

    # Add hg19 coords for linsight bw
    track["pos1_hg19"] = track.apply(lambda x: convert_coord(x, lo), axis=1)

    # Initialize database to store the scores
    score_dbase = []

    print("Fetching scores...")
    # Loop over variants (row) and gets scores
    for row in track.itertuples():

        # Row information
        chrom = row.chrm
        start = row.pos0
        end = row.pos1
        pid = row.pid
        mid = row.mid

        # search phastcons with regular chrom
        phast = phast_bw.values(chrom, start, end)[0]

        # Search phylop with regular chrom
        phylop = phylop_bw.values(chrom, start, end)[0]

        # Search linsight with regular chrom and hg19 coords
        if row.pos1_hg19 != "NaN":
            end_hg19 = int(row.pos1_hg19)
            start_hg19 = end_hg19 - 1
            try:
                linsight = linsight_bw.values(chrom, start_hg19, end_hg19)[0]
            except RuntimeError:
                linsight = np.nan
            try:
                gerp = gerp_bw.values(chrom, start_hg19, end_hg19)[0]
            except RuntimeError:
                gerp = np.nan
        else:
            gerp = np.nan
            linsight = np.nan

        # Store results
        score_dbase.append([pid, gerp, phast, linsight, phylop])

    # Update fields
    score_dbase = pd.DataFrame(
        score_dbase,
        columns=[
            "pid",
            "gerp++",
            "phastcons100",
            "linsight",
            "phylop100",
        ],
    )
    
    # Tidy digits
    score_dbase = score_dbase.round(4)
    
    # Write out
    score_dbase.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
