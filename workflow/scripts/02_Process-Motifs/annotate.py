import pandas as pd
import pybedtools as pbt

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
ACTIVITY = snakemake.input[1]  # type: ignore
MATRIX_OP = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_bed(filepath: str) -> pd.DataFrame:
    """Reads bed file"""
    return pd.read_csv(filepath, sep="\t", engine="c")


def main():
    """d"""
    # Setup bedtools
    track = read_bed(TRACK)
    activity = read_bed(ACTIVITY)

    # Record cols
    track_cols = track.columns.to_list()
    activity_cols = activity.columns.to_list()

    # adjust activity cols to avoid redundancy
    activity_cols = [
        f"activity_{i}" if i in ["chrm", "pos0", "pos1"] else i for i in activity_cols
    ]

    # Make bedtools
    track = pbt.BedTool.from_dataframe(track)
    activity = pbt.BedTool.from_dataframe(activity)

    # Ix
    ix = track.intersect(activity, loj=True).to_dataframe(header=None)
    ix.columns = track_cols + activity_cols
    
    # Records mids that either do not intersect a rDHS
    mids_no_rdhs = ix[ix["activity_chrm"] == "."]["mid"]

    # Drop mids that do not intersect a rDHS
    ix = ix[~ix["mid"].isin(mids_no_rdhs)]

    # Clean up

    # Drop cols
    ix.drop(columns=["activity_chrm", "activity_pos0", "activity_pos1"], inplace=True)

    # Replace "." with "NaN"
    ix = ix.replace(".", "NaN")

    # Round last 2 cols (activity scores)
    ix.iloc[:, -2:] = ix.iloc[:, -2:].astype(float).round(2)

    # Save dropped motifs - no ix
    #mids_no_rdhs.to_csv(DROPPED_OP, sep="\t", index=False)

    # Save matrix
    ix.to_csv(MATRIX_OP, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
