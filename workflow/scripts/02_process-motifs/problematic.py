import pandas as pd
import pybedtools as pbt


# Snakemake
TRACK = snakemake.input[0]  # type:ignore
MASK = snakemake.input[1]  # type:ignore
OUTPUT = snakemake.output[0]  # type:ignore

# ------------- #
# Functions     #
# ------------- #


def read_annotation(filepath: str) -> pd.DataFrame:
    """Reads annotation"""
    return pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrm", "start", "end"],
        engine="c",
    )


def main():
    ###
    # Main
    ###

    # Make bedtools
    track_bt = pbt.BedTool(TRACK).sort()
    mask_bt = pbt.BedTool(MASK).sort()

    # Flag overlaps
    ix = track_bt.intersect(b=mask_bt, c=True)

    # Cluster to assess tiled
    ix = ix.cluster().to_dataframe(header=None)

    # Set columns
    ix.columns = [
        "chrm",
        "pos0",
        "pos1",
        "tf_name",
        "score",
        "pval",
        "strand",
        "exclude_count",
        "clus_id",
    ]

    # Add motif id (mid)
    ix["mid"] = (
        ix["chrm"]
        + ":"
        + ix["pos0"].astype(str)
        + ":"
        + ix["pos1"].astype(str)
        + ":"
        + ix["strand"]
        + ":"
        + ix["tf_name"]
    )

    ###
    # Filtering
    ###

    # Flag overlaps with exclude files
    ix = ix[ix["exclude_count"] == 0]

    # Flag non-primary contigs
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    ix = ix[ix["chrm"].isin(chroms)]

    # Drop tiled - select highest
    ix = ix.sort_values(["clus_id", "score"], ascending=[True, False]).drop_duplicates(
        "clus_id", keep="first"
    )

    # Drop cols
    ix.drop(columns=["exclude_count", "clus_id"], inplace=True)

    # Save, no header for processing later
    ix.to_csv(OUTPUT, sep="\t", index=False, header=None)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
