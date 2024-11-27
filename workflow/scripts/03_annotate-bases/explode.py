import pandas as pd

# Snakemake
IP_TRACK = snakemake.input[0]  # type: ignore
OP_TRACK = snakemake.output[0]  # type: ignore

# Params
MOTIF_LENGTH = snakemake.params["motif_length"]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    """Returns track as DataFrame"""
    fields = [0, 1, 2, 7]
    dtypes = {"chrm": str, "pos0": int, "pos1": int, "mid": str}
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        usecols=fields,
        names=dtypes.keys(),
        dtype=dtypes,
    )


def main():
    """d"""
    # Read track
    track = read_track(IP_TRACK)

    # Get beginning of 1-based positions in motif
    track["motif_start"] = (track["pos1"] - MOTIF_LENGTH) + 1

    # Get range of nucleotides in motifs
    track["range"] = track.apply(
        lambda row: range(row["motif_start"], row["pos1"] + 1), axis=1
    )

    # Get motif index
    track["idx"] = track["range"].apply(
        lambda x: [idx for idx, _ in enumerate(x, start=1)]
    )

    # Explode range
    track = track.explode(["range", "idx"]).rename(columns={"range": "pos_end"})
    track.reset_index(drop=True, inplace=True)
    track.insert(2, "pos_start", track["pos_end"] - 1)
    track.drop(columns=["pos0", "pos1"], inplace=True)

    # Rename columns
    track.rename(columns={"pos_start": "pos0", "pos_end": "pos1"}, inplace=True)

    # Add position id
    track.insert(3, "pid", track["chrm"] + "-" + track["pos1"].astype(str))

    # Save
    fields = ["chrm", "pos0", "pos1", "pid", "mid", "idx"]
    track[fields].to_csv(
        OP_TRACK,
        sep="\t",
        header=False,
        index=False,
    )


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
