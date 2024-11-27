import pandas as pd

# Snakemake
TRACK = snakemake.input[0]  # type: ignore
ACTIVITY = snakemake.input[1]  # type: ignore
CONSERVATION = snakemake.input[2]  # type: ignore
VARIANTS = snakemake.input[3]  # type: ignore
ATSNP = snakemake.input[4]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_base_track(filepath: str) -> pd.DataFrame:
    """Returns track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        usecols=[3, 4],
        names=["pid", "mid"],
        dtype={"pid": str, "mid": str},
    )
    
def read_activity(filepath: str) -> pd.DataFrame:
    """Returns activity track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        usecols=[4, 5, 6, 7, 8, 9, 10],
        names=["pwm_score", "pwm_pval", "pwm_strand", "mid", "rdhs", "activity", "activity_masked"],
    )
    
def read_annotation(filepath: str) -> pd.DataFrame:
    """Returns track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
    )


def main():
    """d"""

    print("Combining annotations...")

    # Read in positions track
    track = read_base_track(TRACK)

    # Read in annotations
    activity = read_activity(ACTIVITY)
    variants = read_annotation(VARIANTS)
    conservation = read_annotation(CONSERVATION)
    atsnp = read_annotation(ATSNP)

    # Start with variants as left - this is all possible observations due to multiallelic
    print("Merging variants...")
    track = pd.merge(variants, track, on=["pid", "mid"], how="left")

    # Merge with activity
    print("Merging activity...")
    track = pd.merge(track, activity, on=["mid"], how="left")

    # Merge with conservation
    print("Merging conservation...")
    conservation = conservation.round(6)
    track = pd.merge(track, conservation, on=["pid", "mid"], how="left")

    # Merge with atsnp
    print("Merging atsnp...")
    track = pd.merge(track, atsnp, on=["vid"], how="left")

    ###
    # Clean up
    ###

    # Fill in missing values
    track.fillna("NaN", inplace=True)

    # Write out
    track.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
