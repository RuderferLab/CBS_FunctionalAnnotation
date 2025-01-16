import pandas as pd

# Snakemake
VARIANTS = snakemake.input[0]  # type: ignore
ATSNP = snakemake.input[1]  # type: ignore
CONSERVATION = snakemake.input[2]  # type: ignore
POSITIONS = snakemake.input[3]  # type: ignore
CBSS = snakemake.input[4]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_variants(filepath: str) -> pd.DataFrame:
    """Returns variants as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        dtype={
            "pid": str,
            "vid": str,
            "af": float,
            "ac": float,
            "an": float,
            "caddPhred": float,
            "caddRaw": float,
            "singleton": float,
            "isvar": int,
            "context": str,
        },
    )


def read_atsnp(filepath: str) -> pd.DataFrame:
    """Returns atsnp as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        dtype={
            "vid": str,
            "dpwm": float,
            "dpwm_log10p": float,
            "dpwm_class": str,
        },
    )


def read_conservation(filepath: str) -> pd.DataFrame:
    """Returns conservation track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        dtype={"pid": str, "gerp++": float, "phastcons100": float, "phylop100": float},
    )


def read_positions(filepath: str) -> pd.DataFrame:
    """Returns positions track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        names=["pid", "mid", "mid_idx"],
        dtype={
            "pid": str,
            "mid": str,
            "mid_idx": int,
        },
    )


def read_cbss(filepath: str) -> pd.DataFrame:
    """Returns CBSs as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        usecols=[
            "score",
            "pval",
            "strand",
            "mid",
            "rDHS",
            "activity",
            "quantile_rdhswide",
            "cid",
            "quantile_cbswide",
        ],
        dtype={
            "score": int,
            "pval": int,
            "strand": str,
            "mid": str,
            "rDHS": str,
            "activity": float,
            "quantile_rdhswide": int,
            "cid": str,
            "quantile_cbswide": int,
        },
    )


def main():
    """Main"""
    print("Combining annotations...")

    # Read in variants
    variants = read_variants(VARIANTS)

    # Read in atsnp track
    atsnp = read_atsnp(ATSNP)

    # Read in conservation track
    conservation = read_conservation(CONSERVATION)

    # Read in positions
    positions = read_positions(POSITIONS)

    # Read in cbss track
    cbss = read_cbss(CBSS)

    ###
    # Merge annotations
    ###

    # Merge atsnp
    variants = variants.merge(atsnp, on="vid", how="left")

    # Merge conservation
    variants = variants.merge(conservation, on="pid", how="left")

    # Merge positions
    variants = variants.merge(positions, on="pid", how="left")

    # Merge cbss
    variants = variants.merge(cbss, on="mid", how="left")

    # Write out
    variants.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
