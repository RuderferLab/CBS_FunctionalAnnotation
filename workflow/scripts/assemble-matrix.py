"""assemble-matrix.py"""

import pandas as pd

# Snakemake
VARIANTS = snakemake.input[0]  # type: ignore
ATSNP = snakemake.input[1]  # type: ignore
CONSERVATION = snakemake.input[2]  # type: ignore
MOTIFS = snakemake.input[3]  # type: ignore
CBSS = snakemake.input[4]  # type: ignore
CBS_POSITIONS_MAP = snakemake.input[5]  # type: ignore
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
        compression="gzip",
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
        compression="gzip",
    )


def read_conservation(filepath: str) -> pd.DataFrame:
    """Returns conservation track as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        # Note order of conservation scores
        dtype={"pid": str, "phastcons": str, "phylop": str, "gerp++": str, "linsight": str},
        header=None,
        names=["pid", "phastcons", "phylop", "gerp++", "linsight"],
        compression="gzip",
    )

def read_motifs(filepath: str) -> pd.DataFrame:
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
        ],
        dtype={
            "score": int,
            "pval": int,
            "strand": str,
            "mid": str,
        },
        compression="gzip",
    )

def read_cbss(filepath: str) -> pd.DataFrame:
    """Returns CBSs as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        usecols=[
            "mid",
            "rDHS",
            "activity",
            "quantile_rdhswide",
            "is_cbs",
            "cid"
        ],
        dtype={
            "mid": str,
            "rDHS": str,
            "activity": float,
            "quantile_rdhswide": int,
            "is_cbs": int,
            "cid": str,
        },
    )
    
def read_position_map(filepath: str) -> pd.DataFrame:
    """Returns CBS positions to motif as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        engine="c",
        header=None,
        dtype={
            "pid": str,
            "mid": str,
        },
        names=["pid", "mid"],
    )

def main():
    """Main"""
    print("Combining annotations...")

    # Read in variants
    variants = read_variants(VARIANTS)

    # Read in atsnp track
    atsnp = read_atsnp(ATSNP)

    # Read in conservation track - mask "." with NaN
    conservation = read_conservation(CONSERVATION)
    conservation = conservation.replace(".", np.nan)

    # Read in motifs
    motifs = read_motifs(MOTIFS)

    # Read in cbss track
    cbss = read_cbss(CBSS)

    # Read in cbs position to mid map
    cbs_positions_map = read_position_map(CBS_POSITIONS_MAP)

    ###
    # Merge annotations
    ###

    # Merge atsnp
    variants = variants.merge(atsnp, on="vid", how="left")

    # Merge conservation
    variants = variants.merge(conservation, on="pid", how="left")

    # Merge positions to get motif link
    variants = variants.merge(cbs_positions_map, on="pid", how="left")

    # Merge with motifs to get motif information
    variants = variants.merge(motifs, on="mid", how="left")

    # Merge with cbss to get activity and rDHS
    variants = variants.merge(cbss, on="mid", how="left")

    # Write out
    variants.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
