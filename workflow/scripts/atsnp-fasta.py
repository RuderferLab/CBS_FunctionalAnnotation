"""atsnp-fasta.py"""

import gzip
import pandas as pd
import pyranges as pr


# Snakemake
IP_TRACK = snakemake.input[0]  # type: ignore

# Params
SLOP_SIZE = snakemake.params.slop_size  # type: ignore
HG38_FASTA = snakemake.params.genome_fasta  # type: ignore
OUTDIR = snakemake.params.outdir  # type: ignore
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    return pd.read_csv(
        filepath, sep="\t", engine="c", usecols=["vid"], dtype={"vid": str}
    )


def clean_sequence(seq: str) -> list:
    # NOTE: Update to random, A for now
    random_nuc = "A"
    return seq.upper().replace("N", random_nuc).replace("n", random_nuc)


def mutate_sequence(row: pd.Series) -> str:
    seq = row["refseq_cleaned"]
    alt = row["alt"]
    return seq[:SLOP_SIZE] + alt.upper() + seq[SLOP_SIZE + 1 :]


def main():
    """D"""
    # Read track
    track = read_track(IP_TRACK)

    # Select variants
    track = track[track["vid"] != "."]

    # Expand vid to get coords
    track[["Chromosome", "End", "ref", "alt"]] = track["vid"].str.split(
        "-", expand=True
    )
    track["Chromosome"] = ["chr" + str(i) for i in track["Chromosome"]]
    track["Start"] = track["End"].astype(int) - 1

    # Format to bed
    track = track[["Chromosome", "Start", "End", "vid", "ref", "alt"]]

    # Make pr and extend intervals
    pr_track = pr.PyRanges(track).extend(SLOP_SIZE)

    # Add sequence
    sequences = pr.get_fasta(pr_track, HG38_FASTA)

    # Convert to df and add seq column
    track = pr_track.as_df()
    track["refseq"] = sequences

    # Clean sequences
    track["refseq_cleaned"] = track["refseq"].apply(lambda x: clean_sequence(x))

    # Make alt sequence
    track["altseq_cleaned"] = track.apply(mutate_sequence, axis=1)

    # Write ref and alt sequences as fasta using vid as name
    for chromosome in CHROMOSOMES:
        with gzip.open(f"{OUTDIR}/ref/{chromosome}.fasta.gz", "wb") as ref, gzip.open(
            f"{OUTDIR}/alt/{chromosome}.fasta.gz", "wb"
        ) as alt:
            for i, row in track[track["Chromosome"] == chromosome].iterrows():
                ref.write(f">{row['vid']}\n{row['refseq_cleaned']}\n".encode('utf-8'))
                alt.write(f">{row['vid']}\n{row['altseq_cleaned']}\n".encode('utf-8'))


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
