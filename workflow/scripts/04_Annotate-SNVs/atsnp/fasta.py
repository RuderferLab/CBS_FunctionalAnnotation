import pandas as pd
import pybedtools as pbt  # type: ignore
from Bio import SeqIO
import gzip
from collections import defaultdict


# Snakemake
IP_TRACK = snakemake.input[0]  # type: ignore
OUTPUT_REF = snakemake.output[0]  # type: ignore
OUTPUT_ALT = snakemake.output[1]  # type: ignore

# Params
CHROMOSOME = snakemake.params.chromosome  # type: ignore
SLOP_SIZE = snakemake.params.slop_size  # type: ignore
GENOME_SIZES = snakemake.params.genome_sizes  # type: ignore
GENOME_FASTA = snakemake.params.genome_fasta  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_track(filepath: str) -> pd.DataFrame:
    return pd.read_csv(filepath, sep="\t", engine="c", usecols=["vid"])


def main():
    """D"""

    # Read track
    track = read_track(IP_TRACK)

    # Drop positions with no var
    track = track.dropna(subset=["vid"])

    # Expand vid to get coords
    track[["chrm", "pos1", "ref", "alt"]] = track["vid"].str.split("-", expand=True)
    track["chrm"] = ["chr" + str(i) for i in track["chrm"]]
    track["pos0"] = track["pos1"].astype(int) - 1

    # Format to bed
    track = track[["chrm", "pos0", "pos1", "vid", "ref", "alt"]]

    # Filter to chromosome
    track = track[track["chrm"] == CHROMOSOME]

    # Make pbt
    track = pbt.BedTool.from_dataframe(track)

    # Slop and return to df
    track = track.slop(b=SLOP_SIZE, g=GENOME_SIZES).sort()

    # Get seq
    sequence = track.sequence(fi=GENOME_FASTA, nameOnly=True)

    print("D")

    # Store
    ref_file = []
    alt_file = []

    # Read in seq
    for record in SeqIO.parse(sequence.seqfn, "fasta"):
        # Mutate sequence based on slop size and nuc
        info = record.id.split("-")
        seq = str(record.seq)
        # Check if "n" or "N" in seq and if so replace with random ["A",  C", "G", "T"]
        if "n" in seq or "N" in seq:
            print("FOUND N IN SEQ...")
            for idx, nuc in enumerate(seq):
                if nuc == "n" or nuc == "N":
                    # rand_nuc = str(np.random.choice(["A", "C", "G", "T"]))
                    # Just gonna put A in there for now - need to update to store the rand seq for alt allele that comes after
                    rand_nuc = "A"
                    seq = seq[:idx] + rand_nuc + seq[idx + 1 :]

        # dont have to mutate for ref
        ref_file.append([record.id, seq.upper()])

        # Now handle alt
        # Mutate seq
        nuc = info[3]
        # Mutate using slop size as index
        seq = seq[:SLOP_SIZE] + nuc.upper() + seq[SLOP_SIZE + 1 :]
        alt_file.append([record.id, seq.upper()])

    # Write out

    # REF
    with gzip.open(OUTPUT_REF, "at") as f:
        for record in ref_file:
            f.write(f">{record[0]}\n{record[1]}\n")

    # ALT
    with gzip.open(OUTPUT_ALT, "at") as f:
        for record in alt_file:
            f.write(f">{record[0]}\n{record[1]}\n")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
