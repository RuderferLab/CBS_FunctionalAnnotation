IP_DIR="../../resources/conservation"
OP_DIR="resources/preprocessed/conservation"
CTCF="resources/preprocessed/CTCF-motifs.MA01392.hg38.tsv.processed.gz"

FILES="gerp.hg19.bed.gz linsight.hg19.bed.gz phastcons100way.hg38.bed.gz phylop100way.hg38.bed.gz"

# Loop through all the files
for FILE in $FILES
do
    echo "Processing $FILE"
    basename=$(echo $FILE | cut -d'.' -f1-2)
    # Intersect with CTCF motifs
    zcat $CTCF |
    vawk '{if (NR>1) print $1, $2, $3}' |
    bedtools intersect -sorted -a ${IP_DIR}/${FILE} -b stdin | 
    starch --bzip2 - > ${OP_DIR}/CTCF-motifs.MA01392.positions.${basename}.starch.bzip
done