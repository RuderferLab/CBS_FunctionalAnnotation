library(atSNP)

# Snakemake
JASPAR_PFM <- snakemake@input[[1]]
REF_FASTA <- snakemake@input[[2]]
ALT_FASTA <- snakemake@input[[3]]
OUTPUT <- snakemake@output[[1]]

# ------------- #
# Main          #
# ------------- #

print("Running...")

# Load PFM with pseudocounts
pwm <- LoadMotifLibrary(JASPAR_PFM, tag = ">", skiprows = 1, skipcols = 1, transpose = TRUE, 
 field = 1, sep = c("\t", " ", "\\[", "\\]", ">"),
 pseudocount = 1)

print("Read in PWM.")

# Read in snp - NOTE DEFAULT PAR
snp_info <- LoadFastaData(REF_FASTA, ALT_FASTA, default.par = FALSE)

print("Done reading in snp info.")

# Calculate scores
atsnp.scores <- ComputeMotifScore(pwm, snp_info, ncores = 8)

print("Done calculating scores..")

# Pvalues and results #NOTE: TESTTING.MC=FALSE
atsnp.result <- ComputePValues(motif.lib = pwm, snp.info = snp_info,
                    motif.scores = atsnp.scores$motif.scores, ncores = 8, testing.mc=FALSE)

print("Done calculating pvalues..")

# Subset results
results <- atsnp.result[, c(
  "snpid", "motif", "pval_ref", "pval_snp", "pval_diff", 
  "pval_cond_ref", "pval_cond_snp", "pval_rank", 
  "snpbase", "log_lik_ref", "log_lik_snp", "log_lik_ratio"
)]


print('writing out...')

# Round matrix to 4 decimal places
results[,c("pval_ref", "pval_snp", "pval_diff", "pval_cond_ref", "pval_cond_snp", "pval_rank", "log_lik_ref", "log_lik_snp", "log_lik_ratio")] <- round(results[,c("pval_ref", "pval_snp", "pval_diff", "pval_cond_ref", "pval_cond_snp", "pval_rank", "log_lik_ref", "log_lik_snp", "log_lik_ratio")], 4)

# Write results - tab  sep and compressed
write.table(results, file = gzfile(OUTPUT),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)