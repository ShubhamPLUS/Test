## Script to generate all synthetic test datasets for ProteoForge.
## Run from the package root: Rscript data-raw/generate_synthetic_datasets.R
## Outputs go to inst/extdata/<dataset_name>/

set.seed(1234)

library(data.table)
library(fs)

out_root <- path("inst", "extdata")
dir_create(out_root)

# ── Helpers ──────────────────────────────────────────────────────────────────

.protein_ids <- function(n) sprintf("P%05d", seq_len(n))
.gene_syms   <- function(n) paste0("GENE", seq_len(n))

.sim_lfq <- function(n_proteins, n_samples, n_up, n_down,
                     lfc_true = 2, miss_frac = 0.1,
                     seed = 1234) {
  set.seed(seed)
  base_log2 <- matrix(
    rnorm(n_proteins * n_samples, mean = 22, sd = 2),
    nrow = n_proteins, ncol = n_samples
  )
  rownames(base_log2) <- .protein_ids(n_proteins)
  sample_names <- paste0("S", sprintf("%02d", seq_len(n_samples)))
  colnames(base_log2) <- sample_names

  half <- n_samples %/% 2
  grp  <- c(rep("Control", half), rep("Treatment", n_samples - half))

  # Spike in true DE
  up_idx   <- seq_len(n_up)
  down_idx <- seq(n_up + 1, n_up + n_down)
  for (j in which(grp == "Treatment")) {
    base_log2[up_idx,   j] <- base_log2[up_idx,   j] + lfc_true
    base_log2[down_idx, j] <- base_log2[down_idx, j] - lfc_true
  }

  # Introduce missingness
  n_miss <- round(n_proteins * n_samples * miss_frac)
  miss_idx <- sample(length(base_log2), n_miss)
  base_log2[miss_idx] <- NA

  list(
    matrix   = base_log2,
    group    = grp,
    up_idx   = up_idx,
    down_idx = down_idx,
    samples  = sample_names
  )
}

# ── 1. generic_lfq_small ─────────────────────────────────────────────────────
# 300 proteins, 12 samples (6 Control, 6 Treatment), 30 up + 30 down

cat("Generating generic_lfq_small ...\n")
ds1 <- .sim_lfq(300, 12, 30, 30)

mat1 <- 2^ds1$matrix  # back to linear scale for the "raw" file

quant_dt <- data.table(
  ProteinID   = rownames(ds1$matrix),
  GeneSymbol  = .gene_syms(300),
  Description = paste("Protein", .protein_ids(300)),
  as.data.table(mat1)
)

smd1 <- data.table(
  sample_id   = ds1$samples,
  raw_file    = paste0(ds1$samples, "_DIA.raw"),
  condition   = ds1$group,
  replicate   = c(1,2,3,4,5,6,1,2,3,4,5,6),
  batch       = rep("batch1", 12),
  instrument  = "Orbitrap_Astral"
)

truth1 <- data.table(
  ProteinID = rownames(ds1$matrix),
  true_log2fc_Treatment_vs_Control = ifelse(
    seq_len(300) %in% ds1$up_idx, 2,
    ifelse(seq_len(300) %in% ds1$down_idx, -2, 0)
  ),
  true_de = seq_len(300) %in% c(ds1$up_idx, ds1$down_idx)
)

dir1 <- path(out_root, "generic_lfq_small")
dir_create(dir1)
fwrite(quant_dt, path(dir1, "protein_matrix.tsv"), sep = "\t")
fwrite(smd1,     path(dir1, "sample_metadata.tsv"), sep = "\t")
fwrite(truth1,   path(dir1, "ground_truth.tsv"),    sep = "\t")
cat("  Wrote", dir1, "\n")

# ── 2. generic_lfq_batch ─────────────────────────────────────────────────────
# 2 conditions × 2 batches, controlled batch effect

cat("Generating generic_lfq_batch ...\n")
set.seed(42)
n_prot2 <- 200; n_samp2 <- 16
base2 <- matrix(rnorm(n_prot2 * n_samp2, 22, 2), n_prot2, n_samp2)
rownames(base2) <- .protein_ids(n_prot2)
s2_names <- paste0("S", sprintf("%02d", seq_len(n_samp2)))
colnames(base2) <- s2_names

grp2   <- rep(c("Control","Treatment"), each = n_samp2/2)
batch2 <- rep(c("batch1","batch2"), times = n_samp2/2)

# True DE: first 20 proteins up in Treatment
for (j in which(grp2 == "Treatment")) base2[1:20, j] <- base2[1:20, j] + 2

# Batch effect: batch2 samples are globally shifted +1.5 on select proteins
for (j in which(batch2 == "batch2")) base2[21:80, j] <- base2[21:80, j] + 1.5

mat2 <- 2^base2
quant2 <- data.table(
  ProteinID   = rownames(base2),
  GeneSymbol  = .gene_syms(n_prot2),
  Description = paste("Protein", .protein_ids(n_prot2)),
  as.data.table(mat2)
)
smd2 <- data.table(
  sample_id   = s2_names,
  raw_file    = paste0(s2_names, "_DIA.raw"),
  condition   = grp2,
  replicate   = rep(1:4, times=4),
  batch       = batch2,
  instrument  = "Orbitrap_Eclipse"
)

dir2 <- path(out_root, "generic_lfq_batch")
dir_create(dir2)
fwrite(quant2, path(dir2, "protein_matrix.tsv"), sep = "\t")
fwrite(smd2,   path(dir2, "sample_metadata.tsv"), sep = "\t")
cat("  Wrote", dir2, "\n")

# ── 3. paired_design_small ────────────────────────────────────────────────────
# 6 paired subjects, Before/After treatment

cat("Generating paired_design_small ...\n")
set.seed(99)
n_prot3 <- 250; n_samp3 <- 12
base3 <- matrix(rnorm(n_prot3 * n_samp3, 22, 2), n_prot3, n_samp3)
rownames(base3) <- .protein_ids(n_prot3)
s3_names <- c(paste0("Before_P", 1:6), paste0("After_P", 1:6))
colnames(base3) <- s3_names

# Subject-level random effects
for (i in 1:6) {
  subj_eff <- rnorm(1, 0, 1)
  base3[, i]     <- base3[, i]     + subj_eff
  base3[, i + 6] <- base3[, i + 6] + subj_eff
}
# True DE: first 25 proteins go up after treatment
for (j in 7:12) base3[1:25, j] <- base3[1:25, j] + 2

mat3 <- 2^base3
quant3 <- data.table(
  ProteinID   = rownames(base3),
  GeneSymbol  = .gene_syms(n_prot3),
  Description = paste("Protein", .protein_ids(n_prot3)),
  as.data.table(mat3)
)
smd3 <- data.table(
  sample_id   = s3_names,
  raw_file    = paste0(s3_names, "_DIA.raw"),
  condition   = c(rep("Before", 6), rep("After", 6)),
  replicate   = c(1:6, 1:6),
  subject_id  = c(paste0("P", 1:6), paste0("P", 1:6)),
  batch       = "batch1",
  instrument  = "Orbitrap_Astral"
)

dir3 <- path(out_root, "paired_design_small")
dir_create(dir3)
fwrite(quant3, path(dir3, "protein_matrix.tsv"), sep = "\t")
fwrite(smd3,   path(dir3, "sample_metadata.tsv"), sep = "\t")
cat("  Wrote", dir3, "\n")

# ── 4. DIA-NN-like example ─────────────────────────────────────────────────
# Minimal pg_matrix and report.tsv schemas

cat("Generating diann_example ...\n")
set.seed(77)
n_pg <- 150; n_samp4 <- 8
samples4 <- paste0("Run", seq_len(n_samp4))
grp4 <- c(rep("Ctrl",4), rep("Treat",4))

pg_mat <- data.table(
  `Protein.Group`    = .protein_ids(n_pg),
  `Protein.Ids`      = .protein_ids(n_pg),
  `Protein.Names`    = .gene_syms(n_pg),
  `Genes`            = .gene_syms(n_pg),
  `First.Protein.Description` = paste("Description", seq_len(n_pg))
)
for (s in samples4) {
  pg_mat[[s]] <- 2^rnorm(n_pg, 22, 2)
}

# Minimal report.tsv (precursor-level)
n_prec <- n_pg * 5
report_dt <- data.table(
  `Run`                    = sample(samples4, n_prec, replace = TRUE),
  `Protein.Group`          = sample(.protein_ids(n_pg), n_prec, replace = TRUE),
  `Protein.Ids`            = sample(.protein_ids(n_pg), n_prec, replace = TRUE),
  `Protein.Names`          = sample(.gene_syms(n_pg), n_prec, replace = TRUE),
  `Genes`                  = sample(.gene_syms(n_pg), n_prec, replace = TRUE),
  `Precursor.Id`           = paste0("PEP", sprintf("%05d", seq_len(n_prec))),
  `Modified.Sequence`      = paste0("PEPT", seq_len(n_prec)),
  `Stripped.Sequence`      = paste0("PEPT", seq_len(n_prec)),
  `Precursor.Charge`       = sample(2:4, n_prec, replace = TRUE),
  `Precursor.Normalised`   = 2^rnorm(n_prec, 20, 2),
  `PG.MaxLFQ`              = 2^rnorm(n_prec, 22, 2),
  `Global.Q.Value`         = runif(n_prec, 0, 0.02),
  `PG.Q.Value`             = runif(n_prec, 0, 0.02),
  `Lib.Q.Value`            = runif(n_prec, 0, 0.02),
  `Lib.PG.Q.Value`         = runif(n_prec, 0, 0.02)
)

smd4 <- data.table(
  sample_id  = samples4,
  raw_file   = paste0(samples4, ".raw"),
  condition  = grp4,
  replicate  = c(1,2,3,4,1,2,3,4),
  batch      = "batch1"
)

dir4 <- path(out_root, "diann_example")
dir_create(dir4)
fwrite(pg_mat,    path(dir4, "report.pg_matrix.tsv"), sep = "\t")
fwrite(report_dt, path(dir4, "report.tsv"),           sep = "\t")
fwrite(smd4,      path(dir4, "sample_metadata.tsv"),  sep = "\t")
cat("  Wrote", dir4, "\n")

# ── 5. Spectronaut-like example ────────────────────────────────────────────
cat("Generating spectronaut_example ...\n")
set.seed(88)
n_pg5 <- 120; n_samp5 <- 8
samples5 <- paste0("File", seq_len(n_samp5))
grp5 <- c(rep("Ctrl",4), rep("Treat",4))
n_rows5 <- n_pg5 * 10

sn_report <- data.table(
  `R.FileName`              = sample(samples5, n_rows5, replace = TRUE),
  `PG.ProteinGroups`        = sample(.protein_ids(n_pg5), n_rows5, replace = TRUE),
  `PG.Genes`                = sample(.gene_syms(n_pg5), n_rows5, replace = TRUE),
  `EG.IsDecoy`              = FALSE,
  `EG.Qvalue`               = runif(n_rows5, 0, 0.02),
  `PG.Qvalue`               = runif(n_rows5, 0, 0.02),
  `PG.Quantity`             = 2^rnorm(n_rows5, 22, 2),
  `PEP.Quantity`            = 2^rnorm(n_rows5, 20, 2),
  `FG.MS2Quantity`          = 2^rnorm(n_rows5, 18, 2),
  `EG.ModifiedSequence`     = paste0("[Prot]PEPT", seq_len(n_rows5)),
  `PEP.StrippedSequence`    = paste0("PEPT", seq_len(n_rows5)),
  `FG.Charge`               = sample(2:4, n_rows5, replace = TRUE)
)

smd5 <- data.table(
  sample_id  = samples5,
  raw_file   = paste0(samples5, ".raw"),
  condition  = grp5,
  replicate  = c(1,2,3,4,1,2,3,4),
  batch      = "batch1"
)

dir5 <- path(out_root, "spectronaut_example")
dir_create(dir5)
fwrite(sn_report, path(dir5, "spectronaut_report.tsv"), sep = "\t")
fwrite(smd5,      path(dir5, "sample_metadata.tsv"),    sep = "\t")
cat("  Wrote", dir5, "\n")

# ── 6. FragPipe-like example ───────────────────────────────────────────────
cat("Generating fragpipe_example ...\n")
set.seed(55)
n_prot6 <- 130; n_samp6 <- 8
samples6 <- paste0("Sample", seq_len(n_samp6))
grp6 <- c(rep("Ctrl",4), rep("Treat",4))

fp_mat <- data.table(
  `Protein`              = .protein_ids(n_prot6),
  `Protein ID`           = .protein_ids(n_prot6),
  `Entry Name`           = .gene_syms(n_prot6),
  `Gene`                 = .gene_syms(n_prot6),
  `Protein Description`  = paste("FragPipe protein", seq_len(n_prot6)),
  `Protein Probability`  = runif(n_prot6, 0.99, 1),
  `Top Peptide Probability` = runif(n_prot6, 0.9, 1),
  `Unique Spectral Count`   = sample(5:50, n_prot6, replace = TRUE),
  `Razor Spectral Count`    = sample(5:50, n_prot6, replace = TRUE),
  `Total Spectral Count`    = sample(10:100, n_prot6, replace = TRUE),
  `Unique Peptide Count`    = sample(2:20, n_prot6, replace = TRUE),
  `Razor Peptide Count`     = sample(2:20, n_prot6, replace = TRUE),
  `Total Peptide Count`     = sample(5:30, n_prot6, replace = TRUE)
)
for (s in samples6) {
  fp_mat[[paste0(s, " MaxLFQ Intensity")]] <- 2^rnorm(n_prot6, 22, 2)
}

smd6 <- data.table(
  sample_id  = samples6,
  raw_file   = paste0(samples6, ".raw"),
  condition  = grp6,
  replicate  = c(1,2,3,4,1,2,3,4),
  batch      = "batch1"
)

dir6 <- path(out_root, "fragpipe_example")
dir_create(dir6)
fwrite(fp_mat, path(dir6, "combined_protein.tsv"), sep = "\t")
fwrite(smd6,   path(dir6, "sample_metadata.tsv"),  sep = "\t")
cat("  Wrote", dir6, "\n")

# ── 7. MaxQuant-like example ───────────────────────────────────────────────
cat("Generating maxquant_example ...\n")
set.seed(33)
n_prot7 <- 140; n_samp7 <- 8
samples7 <- paste0("Experiment", seq_len(n_samp7))
grp7 <- c(rep("Ctrl",4), rep("Treat",4))

mq_mat <- data.table(
  `Protein IDs`         = .protein_ids(n_prot7),
  `Majority protein IDs`= .protein_ids(n_prot7),
  `Gene names`          = .gene_syms(n_prot7),
  `Protein names`       = paste("MaxQuant protein", seq_len(n_prot7)),
  `Number of proteins`  = sample(1:4, n_prot7, replace = TRUE),
  `Peptides`            = sample(2:30, n_prot7, replace = TRUE),
  `Razor + unique peptides` = sample(2:25, n_prot7, replace = TRUE),
  `Unique peptides`     = sample(2:20, n_prot7, replace = TRUE),
  `Mol. weight [kDa]`   = runif(n_prot7, 10, 200),
  `Q-value`             = runif(n_prot7, 0, 0.01),
  `Score`               = runif(n_prot7, 50, 300),
  `Reverse`             = c(rep("", n_prot7 - 5), rep("+", 5)),
  `Potential contaminant` = c(rep("", n_prot7 - 3), rep("+", 3)),
  `Only identified by site` = c(rep("", n_prot7 - 2), rep("+", 2))
)
for (s in samples7) {
  mq_mat[[paste0("LFQ intensity ", s)]] <- 2^rnorm(n_prot7, 22, 2)
  mq_mat[[paste0("iBAQ ", s)]]          <- 2^rnorm(n_prot7, 20, 2)
  mq_mat[[paste0("Intensity ", s)]]     <- 2^rnorm(n_prot7, 23, 2)
}

smd7 <- data.table(
  sample_id   = samples7,
  raw_file    = paste0(samples7, ".raw"),
  condition   = grp7,
  replicate   = c(1,2,3,4,1,2,3,4),
  batch       = "batch1",
  instrument  = "Orbitrap_Eclipse"
)

dir7 <- path(out_root, "maxquant_example")
dir_create(dir7)
fwrite(mq_mat, path(dir7, "proteinGroups.txt"), sep = "\t")
fwrite(smd7,   path(dir7, "sample_metadata.tsv"), sep = "\t")
cat("  Wrote", dir7, "\n")

# ── 8. PD-like example ────────────────────────────────────────────────────
cat("Generating pd_example ...\n")
set.seed(21)
n_prot8 <- 110; n_samp8 <- 8
samples8 <- paste0("F", seq_len(n_samp8))
grp8 <- c(rep("Ctrl",4), rep("Treat",4))

pd_mat <- data.table(
  `Accession`            = .protein_ids(n_prot8),
  `Description`          = paste("PD protein", seq_len(n_prot8)),
  `Gene Symbol`          = .gene_syms(n_prot8),
  `Master`               = c(rep("IsMasterProtein", n_prot8 - 5),
                              rep("", 5)),
  `Protein FDR Confidence: Combined` = c(rep("High", n_prot8 - 8),
                                          rep("Medium", 5), rep("Low", 3)),
  `Contaminant`          = c(rep(FALSE, n_prot8 - 3), rep(TRUE, 3)),
  `# Peptides`           = sample(2:25, n_prot8, replace = TRUE),
  `# Unique Peptides`    = sample(2:20, n_prot8, replace = TRUE),
  `Coverage [%]`         = runif(n_prot8, 5, 80)
)
for (i in seq_len(n_samp8)) {
  pd_mat[[sprintf("Abundance: F%d: Sample, %s", i, grp8[i])]] <-
    2^rnorm(n_prot8, 22, 2)
}

smd8 <- data.table(
  sample_id   = samples8,
  raw_file    = paste0("F", seq_len(n_samp8), "_raw"),
  condition   = grp8,
  replicate   = c(1,2,3,4,1,2,3,4),
  batch       = "batch1"
)

dir8 <- path(out_root, "pd_example")
dir_create(dir8)
fwrite(pd_mat, path(dir8, "proteins.tsv"), sep = "\t")
fwrite(smd8,   path(dir8, "sample_metadata.tsv"), sep = "\t")
cat("  Wrote", dir8, "\n")

# ── 9. Phospho-LFQ small ──────────────────────────────────────────────────
cat("Generating phospho_lfq_small ...\n")
set.seed(66)
n_sites <- 500; n_samp9 <- 12
samples9 <- paste0("S", sprintf("%02d", seq_len(n_samp9)))
grp9 <- c(rep("Ctrl",6), rep("Treat",6))
parent_prots9 <- sample(.protein_ids(100), n_sites, replace = TRUE)

phospho_mat <- data.table(
  SiteID            = paste0(parent_prots9, "_",
                             sample(c("S","T","Y"), n_sites, replace=TRUE),
                             sample(100:900, n_sites, replace=FALSE)),
  ProteinID         = parent_prots9,
  GeneSymbol        = paste0("GENE", match(parent_prots9, unique(parent_prots9))),
  Residue           = sample(c("S","T","Y"), n_sites, replace = TRUE),
  Position          = sample(100:900, n_sites, replace = TRUE),
  LocalizationProb  = runif(n_sites, 0.5, 1),
  Multiplicity      = sample(1:3, n_sites, replace = TRUE)
)
# Add quant columns
for (s in samples9) {
  phospho_mat[[s]] <- 2^rnorm(n_sites, 20, 2)
}
# True DE: first 50 sites up, next 50 down
for (j in which(grp9 == "Treat")) {
  phospho_mat[1:50,  (samples9[j]) := phospho_mat[[samples9[j]]][1:50] * 4]
  phospho_mat[51:100,(samples9[j]) := phospho_mat[[samples9[j]]][51:100] / 4]
}

# Paired proteome subset
prot_parent <- data.table(
  ProteinID  = unique(parent_prots9),
  GeneSymbol = paste0("GENE", seq_along(unique(parent_prots9)))
)
for (s in samples9) {
  prot_parent[[s]] <- 2^rnorm(nrow(prot_parent), 22, 2)
}

phospho_smd <- data.table(
  sample_id  = samples9,
  raw_file   = paste0(samples9, ".raw"),
  condition  = grp9,
  replicate  = c(1:6, 1:6),
  batch      = "batch1"
)

dir9 <- path(out_root, "phospho_lfq_small")
dir_create(dir9)
fwrite(phospho_mat,  path(dir9, "phospho_sites.tsv"),    sep = "\t")
fwrite(prot_parent,  path(dir9, "protein_matrix.tsv"),   sep = "\t")
fwrite(phospho_smd,  path(dir9, "sample_metadata.tsv"),  sep = "\t")
cat("  Wrote", dir9, "\n")

cat("\nAll synthetic datasets generated successfully.\n")
