#' Import MaxQuant proteinGroups.txt
#'
#' Reads `proteinGroups.txt`, removes reverse hits, contaminants, and
#' only-identified-by-site entries, extracts LFQ intensity columns, and
#' (for phospho) optionally processes `Phospho (STY)Sites.txt` with proper
#' multiplicity expansion.
#'
#' @param quant_file    Path to `proteinGroups.txt`.
#' @param metadata_file Path to sample metadata.
#' @param config        ProteoForge config list.
#' @param sites_file    Optional path to `Phospho (STY)Sites.txt` for
#'   phosphosite-level import.
#'
#' @return A `ProteoForgeData` object.
#'
#' @examples
#' qf <- system.file("extdata","maxquant_example","proteinGroups.txt",
#'                   package="proteoforge")
#' mf <- system.file("extdata","maxquant_example","sample_metadata.tsv",
#'                   package="proteoforge")
#' pfd <- import_maxquant(qf, mf)
#'
#' @export
import_maxquant <- function(quant_file, metadata_file, config = list(),
                             sites_file = NULL) {
  stopifnot(file.exists(quant_file), file.exists(metadata_file))
  logger::log_info("import_maxquant: {quant_file}")

  smd <- .read_metadata(metadata_file)
  data.table::setDT(smd)
  cfg <- config$filtering %||% list()

  dt <- .read_delim_auto(quant_file)

  # ── Hard filters ───────────────────────────────────────────────────────────
  if (isTRUE(cfg$remove_reverse_decoys %||% TRUE) &&
      "Reverse" %in% names(dt)) {
    dt <- dt[Reverse != "+" | is.na(Reverse)]
  }
  if (isTRUE(cfg$remove_contaminants %||% TRUE) &&
      "Potential contaminant" %in% names(dt)) {
    dt <- dt[`Potential contaminant` != "+" | is.na(`Potential contaminant`)]
  }
  if (isTRUE(cfg$remove_only_identified_by_site %||% TRUE) &&
      "Only identified by site" %in% names(dt)) {
    dt <- dt[`Only identified by site` != "+" | is.na(`Only identified by site`)]
  }

  # ── LFQ intensity columns ──────────────────────────────────────────────────
  lfq_cols    <- grep("^LFQ intensity ", names(dt), value=TRUE)
  sample_names <- gsub("^LFQ intensity ", "", lfq_cols)

  if (length(lfq_cols) == 0) {
    # Fallback to iBAQ or raw Intensity
    lfq_cols    <- grep("^iBAQ ", names(dt), value=TRUE)
    sample_names <- gsub("^iBAQ ", "", lfq_cols)
  }
  if (length(lfq_cols) == 0) {
    lfq_cols    <- grep("^Intensity ", names(dt), value=TRUE)
    sample_names <- gsub("^Intensity ", "", lfq_cols)
  }
  if (length(lfq_cols) == 0) {
    stop("No LFQ/iBAQ/Intensity columns found in MaxQuant proteinGroups.", call.=FALSE)
  }

  shared_samples <- intersect(smd$sample_id, sample_names)
  if (length(shared_samples) == 0) {
    stop("No sample IDs match MaxQuant column names. ",
         "Check that sample_id in metadata matches the experiment names in MaxQuant.",
         call.=FALSE)
  }
  matched_cols <- lfq_cols[sample_names %in% shared_samples]

  mat <- as.matrix(dt[, matched_cols, with=FALSE])
  mode(mat) <- "numeric"
  mat[mat == 0] <- NA
  colnames(mat) <- sample_names[sample_names %in% shared_samples]

  id_col <- intersect(c("Protein IDs","Majority protein IDs"), names(dt))[1]
  rownames(mat) <- as.character(dt[[id_col]])

  # ── Feature metadata ───────────────────────────────────────────────────────
  meta_keep <- intersect(c("Protein IDs","Majority protein IDs","Gene names",
                             "Protein names","Peptides",
                             "Razor + unique peptides","Unique peptides",
                             "Q-value","Score"),
                          names(dt))
  fmd <- dt[, meta_keep, with=FALSE]
  setnames_if <- function(dt, old, new) {
    if (old %in% names(dt)) data.table::setnames(dt, old, new)
  }
  fmd[, feature_id := dt[[id_col]]]
  setnames_if(fmd, "Gene names",            "gene_symbol")
  setnames_if(fmd, "Protein names",         "description")
  setnames_if(fmd, "Unique peptides",       "n_unique_peptides")
  setnames_if(fmd, "Razor + unique peptides","n_razor_peptides")
  setnames_if(fmd, "Majority protein IDs",  "protein_id")

  smd2 <- smd[sample_id %in% shared_samples]
  long <- data.table::melt(
    data.table::data.table(feature_id=rownames(mat), mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  quant_type <- if (grepl("^LFQ intensity", lfq_cols[1])) "lfq" else
                if (grepl("^iBAQ",          lfq_cols[1])) "ibaq" else "intensity"

  logger::log_info("import_maxquant: {nrow(mat)} proteins × {ncol(mat)} samples [{quant_type}]")

  pfd <- ProteoForgeData(
    raw_long        = long,
    raw_matrix      = mat,
    feature_metadata= fmd,
    sample_metadata = smd2,
    analysis_level  = config$input$analysis_level %||% "protein",
    quant_type      = quant_type,
    source_software = "maxquant",
    source_files    = c(quant_file, metadata_file),
    parameters      = list(quant_col=lfq_cols[1])
  )

  # ── Phospho sites (optional) ───────────────────────────────────────────────
  if (!is.null(sites_file) && file.exists(sites_file)) {
    pfd <- .import_mq_phospho(pfd, sites_file, smd2, config)
  }

  pfd
}

#' Expand MaxQuant Phospho (STY) Sites with multiplicity
#'
#' The most common bug in custom MaxQuant pipelines is not expanding the
#' `___1`, `___2`, `___3` multiplicity columns into separate rows.
#'
#' @keywords internal
.import_mq_phospho <- function(pfd, sites_file, smd, config) {
  logger::log_info("import_maxquant phospho: expanding multiplicity columns")

  sites <- .read_delim_auto(sites_file)

  # Remove decoys/contaminants
  for (col in intersect(c("Reverse","Potential contaminant"), names(sites))) {
    sites <- sites[get(col) != "+" | is.na(get(col))]
  }

  loc_cutoff <- config$phospho$localization_prob_cutoff %||% 0.75

  # Identify multiplicity columns (___1, ___2, ___3)
  mult_suffixes <- unique(gsub(".*(___(\\d+))$","___\\2",
                                grep("___\\d+$", names(sites), value=TRUE)))

  all_rows <- list()

  for (suf in mult_suffixes) {
    # Localisation probability column for this multiplicity
    loc_col <- grep(paste0("Localization prob.*", suf, "$"), names(sites),
                    value=TRUE, ignore.case=TRUE)[1]
    if (is.na(loc_col)) next

    quant_prefix <- "Intensity"
    quant_cols   <- grep(paste0("^", quant_prefix, ".*", suf, "$"), names(sites),
                          value=TRUE)
    sample_from_quant <- gsub(paste0("^", quant_prefix, " | ", suf, "$"), "",
                               quant_cols)

    if (length(quant_cols) == 0) next

    sub <- data.table::copy(sites)
    sub <- sub[!is.na(get(loc_col)) &
                 suppressWarnings(as.numeric(get(loc_col))) >= loc_cutoff]
    if (nrow(sub) == 0) next

    # Build site ID: UniProt_residue_position___multiplicity
    prot_col <- intersect(c("Protein","Proteins","Leading proteins"), names(sub))[1]
    res_col  <- intersect(c("Amino acid","Amino.acid"), names(sub))[1]
    pos_col  <- intersect(c("Position in peptide","Position","Position in protein"),
                           names(sub))[1]
    if (is.na(prot_col) || is.na(res_col)) next

    sub[, site_id := paste0(
      get(prot_col %||% "Protein"), "_",
      get(res_col  %||% "Amino acid"),
      if (!is.na(pos_col)) get(pos_col) else seq_len(.N),
      gsub("^___","_m", suf)
    )]
    sub[, multiplicity := as.integer(gsub("^___","", suf))]

    # Quant matrix
    shared <- intersect(smd$sample_id, sample_from_quant)
    if (length(shared) == 0) next
    matched <- quant_cols[sample_from_quant %in% shared]
    q_mat   <- as.matrix(sub[, matched, with=FALSE])
    q_mat[q_mat == 0] <- NA
    mode(q_mat) <- "numeric"
    colnames(q_mat) <- sample_from_quant[sample_from_quant %in% shared]
    rownames(q_mat)  <- sub$site_id

    all_rows[[suf]] <- list(mat=q_mat, fmd=sub[, .(feature_id=site_id,
                                                     multiplicity,
                                                     localization_prob=get(loc_col))])
  }

  if (length(all_rows) == 0) return(pfd)

  # Combine matrices
  all_ids <- unique(unlist(lapply(all_rows, function(x) rownames(x$mat))))
  all_sams <- colnames(all_rows[[1]]$mat)
  combined_mat <- matrix(NA, nrow=length(all_ids), ncol=length(all_sams),
                          dimnames=list(all_ids, all_sams))
  for (x in all_rows) {
    combined_mat[rownames(x$mat), colnames(x$mat)] <- x$mat
  }

  fmd_combined <- data.table::rbindlist(lapply(all_rows, `[[`, "fmd"), fill=TRUE)
  logger::log_info("MaxQuant phospho: {nrow(combined_mat)} Class-1 sites expanded from {length(mult_suffixes)} multiplicities")

  pfd@raw_matrix        <- combined_mat
  pfd@feature_metadata  <- fmd_combined
  pfd@analysis_level    <- "phosphosite"
  pfd@long              <- data.table::melt(
    data.table::data.table(feature_id=all_ids, combined_mat),
    id.vars="feature_id", variable.name="sample_id",
    value.name="intensity", variable.factor=FALSE)

  pfd
}
