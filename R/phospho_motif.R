#' Phosphomotif enrichment analysis
#'
#' Runs `rmotifx` on up-regulated and down-regulated phosphosites using
#' the detected phosphoproteome as the background. A ±7 amino-acid window
#' is used by default.
#'
#' @param phospho_result_dt `data.table` with site-level DE results.
#' @param fasta_file        Optional path to a FASTA file for sequence
#'   extraction. If `NULL`, sequences must be in `feature_metadata`.
#' @param window            Half-window size in amino acids. Default 7.
#' @param min_occurrences   Minimum occurrences for a motif. Default 5.
#' @param p_cutoff          p-value threshold. Default 0.000001.
#'
#' @return A named list: `up_motifs`, `down_motifs` (data.tables), or `NULL`.
#'
#' @examples
#' \dontrun{
#' motifs <- run_phospho_motif(phospho_result_dt)
#' }
#'
#' @export
run_phospho_motif <- function(phospho_result_dt,
                               fasta_file      = NULL,
                               window          = 7L,
                               min_occurrences = 5L,
                               p_cutoff        = 1e-6) {
  if (!requireNamespace("rmotifx", quietly=TRUE)) {
    collect_warning("Package 'rmotifx' not installed. Motif analysis skipped.",
                    "run_phospho_motif")
    return(NULL)
  }

  seq_col <- intersect(c("sequence_window","flanking_seq","Sequence window",
                          "flanking_sequence"),
                        names(phospho_result_dt))[1]
  if (is.na(seq_col)) {
    collect_warning("No sequence window column found. Motif analysis skipped.",
                    "run_phospho_motif")
    return(NULL)
  }

  .do_motif <- function(seqs, bg_seqs, label) {
    if (length(seqs) < min_occurrences) return(NULL)
    tryCatch({
      res <- rmotifx::motifx(seqs, bg_seqs,
                              central.res  = "S|T|Y",
                              min.seqs     = min_occurrences,
                              pval.cutoff  = p_cutoff,
                              verbose      = FALSE)
      if (!is.null(res) && nrow(res) > 0) {
        dt <- data.table::as.data.table(res)
        dt[, group := label]
        dt
      }
    }, error=function(e) NULL)
  }

  all_seqs  <- na.omit(phospho_result_dt[[seq_col]])
  up_seqs   <- na.omit(phospho_result_dt[direction=="up",   get(seq_col)])
  down_seqs <- na.omit(phospho_result_dt[direction=="down", get(seq_col)])

  list(
    up_motifs   = .do_motif(up_seqs,   all_seqs, "up"),
    down_motifs = .do_motif(down_seqs, all_seqs, "down")
  )
}
