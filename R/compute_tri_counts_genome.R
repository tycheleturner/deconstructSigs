#' Compute trinucleotide counts for an arbitrary genome
#'
#' This helper computes the counts of every trinucleotide (e.g. "ACA")
#' in a BSgenome, restricted to contexts with a central base in C/T
#' (32 contexts total), matching the format of tri.counts.genome.
#'
#' @param bsg A BSgenome object.
#' @param include.seqnames Character vector of seqnames to include.
#'   Defaults to all seqnames(bsg). You can subset, e.g. only autosomes.
#' @param central.bases Vector of central bases to consider (default c("C","T")).
#'   Change only if you know what you’re doing; the deconstructSigs
#'   ecosystem assumes C/T as the reference bases for 96-context signatures.
#'
#' @return A data.frame with 32 rows and 1 column, rownames are trinucleotides
#'   like "ACA", "ACC", etc., suitable for use as tri.counts.method in
#'   whichSignatures().
#' @importFrom GenomeInfoDb seqnames
#' @export
compute_tri_counts_genome <- function(bsg, include.seqnames = GenomeInfoDb::seqnames(bsg), central.bases = c("C", "T")) {
  if (!inherits(bsg, "BSgenome")) {
    stop("bsg must be a BSgenome object.")
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required to compute trinucleotide counts. ",
         "Install it with BiocManager::install(\"Biostrings\").")
  }

  # Only keep seqnames that actually exist in the BSgenome
  include.seqnames <- intersect(include.seqnames, GenomeInfoDb::seqnames(bsg))
  if (length(include.seqnames) == 0L) {
    stop("No seqnames in 'include.seqnames' are present in 'bsg'.")
  }

  # Get sequences for selected chromosomes/contigs
  dna <- BSgenome::getSeq(bsg, include.seqnames)

  # Count all 64 possible 3-mers across the selected sequences
  # rows: sequences, cols: 3-mers; colnames like "AAA", "AAC", ...
  freqs <- Biostrings::oligonucleotideFrequency(dna, width = 3L, step = 1L, as.prob = FALSE)

  # Sum across sequences to get genome-wide counts
  genome_counts <- colSums(freqs)

  # Keep only contexts whose central base is in central.bases (default C/T)
  all_trimers <- names(genome_counts)
  mid_base <- substr(all_trimers, 2L, 2L)
  keep <- mid_base %in% central.bases
  genome_counts <- genome_counts[keep]

  # Re-order to a canonical 32-context order:
  # left in A,C,G,T; mid in central.bases; right in A,C,G,T
  base <- c("A", "C", "G", "T")
  contexts <- expand.grid(left = base, mid = central.bases, right = base, stringsAsFactors = FALSE)
  tri_names <- paste0(contexts$left, contexts$mid, contexts$right)

  # Make sure we actually have counts for all these trimers
  missing <- setdiff(tri_names, names(genome_counts))
  if (length(missing) > 0L) {
    stop("Missing trinucleotides in computed counts: ",
         paste(missing, collapse = ", "))
  }

  genome_counts <- genome_counts[tri_names]

  # Return data frame in the same shape as tri.counts.genome:
  # 32 rows x 1 column, rownames are contexts
  out <- data.frame(tri.counts = as.numeric(genome_counts), row.names  = tri_names, stringsAsFactors = FALSE)

  out
}
