# Helper functions for gene_finder.R: sampling gene models from a TxDb,
# building labeled segments, and turning DNA sequence into the hidden
# states X_t and the observed 3-mers Y_t for an HMM.

# The hidden states of the model. Coding sequence (CDS) is read in
# codons, so there are three coding states c1, c2, c3 for the position
# within the codon, and the start and stop codons each get three
# states. Introns can fall after any codon position, so to keep track
# of the reading frame there are three copies of the intron states:
# the phase p is the codon position (c1, c2, c3) just before the
# intron, and after the acceptor the CDS resumes at the next codon
# position.
#
# Each intron starts with a chain of n_donor donor states and ends with
# a chain of n_acceptor acceptor states, one state per position, so
# that the emissions can learn the splice site motifs beyond 'GT' and
# 'AG'. The "simple" model has n_donor = n_acceptor = 1.
#
# Returns a list with:
#   states          names of the states
#   allowed         allowed transitions from each state (all others
#                   have probability 0)
#   emission_group  states in the same group share emission
#                   probabilities, e.g. donor1_2, donor2_2, donor3_2
#                   all emit as "donor_2"

make_model <- function(n_donor=1, n_acceptor=1) {
  intron_states <- function(p) {
    c(paste0("donor", p, "_", 1:n_donor),
      paste0("intronic", p),
      paste0("acceptor", p, "_", 1:n_acceptor))
  }
  states <- c("intergenic",
              "start1", "start2", "start3",
              "c1", "c2", "c3",
              intron_states(1), intron_states(2), intron_states(3),
              "stop1", "stop2", "stop3")

  allowed <- list(
    intergenic = c("intergenic", "start1"),
    start1 = "start2", start2 = "start3", start3 = "c1",
    c1 = c("c2", "donor1_1"),
    c2 = c("c3", "donor2_1"),
    c3 = c("c1", "donor3_1", "stop1"),
    stop1 = "stop2", stop2 = "stop3", stop3 = "intergenic"
  )
  for (p in 1:3) {
    # each intron state goes to the next one in the chain,
    # and the intronic state can also stay put
    chain <- c(intron_states(p), paste0("c", p %% 3 + 1))
    for (i in seq_len(length(chain) - 1)) {
      allowed[[chain[i]]] <- chain[i + 1]
    }
    intronic <- paste0("intronic", p)
    allowed[[intronic]] <- c(intronic, allowed[[intronic]])
  }

  emission_group <- sub("^(donor|intronic|acceptor)[123]", "\\1", states)

  list(states=states, allowed=allowed, emission_group=emission_group)
}

# Sample + strand protein-coding transcripts (one per gene) with fewer
# than max_cds CDS pieces. Keep only canonical CDS (starting with ATG,
# ending with a stop codon, length a multiple of 3), with CDS pieces of
# at least 6 bp and introns of at least 20 bp (room for the donor and
# acceptor chains), and that, padded by 'pad' bp on either side,
# overlap no other gene.
# Returns the CDS span (start codon to stop codon) of each transcript.

sample_genes <- function(txdb, genome, n_genes, n_sample, max_cds, pad) {
  # all genes, for checking overlaps with other genes
  all_genes <- genes(txdb) |>
    GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")

  # CDS pieces and introns of each transcript, named by tx_id
  cds_by_tx <- cdsBy(txdb, by="tx")
  intr_by_tx <- intronsByTranscript(txdb)

  # positive strand protein-coding transcripts with a single gene ID
  tx <- transcripts(txdb, columns=c("tx_id","gene_id")) |>
    GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse") |>
    filter(strand == "+", lengths(gene_id) == 1,
           tx_id %in% names(cds_by_tx)) |>
    mutate(gene_id = unlist(gene_id))

  candidates <- tx |>
    slice_sample(n=n_sample) |>
    filter(!duplicated(gene_id), # one transcript per gene
           gene_id %in% all_genes$gene_id)

  cds <- cds_by_tx[as.character(candidates$tx_id)]
  cds_seq <- extractTranscriptSeqs(genome, cds)
  first_codon <- as.character(subseq(cds_seq, 1, 3))
  last_codon <- as.character(subseq(cds_seq, width(cds_seq) - 2))

  # shrink each transcript to its CDS span, and check the CDS
  candidates <- candidates |>
    mutate(start = min(start(cds)), end = max(end(cds)),
           n_cds = lengths(cds),
           canonical = first_codon == "ATG" &
                       last_codon %in% c("TAA","TAG","TGA") &
                       width(cds_seq) %% 3 == 0,
           cds_6bp = all(width(cds) >= 6))

  # introns within the CDS span that are shorter than 20 bp
  short_introns <- stack(intr_by_tx[as.character(candidates$tx_id)], "tx_id") |>
    filter_by_overlaps(candidates) |>
    filter(width < 20)

  candidates <- candidates |>
    filter(n_cds < max_cds, canonical, cds_6bp,
           !tx_id %in% short_introns$tx_id)

  # padded CDS span should overlap only its own gene (ignoring strand)
  sampled <- candidates |>
    mutate(start = start - pad, end = end + pad)
  sampled <- sampled %>%
    mutate(n_overlaps = count_overlaps(., all_genes)) |>
    filter(n_overlaps == 1)

  stopifnot(length(sampled) >= n_genes)
  sampled <- head(sampled, n_genes)

  # return the original (un-padded) CDS spans
  candidates |> filter(tx_id %in% sampled$tx_id)
}

# Labeled segments (intergenic, cds, intronic) for each transcript,
# sorted along the genome. Only introns within the CDS span are used.
# Introns longer than max_intron are trimmed to their first
# (max_intron - 10) bp plus their last 10 bp, so that the splice donor
# and acceptor sites are kept.

make_segments <- function(txdb, genes_tx, pad, max_intron) {
  ids <- as.character(genes_tx$tx_id)

  introns <- stack(intronsByTranscript(txdb)[ids], "tx_id") |>
    filter_by_overlaps(genes_tx)
  long_introns <- introns |> filter(width > max_intron)
  introns <- bind_ranges(
    introns |> filter(width <= max_intron),
    long_introns |> mutate(width = max_intron - 10), # start of intron
    long_introns |> anchor_end() |> mutate(width = 10) # end of intron
  )

  bind_ranges(
    intergenic = bind_ranges(
      genes_tx |> flank_upstream(pad),
      genes_tx |> flank_downstream(pad)
    ) |> select(tx_id),
    cds = stack(cdsBy(txdb, by="tx")[ids], "tx_id") |> select(tx_id),
    intronic = introns |> select(tx_id),
    .id = "feature"
  ) |>
    mutate(tx_id = as.integer(as.character(tx_id))) |>
    arrange(seqnames, start)
}

# Hidden state for each base of the concatenated sequence, as an
# integer indexing model$states.
#   intergenic             label of the segment
#   c1, c2, c3             position in the codon, counting CDS bases
#                          along the gene (so the frame is correct
#                          across introns)
#   start1-3, stop1-3      the first and last three CDS bases
#   donor, intronic,       first n_donor bases, middle, and last
#   acceptor               n_acceptor bases of an intron, with the
#                          phase (1, 2, 3) of the last CDS base before
#                          the intron

segments_to_states <- function(segments, model) {
  n_donor <- sum(grepl("^donor1_", model$states))
  n_acceptor <- sum(grepl("^acceptor1_", model$states))
  feature <- as.character(segments$feature)
  w <- width(segments)
  x <- rep(feature, w)
  gene <- rep(segments$tx_id, w)

  # k: position along the CDS of each gene, n_cds: CDS length
  is_cds <- x == "cds"
  k <- ave(as.integer(is_cds), gene, FUN=cumsum)
  n_cds <- ave(as.integer(is_cds), gene, FUN=sum)
  codon_pos <- (k - 1) %% 3 + 1
  stop_k <- k - (n_cds - 3) # 1, 2, 3 for the last three CDS bases

  # intron phase: codon position of the last CDS base before each base
  last_cds <- cummax(ifelse(is_cds, seq_along(x), 0))
  phase <- codon_pos[pmax(last_cds, 1)]

  # splice sites, from the features of neighboring segments
  prev <- c(NA, head(feature, -1))
  nxt <- c(tail(feature, -1), NA)
  seg_end <- cumsum(w) # position in concatenated sequence
  seg_start <- seg_end - w + 1
  intron_start <- seg_start[feature == "intronic" & prev == "cds"]
  intron_end <- seg_end[feature == "intronic" & nxt == "cds"]
  # positions of the donor and acceptor chains, and the step in the chain
  donor <- rep(intron_start, each=n_donor) + 0:(n_donor - 1)
  donor_k <- rep(1:n_donor, length(intron_start))
  acceptor <- rep(intron_end, each=n_acceptor) - (n_acceptor - 1):0
  acceptor_k <- rep(1:n_acceptor, length(intron_end))

  is_intron <- x == "intronic"
  x[is_cds] <- paste0("c", codon_pos[is_cds])
  first3 <- is_cds & k <= 3
  x[first3] <- paste0("start", k[first3])
  last3 <- is_cds & stop_k >= 1
  x[last3] <- paste0("stop", stop_k[last3])
  x[is_intron] <- paste0("intronic", phase[is_intron])
  x[donor] <- paste0("donor", phase[donor], "_", donor_k)
  x[acceptor] <- paste0("acceptor", phase[acceptor], "_", acceptor_k)

  match(x, model$states)
}

# Hidden states X and centered 3-mers Y for the concatenated sequence
# of a set of segments (see dna_to_kmers for why X is trimmed).

segments_to_xy <- function(segments, genome, model) {
  dna <- getSeq(genome, segments) |> unlist()
  x <- segments_to_states(segments, model)
  list(x = x[-c(1, length(x))], y = dna_to_kmers(dna))
}

# Centered 3-mers: Y_t is the 3-mer made of bases t-1, t, t+1. The
# first and last bases have no centered 3-mer, so Y has length
# length(dna) - 2 and lines up with X[-c(1, length(X))].
# Returns a factor with the 64 3-mers as levels.

dna_to_kmers <- function(dna) {
  bases <- c("A","C","G","T")
  b <- match(strsplit(as.character(dna), "")[[1]], bases)
  stopifnot(!anyNA(b)) # no N's
  L <- length(b)
  left <- b[1:(L-2)]
  mid <- b[2:(L-1)]
  right <- b[3:L]
  kmers <- mkAllStrings(bases, 3) # AAA, AAC, AAG, ..., TTT
  factor(kmers[16 * (left - 1) + 4 * (mid - 1) + right], levels=kmers)
}

# Parameters estimated by counting, when the hidden states X are known.
# Emissions are counted by model$emission_group, so e.g. the three
# copies of each intron state share the same emission probabilities.

count_params <- function(x, y, model) {
  n <- length(model$states)
  pi <- tabulate(x[1], n)
  A <- table(factor(head(x, -1), 1:n), factor(tail(x, -1), 1:n))
  A <- unclass(A / rowSums(A))
  dimnames(A) <- list(model$states, model$states)
  group <- model$emission_group
  counts <- table(factor(group[x], unique(group)), y)
  prob <- unclass(counts / rowSums(counts))[group,]
  dimnames(prob) <- list(model$states, levels(y))
  list(pi=pi, A=A, B=list(prob=prob))
}

# Coarser labels for evaluating a state path: intergenic and intronic
# are both "noncoding", and the donor and acceptor states are merged

state_group <- function(x, model) {
  s <- model$states[x]
  s[s %in% c("c1","c2","c3")] <- "coding"
  s[grepl("^(intergenic|intronic)", s)] <- "noncoding"
  s[grepl("^start", s)] <- "start"
  s[grepl("^stop", s)] <- "stop"
  s[grepl("^donor", s)] <- "donor"
  s[grepl("^acceptor", s)] <- "acceptor"
  factor(s, c("noncoding","start","coding","donor","acceptor","stop"))
}
