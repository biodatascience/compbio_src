# Use an HMM to find coding sequence, introns, and the boundaries
# between them in DNA sequence from human genes.
#
# We sample 100 protein-coding genes from hg38 (UCSC knownGene) on the
# + strand, with fewer than 15 CDS pieces and a canonical CDS (ATG ...
# stop codon). UTRs are removed: each gene runs from the start codon to
# the stop codon, padded with 100 bp of intergenic sequence on either
# side, and does not overlap any other gene. Introns are trimmed to
# 200 bp.
#
# Half of the genes are used for training: since we know their hidden
# states, we estimate the emission probabilities B by counting. On the
# other half (test), we hold B fixed and use Baum-Welch to estimate pi
# and A, then find the Viterbi path and compare it to the truth.
#
# Hidden states X_t (see make_model in gene_finder_helpers.R):
#   intergenic, start1-3, c1-3 (codon position), a chain of donor
#   states, intronic, a chain of acceptor states (three copies each, for
#   the intron phase), stop1-3
# Observations Y_t: the 3-mer centered at position t (bases t-1, t, t+1)

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(plyranges)
})
source("hmm_wrapper.R")
source("gene_finder_helpers.R")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

pad <- 100 # bp of intergenic sequence before and after each gene
max_intron <- 200 # longer introns are trimmed to this length

# simple model: one donor and one acceptor state (19 states)
# advanced model: chains of 4 donor and 6 acceptor states (43 states)
model <- make_model(n_donor=1, n_acceptor=1)
# model <- make_model(n_donor=4, n_acceptor=6)
states <- model$states
n <- length(states)

# Sample genes and build labeled segments ----------------------------

set.seed(1)
genes_tx <- sample_genes(txdb, Hsapiens, n_genes=100, n_sample=500,
                         max_cds=15, pad=pad)
segments <- make_segments(txdb, genes_tx, pad=pad, max_intron=max_intron)
table(segments$feature)

# split the genes in half, into training and test sets
train_ids <- sample(genes_tx$tx_id, 50)
train <- segments |> filter(tx_id %in% train_ids) |>
  segments_to_xy(Hsapiens, model)
test <- segments |> filter(!tx_id %in% train_ids) |>
  segments_to_xy(Hsapiens, model)
table(factor(states[train$x], states))

# Train: parameters by counting on the labeled training genes --------

truth <- count_params(train$x, train$y, model)

# The most likely 3-mers for each emission group: 'ATG' at start2,
# 'TAA', 'TAG', 'TGA' at stop2, '?GT' at the first donor state, 'AG?'
# at the last acceptor state
groups <- unique(model$emission_group)
B_groups <- truth$B$prob[match(groups, model$emission_group),]
rownames(B_groups) <- groups
t(apply(B_groups, 1, function(p) names(sort(p, decreasing=TRUE))[1:3]))

# Test: Baum-Welch for pi and A, with B fixed from training ----------

# start in the intergenic state
init_pi <- function(n) c(1, rep(0, n-1))

# only the transitions in 'allowed' are possible, the rest are zero
# and stay zero in EM. Long runs of intergenic and intronic, and the
# codon cycle c1 -> c2 -> c3 -> c1, start with higher weight.
init_A <- function(n) {
  A <- matrix(0, nrow=n, ncol=n, dimnames=list(states, states))
  for (s in names(model$allowed)) {
    A[s, model$allowed[[s]]] <- 1
  }
  A["intergenic", "intergenic"] <- 99
  A[cbind(paste0("intronic", 1:3), paste0("intronic", 1:3))] <- 99
  A["c1", "c2"] <- A["c2", "c3"] <- A["c3", "c1"] <- 99
  A / rowSums(A)
}

# Down-weight the emissions: each base appears in three overlapping
# 3-mers, so the HMM would count its evidence about three times over,
# making the emissions overconfident relative to the transitions. Then
# the path switches states too often, in short runs, and Baum-Welch
# estimates switching probabilities in A that are too high. Raising B
# to the power 1/2 fixes much of this (see emission_weight in
# hmm_wrapper.R). Try emission_weight=1 to compare.
fit <- hmm_fit(test$y, n_states=n, family="multinomial",
               init_pi=init_pi, init_A=init_A,
               init_B=function(y, n) truth$B, update_B=FALSE,
               emission_weight=1/2, verbose=TRUE)

plot(fit$loglik, type="b", xlab="iteration", ylab="log likelihood")

# compare A from Baum-Welch on the test genes with A counted on the
# training genes, for the states with more than one possible transition
dimnames(fit$A) <- list(states, states)
branching <- names(model$allowed)[lengths(model$allowed) > 1]
for (s in branching) {
  to <- model$allowed[[s]]
  cat("from", s, "\n")
  print(round(rbind(baum_welch = fit$A[s, to], counted = truth$A[s, to]), 4))
}

# Compare the Viterbi path to the true states of the test genes ------

# intergenic and intronic merged as noncoding
truth_group <- state_group(test$x, model)
viterbi_group <- state_group(fit$viterbi, model)
table(truth=truth_group, viterbi=viterbi_group)

# fraction of positions with the right group
mean(truth_group == viterbi_group)

# among positions correctly called coding, the right codon position
is_coding <- truth_group == "coding" & viterbi_group == "coding"
mean(test$x[is_coding] == fit$viterbi[is_coding])
