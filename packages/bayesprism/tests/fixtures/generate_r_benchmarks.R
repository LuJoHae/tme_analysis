# R script to generate ground-truth benchmark outputs from original BayesPrism R package
.libPaths(c("scratch/R_lib", .libPaths()))
devtools::load_all("scratch/BayesPrism/BayesPrism")

output_dir <- "packages/bayesprism/tests/fixtures/r_benchmark"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

load("scratch/BayesPrism/tutorial.dat/tutorial.gbm.rdata")

set.seed(123)

# Sample a balanced subset across all cell types for fast deterministic testing
selected_idx <- c(
  which(cell.type.labels == "tumor")[1:40],
  which(cell.type.labels == "myeloid")[1:30],
  which(cell.type.labels == "oligo")[1:20],
  which(cell.type.labels == "tcell")[1:20],
  which(cell.type.labels == "endothelial")[1:20],
  which(cell.type.labels == "pericyte")[1:20]
)

shared_genes <- intersect(colnames(sc.dat), colnames(bk.dat))[1:500]

sc_sub <- sc.dat[selected_idx, shared_genes]
cell_type_sub <- cell.type.labels[selected_idx]
cell_state_sub <- cell.state.labels[selected_idx]
bk_sub <- bk.dat[1:5, shared_genes]

# Run new.prism
myPrism <- new.prism(
  reference = sc_sub,
  input.type = "count.matrix",
  cell.type.labels = cell_type_sub,
  cell.state.labels = cell_state_sub,
  key = "tumor",
  mixture = bk_sub,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# Export collapsed and normalized reference phi
write.csv(myPrism@phi_cellState@phi, file.path(output_dir, "phi_cellState.csv"))
write.csv(myPrism@phi_cellType@phi, file.path(output_dir, "phi_cellType.csv"))
write.csv(myPrism@mixture, file.path(output_dir, "mixture_filtered.csv"))

# Run run.prism with small MCMC chain length for test benchmark
bp_res <- run.prism(
  prism = myPrism,
  n.cores = 1,
  update.gibbs = TRUE,
  gibbs.control = list(chain.length = 50, burn.in = 10, thinning = 2, seed = 123),
  opt.control = list(maxit = 100, optimizer = "MAP")
)

# Export initial and final posterior results
write.csv(bp_res@posterior.initial.cellType@theta, file.path(output_dir, "theta_initial_cellType.csv"))
write.csv(bp_res@posterior.theta_f@theta, file.path(output_dir, "theta_final_cellType.csv"))

cat("R benchmark fixtures generated successfully in", output_dir, "\n")
