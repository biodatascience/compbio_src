# builds the four biol_signal/tech_signal recovery plots (a PC1-PC2 scatter
# colored by biol_signal, a PC1 boxplot by effect_size, a PC1-PC2 scatter
# colored by tech_signal, and a PC1 vs tech_signal scatter) for both the
# "counts" and "vst" assays of 'se', and lays them out in a 4x2 grid with
# counts on the left, vst on the right
myFourPlots <- function(se) {
  make_col <- function(assay_name) {
    se_pca <- myPcaFunction(se, assay_name = assay_name, num_pcs = 4)

    p1 <- ggplot(as.data.frame(colData(se_pca)), aes(PC1, PC2, fill = biol_signal)) +
      geom_point(size = 3, shape = 21, color = "black") +
      scale_fill_distiller(palette = "RdBu")

    p2 <- ggplot(as.data.frame(rowData(se_pca)), aes(factor(effect_size), PC1)) +
      geom_boxplot()

    p3 <- ggplot(as.data.frame(colData(se_pca)), aes(PC1, PC2, fill = tech_signal)) +
      geom_point(size = 3, shape = 21, color = "black") +
      scale_fill_distiller(palette = "PRGn")

    p4 <- ggplot(as.data.frame(colData(se_pca)), aes(tech_signal, PC1)) +
      geom_point(size = 3)

    list(p1, p2, p3, p4)
  }

  counts_plots <- make_col("counts")
  vst_plots <- make_col("vst")

  wrap_plots(
    counts_plots[[1]], vst_plots[[1]],
    counts_plots[[2]], vst_plots[[2]],
    counts_plots[[3]], vst_plots[[3]],
    counts_plots[[4]], vst_plots[[4]],
    ncol = 2
  )
}

myFourPlots(se_pois)
myFourPlots(se_binom)
myFourPlots(se_nb)
