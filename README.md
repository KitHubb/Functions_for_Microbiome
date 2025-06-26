# Microbiome Analysis Toolkit (R/phyloseq-based)

This repository provides a comprehensive R-based toolkit for microbiome data analysis using the `phyloseq` package. It covers alpha/beta diversity, taxonomic analysis, statistical testing, and publication-ready visualization including volcano plots and bar plots. Below is a structured README describing each function module and usage example.

---

## Requirements

### R packages

```
install.packages(c("phyloseq", "ggplot2", "dplyr", "ggtext", "ggprism", "ggrepel", "patchwork", "broom", "vegan", "glue", "ggpubr", "tidyr", "purrr", "rvg", "officer", "combinat", "RColorBrewer", "GGally"))
```

```
source("./250429_microbiome_analysis_functions.R")
```

## 1. Alpha Diversity

### Function: `create_alpha_diversity_df()`

* Calculates richness metrics: Observed, Shannon, Chao1, Inverse Simpson.
* Optional: Faith's PD (requires tree).

```r
alpha_df <- create_alpha_diversity_df(physeq_obj, tree = TRUE)
```

### Function: `alpha_plot()`

* Plots violin + boxplot + significance tests (Kruskal, Wilcoxon).

```
alpha_plot(alpha_df, x_ = "Group", y_ = "Shannon")
```

### Function: `alpha_df()`

* Outputs CI & p-values for each index using Kruskal-Wallis and Wilcoxon.

```
res <- alpha_df(df = alpha_df, group = "Group")
```

---

## 2. Beta Diversity

### Function: `beta_plot()`

* Computes PERMANOVA, ANOSIM, PERMDISP.
* Visualizes PCoA with ellipse.

```
beta_out <- beta_plot(phyloseq, type = "Group", type_col = c("red", "blue"))
```

### Function: `beta_plot.envfit()`

* Extension of `beta_plot()` with `vegan::envfit()` support and global PERMANOVA.

```
beta_env <- beta_plot.envfit(phyloseq, type = "Group", envfit = TRUE, formula = "Age + BMI", PERMANOVA = TRUE)
```

---

## 3. Taxonomic Abundance & Visualization

### Function: `Abund_cal()`

* Computes mean, sd, se, CI per taxon and group.

```
Abund_cal(ps.glom, tax_level = "Genus", group = "Group", path = "./results/")
```

### Function: `taxa_plot()` / `tax_plot.sp()`

* Visualizes top taxa at any level or within a genus.
* Automatically color-coded by phylum.

```
tax_plot_out <- taxa_plot(melt, taxa = "Genus", tax_otu = otu_list, x_axis = "SampleID")
```

---

## 4. Mean Abundance Comparison (with Statistics)

### Function: `mean_plot()` / `mean_plot.sp()`

* Generates comparative plots: bar, log2FC, p-value, heatmap.
* Returns `patchwork` plot and statistics.

```
res <- mean_plot(melt, tax_level = "Genus", tax_otu = top_genera, type = "Group", palet = c("#E31A1C", "#1F78B4"))
```

---

## 5. Volcano Plot

### Function: `volcano_plot()`

* Computes Wilcoxon p-values, log2FC, and highlights significant taxa.

```
vol_out <- volcano_plot(ph_glom, tax_rank = "Species", comparision = "Group", compar1 = "A", compar2 = "B")
```

---

## 6. BLAST Annotation Processing

### Function: `blast_arrange()`

* Extracts Genus, Species names from BLAST `sscinames` output.

```
blast_table <- blast_arrange(blast_output, ASV = "Query")
```

---

## 7. Prevalence & Abundance Filtering

### Function: `filter_phyloseq_by_prevalence_abundance()`

* Filters OTUs based on prevalence and average abundance thresholds.

```
ps_filtered <- filter_phyloseq_by_prevalence_abundance(ps, 0.001, 0.1)
```

---

## 8. Miscellaneous Utilities

* `save_ggplot_to_pptx()`: Saves ggplot as editable PowerPoint.
* `ps_tax_update()`: Updates ASV ID based on abundance and reassigns taxonomy table.
* `seq_to_fas()`: Converts sequence list to FASTA.

---

## Author

This script was developed for microbiome analysis pipelines based on `phyloseq`. Contact the author for details or collaboration.

## License

MIT License.


