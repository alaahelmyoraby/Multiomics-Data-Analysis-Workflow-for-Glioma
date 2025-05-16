---
author: Alaa Oraby
date: 2025-04-15
output:
  html_document: default
  pdf_document: default
title: code_markdown
---

`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)`

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

`{r cars} summary(cars)`

## Including Plots

You can also embed plots, for example:

`{r pressure, echo=FALSE} plot(pressure)`

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

``` {r}
library(readxl)
library(minfi)  # For methylation data analysis
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # Annotation for EPIC array
library(IlluminaHumanMethylationEPICmanifest)  # Manifest for EPIC array
library(FlowSorted.Blood.EPIC)  # Reference for blood cell types
library(GeneNet)  # For network analysis
library(limma)  # For linear modeling and normalization
library(GenomicRanges)
library(data.table)  # For data table operations
library(tidyverse)  # For data wrangling and visualization
library(zoo)  # For time series operations
library(lattice)  # For creating trellis graphics
library(ggplot2)  # For data visualization
library(reshape)  # For reshaping data
library(stringr)  # For string operations
library(RColorBrewer)  # For color palettes
library(mixOmics)  # For multivariate analysis (PLS-DA, etc.)
library(pROC)  # For ROC analysis
library(biomaRt)  # For gene annotation from Ensembl
library(org.Hs.eg.db)  # For human gene annotations
library(entropy)
library(data.table)
library(MOFA2)
library(impute)
library(psych)
library(basilisk.utils)
library(reticulate)
```

``` {r}
trans <- fread("num.mat_trans.imputed.logged.csv", data.table = FALSE)
epi <- fread("num.mat_epi.imputed.logged.csv", data.table = FALSE)
miRNA<- fread("num.mat_miRNA.imputed.logged.csv", data.table = FALSE)
annotation <- read.delim("annotation.txt", row.names=1)
```

``` {r}
rownames(trans)=trans$V1
trans=trans[,-1]

rownames(miRNA)=miRNA$V1
miRNA=miRNA[,-1]

rownames(epi)=epi$V1
epi=epi[,-1]

common_cols <- intersect(colnames(epi), colnames(annotation))
annotation <- annotation[, common_cols]

trans<-trans[,common_cols]
miRNA<-miRNA[,common_cols]
epi<- epi[,common_cols]

trans <- as.matrix(trans, data.table = FALSE)
miRNA <- as.matrix(miRNA, data.table = FALSE)
epi <- as.matrix(epi, data.table = FALSE)



my_data <- list(trans = trans, miRNA = miRNA, epi = epi)
save(my_data, file = "my_data.rda")

my_data_clean <- lapply(my_data, function(g) {
  g <- as.matrix(g)
  mode(g) <- "numeric"
  return(g)
})
```

``` {r}
MOFAobject <- create_mofa(my_data_clean)

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)

data_opts$scale_views
data_opts$scale_groups
data_opts$views

model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 15

model_opts$likelihoods

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
# Transpose your annotation
t_annotation <- t(annotation)

```

``` {r}
MOFAobject <- prepare_mofa(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)
```

``` {r}
MOFAobject <- create_mofa(my_data_clean)
MOFAobject <- prepare_mofa(MOFAobject)
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                          training_options = train_opts)
MOFAobject <- run_posterior(MOFAobject)
```

``` {r}
plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=15)
```

``` {r}
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
```

``` {r}
# Convert it to data.frame (this step is essential)
CLL_metadata <- as.data.frame(t_annotation)
# Add sample metadata
CLL_metadata$sample <- rownames(CLL_metadata)

samples_metadata(MOFAobject) <- CLL_metadata

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Tumor_purity","histological_type","gender","radiation_therapy","race","overall_survival","ethnicity"), 
                                  plot="log_pval"
)

```

``` {r}

plot_weights(MOFAobject,
             view = "trans",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_weights(MOFAobject,
             view = "miRNA",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
plot_weights(MOFAobject,
             view = "epi",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "trans",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
plot_top_weights(MOFAobject,
                 view = "epi",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject, 
             view = "miRNA", 
             factor = 1, 
             nfeatures = 10
)
```

``` {r}
plot_weights(MOFAobject,
             view = "trans",
             factor = 3,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_weights(MOFAobject,
             view = "miRNA",
             factor = 3,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
plot_weights(MOFAobject,
             view = "epi",
             factor = 3,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "trans",
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
plot_top_weights(MOFAobject,
                 view = "epi",
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject, 
             view = "miRNA", 
             factor = 3, 
             nfeatures = 10
)
```

``` {r}
# Plot embeddings for Factor 1 and Factor 2, split by view
plot_factor(MOFAobject, 
            factors = 1:3, 
            color_by = "histological_type", 
            dodge = TRUE) +
  ggtitle("Latent Factors 1, 2, and 3 by View")
```

``` {r}
plot_weights(MOFAobject, 
  view = "epi", 
  factor = 3, 
  nfeatures = 10,
  abs = F
)

plot_weights(MOFAobject, 
  view = "miRNA", 
  factor = 3, 
  nfeatures = 10,
  abs = F
)
plot_weights(MOFAobject, 
  view = "trans", 
  factor = 3, 
  nfeatures = 10,
  abs = F
)
```

``` {r}
plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "histological_type")
plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "histological_type",
  dodge = TRUE,
  add_violin = TRUE
)
```

``` {r}
metadata <- samples_metadata(MOFAobject)
metadata2<-metadata
# Create a new categorical column for plotting
metadata2$TumorPurityGroup <- ifelse(as.numeric(metadata2$Tumor_purity) > 0.85, "High", "Low")

# Update the metadata
samples_metadata(MOFAobject) <- metadata2

# Plot with the new grouping
plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "TumorPurityGroup")

plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "TumorPurityGroup",
  dodge = TRUE,
  add_violin = TRUE
)

```

``` {r}
# Plot with the new grouping
plot_factor(MOFAobject, 
            factors = 6, 
            color_by = "TumorPurityGroup")

plot_factor(MOFAobject, 
  factors = 6, 
  color_by = "TumorPurityGroup",
  dodge = TRUE,
  add_violin = TRUE
)
```

``` {r}
plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "AQP5",
  add_violin = TRUE,
  dodge = TRUE
)
# Example 1
plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "NAV2_cg20686479",
  add_violin = TRUE,
  dodge = TRUE
)

# Example 2
plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "LOC100126784_cg20686479",
  add_violin = TRUE,
  dodge = TRUE
)

# Example 3
plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "JAG1_cg26705720",
  add_violin = TRUE,
  dodge = TRUE
)

plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "hsa-mir-196a-1",
  add_violin = TRUE,
  dodge = TRUE
)
```

``` {r}
plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "TERT",
  add_violin = TRUE,
  dodge = TRUE
)
# Example 1
plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "ISM1_cg14060111",
  add_violin = TRUE,
  dodge = TRUE
)

# Example 2
plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "ISM1_cg12664209",
  add_violin = TRUE,
  dodge = TRUE
)

# Example 3
plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "GPR156_cg19093820",
  add_violin = TRUE,
  dodge = TRUE
)

plot_factor(MOFAobject, 
  factors = 3, 
  color_by = "hsa-mir-1262",
  add_violin = TRUE,
  dodge = TRUE
)
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 1,  #Transcriptomics
  features = 6,
  sign = "positive",
  color_by = "NAV2_cg20686479"
) + labs(y="Transcriptomics")

plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 1,  
  features = 6,
  sign = "negative",
  color_by = "NAV2_cg20686479"
) + labs(y="Transcriptomics")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "miRNA",
  factor = 1,  
  features = 6,
  sign = "positive",
  color_by = "AQP5"
) + labs(y="miRNA")

plot_data_scatter(MOFAobject, 
  view = "miRNA",
  factor = 1,  
  features = 6,
  sign = "negative",
  color_by = "AQP5"
) + labs(y="miRNA")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 1,  
  features = 6,
  sign = "positive",
  color_by = "ISM1_cg14060111"
) + labs(y="trans")

plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 1,  
  features = 6,
  sign = "negative",
  color_by = "ISM1_cg14060111"
) + labs(y="trans")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "epi",
  factor = 1,  
  features = 5,
  sign = "positive",
  color_by = "ISM1_cg12664209"
) + labs(y="Methylation")
plot_data_scatter(MOFAobject, 
  view = "epi",
  factor = 1,  
  features = 5,
  sign = "negative",
  color_by = "ISM1_cg12664209"
) + labs(y="Methylation")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "epi",
  factor = 1,  
  features = 5,
  sign = "positive",
  color_by = "GPR156_cg19093820"
) + labs(y="Methylation")

plot_data_scatter(MOFAobject, 
  view = "epi",
  factor = 1,  
  features = 5,
  sign = "negative",
  color_by = "GPR156_cg19093820"
) + labs(y="Methylation")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 3,  
  features = 6,
  sign = "positive",
  color_by = "hsa-mir-1262"
) + labs(y="Transcriptomics")
plot_data_scatter(MOFAobject, 
  view = "trans",
  factor = 3,  
  features = 6,
  sign = "negative",
  color_by = "hsa-mir-1262"
) + labs(y="Transcriptomics")
```

``` {r}
plot_data_scatter(MOFAobject, 
  view = "miRNA",
  factor = 1,  
  features = 6,
  sign = "positive",
  color_by = "hsa-mir-196a-1"
) + labs(y="miRNA")

plot_data_scatter(MOFAobject, 
  view = "miRNA",
  factor = 1,  
  features = 6,
  sign = "negative",
  color_by = "hsa-mir-196a-1"
) + labs(y="miRNA")
```

``` {r}
plot_data_heatmap(MOFAobject, 
  view = "trans",
  factor = 1,  
  features = 25,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  scale = "row"
)
```

``` {r}
plot_data_heatmap(MOFAobject, 
  view = "epi",
  factor = 1,  
  features = 25,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  scale = "row"
)
```

``` {r}
plot_data_heatmap(MOFAobject, 
  view = "miRNA",
  factor = 1,  
  features = 25,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  scale = "row"
)
```

``` {r}
library(circlize)
```

``` {r}
library(circlize)

# Extract top 10 features from each view (reduce clutter)
top_n <- 10
top_features_trans <- get_weights(MOFAobject, view = "trans", factor = 1, as.data.frame = TRUE) %>%
  arrange(desc(abs(value))) %>% head(top_n)

top_features_epi <- get_weights(MOFAobject, view = "epi", factor = 1, as.data.frame = TRUE) %>%
  arrange(desc(abs(value))) %>% head(top_n)

top_features_miRNA <- get_weights(MOFAobject, view = "miRNA", factor = 1, as.data.frame = TRUE) %>%
  arrange(desc(abs(value))) %>% head(top_n)

# Set group for color mapping
top_features_trans$group <- "Transcriptomics"
top_features_epi$group <- "Epigenomics"
top_features_miRNA$group <- "miRNA"

# Combine node metadata
all_nodes <- bind_rows(
  top_features_trans[, c("feature", "group")],
  top_features_epi[, c("feature", "group")],
  top_features_miRNA[, c("feature", "group")]
)

# Create a few clear links only (between trans & epi, trans & miRNA)
circos_data <- expand.grid(
  from = top_features_trans$feature,
  to = c(top_features_miRNA$feature, top_features_epi$feature)
)
circos_data$value <- 1

# Colors
group_colors <- c(
  Transcriptomics = "#99c2ff",
  Epigenomics = "#ffddb3",
  miRNA = "#b3e6cc"
)

# Match feature colors
feature_groups <- setNames(all_nodes$group, all_nodes$feature)
feature_colors <- setNames(group_colors[feature_groups], names(feature_groups))

# Initialize circos plot
circos.clear()
circos.par(gap.after = c(rep(2, length(feature_colors)-1), 8))
circos.initialize(factors = names(feature_colors), xlim = cbind(rep(0, length(feature_colors)), rep(1, length(feature_colors))))

# Add outer colored track
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
                       panel.fun = function(x, y) {
                         sector_name <- get.cell.meta.data("sector.index")
                         circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                                     sector_name, facing = "clockwise", niceFacing = TRUE,
                                     adj = c(0, 0.5), cex = 1)
                       },
                       bg.col = feature_colors, bg.border = NA)

# Plot links (faint blue lines for clarity)
for (i in 1:nrow(circos_data)) {
  circos.link(circos_data$from[i], 0.5,
              circos_data$to[i], 0.5,
              col = "#1f78b4AA", border = NA, lwd = 0.8)
}
```
