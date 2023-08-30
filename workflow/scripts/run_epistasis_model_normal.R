# This script runs a linear model with log-normalized expression to test for
# interaction between ATAC-seq peaks. This is for comparison with the
# Poisson GLM model.
#
# Author: Karthik Guruvayurappan

library(Seurat) # for normalization
library(dplyr) # for data frames

# read in RNA-seq matrix, ATAC-seq matrix, and cell-level metadata
rna <- readRDS(snakemake@input[[1]])
atac <- readRDS(snakemake@input[[2]])
metadata <- readRDS(snakemake@input[[3]])

# filter RNA and ATAC matrices to only include genes and peaks
# with at-least 5% non-zero counts
rna <- rna[(rowSums(rna > 0) / ncol(rna)) >= 0.05, ]
atac <- atac[(rowSums(atac > 0) / ncol(atac)) >= 0.05, ]

# create Seurat object and normalize RNA-seq data
seurat.rna <- CreateSeuratObject(counts = rna)
seurat.rna <- NormalizeData(seurat.rna)
rna <- GetAssayData(seurat.rna, 'data')

# filter RNA and ATAC matrices to only include desired cell type
rna <- rna[, metadata$celltype == snakemake@params[['celltype']]]
atac <- atac[, metadata$celltype == snakemake@params[['celltype']]]

# read in enhancer pairs determined from SCENT
enhancer.pairs <- read.csv(snakemake@input[[4]])

# create dataframe to store epistasis model results
results <- data.frame(matrix(ncol = 10, nrow = nrow(enhancer.pairs)))

# iterate through list of enhancer pairs and run epistasis model
for (i in 1:nrow(enhancer.pairs)) {

    # get name of enhancer 1, enhancer 2, and gene
    enhancer.1 <- enhancer.pairs$enhancer_1[i]
    enhancer.2 <- enhancer.pairs$enhancer_2[i]
    gene <- enhancer.pairs$gene[i]

    # get ATAC and RNA data for peaks and gene
    enhancer.1.atac <- atac[enhancer.1, ]
    enhancer.2.atac <- atac[enhancer.2, ]
    gene.rna <- rna[gene, ]

    # fit linear model
    mdl <- lm(
        gene.rna ~ enhancer.1.atac * enhancer.2.atac
    )

    # store model values and write to data frame
    mdl.values <- summary(mdl)$coefficients

    # isolate individual values of interest
    intercept <- NA
    beta.1.estimate <- NA 
    beta.2.estimate <- NA 
    interaction.estimate <- NA 

    beta.1.pvalue <- NA 
    beta.2.pvalue <- NA 
    interaction.pvalue <- NA 

    # add error checking for possible NA values
    if ("(Intercept)" %in% rownames(mdl.values)) {
        intercept <- mdl.values['(Intercept)', 'Estimate']
    }

    if ("enhancer.1.atac" %in% rownames(mdl.values)) {
        beta.1.estimate <- mdl.values['enhancer.1.atac', 'Estimate']
        beta.1.pvalue <- mdl.values['enhancer.1.atac', 'Pr(>|t|)']


    }

    if ("enhancer.2.atac" %in% rownames(mdl.values)) {
        beta.2.estimate <- mdl.values['enhancer.2.atac', 'Estimate']
        beta.2.pvalue <- mdl.values['enhancer.2.atac', 'Pr(>|t|)']
    }

    if ("enhancer.1.atac:enhancer.2.atac" %in% rownames(mdl.values)) {
        interaction.estimate <- mdl.values[
            'enhancer.1.atac:enhancer.2.atac',
            'Estimate'
        ]
        interaction.pvalue <- mdl.values[
            'enhancer.1.atac:enhancer.2.atac',
            'Pr(>|t|)'
        ]
    }

    # add model results to results data frame
    mdl.vector <- c(gene, enhancer.1, enhancer.2, intercept, beta.1.estimate,
                    beta.1.pvalue, beta.2.estimate, beta.2.pvalue,
                    interaction.estimate, interaction.pvalue)
    results[i, ] <- mdl.vector
}

# write results to output file
write.csv(results, snakemake@output[[1]], row.names = FALSE)
