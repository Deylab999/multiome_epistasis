# This script plots an interaction term which summarizes the interaction terms
# observed from our epistasis models.
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2) # for plotting
library(scales)

# iterate through different model output files
all.epistasis.models <- data.frame()
for (i in 1:32) {

    # create filename
    filename <- paste0('epistasis_models_', i, '.csv')

    # read in epistasis models from batch
    epistasis.models <- read.csv(filename)

    # add model results to combined data frame
    all.epistasis.models <- rbind(all.epistasis.models, epistasis.models)
}

# filter out NA values from model results
epistasis.models <- all.epistasis.models
epistasis.models <- epistasis.models[complete.cases(epistasis.models), ]

# set column names for model results
colnames(epistasis.models) <- c(
    'gene', 'enhancer.1', 'enhancer.2', 'intercept', 'beta.1.estimate',
    'beta.1.pvalue', 'beta.2.estimate', 'beta.2.pvalue', 'interaction.estimate',
    'interaction.pvalue'
)

# compute -log(p-value)
epistasis.models$neglog10p <- -1 * log10(epistasis.models$interaction.pvalue)

# compute the distance between the ATAC-seq peaks (enhancers)
get_start_coord <- function(x) {
    strsplit(x, '-')[[1]][2]
}

enhancer.1.start <- sapply(epistasis.models$enhancer.1, get_start_coord)
enhancer.2.start <- sapply(epistasis.models$enhancer.2, get_start_coord)

names(enhancer.1.start) <- NULL
names(enhancer.2.start) <- NULL

epistasis.models$enhancer.1.start <- as.numeric(enhancer.1.start)
epistasis.models$enhancer.2.start <- as.numeric(enhancer.2.start)

epistasis.models$distance <- abs(
    epistasis.models$enhancer.1.start - 
    epistasis.models$enhancer.2.start
)

# calculate log of distance (for plotting)
epistasis.models$log.dist <- log10(epistasis.models$distance)

# create volcano plot
plot <- ggplot(epistasis.models, aes(
    interaction.estimate,
    neglog10p,
    colour = log.dist
    )) +
    geom_point() +
    scale_colour_gradient2(
        low = muted('red'),
        high = muted('blue'),
        midpoint = floor(min(epistasis.models$log.dist)),
        limits=c(
            min(epistasis.models$log.dist),
            max(epistasis.models$log.dist)
        ),
        breaks = seq(
            floor(min(epistasis.models$log.dist)),
            ceiling(max(epistasis.models$log.dist))
        ),
        labels = seq(
            floor(min(epistasis.models$log.dist)),
            ceiling(max(epistasis.models$log.dist))
        )
    ) +
    theme_classic()

ggsave(
    'volcano_plot.png',
    plot
)
