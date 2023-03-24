## reproducible code for paper 

library(ggplot2)
library (tidyverse)
library(reshape2)
library(gtable)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(patchwork)
library(forcats)
library(ggbreak)
library(svglite)

# CSV containing tally of TP contaminants/mitochondria, FN and FP contaminants
input_scaffolds = "~/Documents/ROCKU/decontam_analyses/decontam_stats_scaffolds.csv" 
input_bases = "~/Documents/ROCKU/decontam_analyses/decontam_stats_bases.csv"


scaffolds <- read.csv(input_scaffolds,
                      header = TRUE,
                      sep = ",")

bases <- read.csv(input_bases, 
                  header = TRUE,
                  sep = ",")

## Checking for expected number of columns and rows. 
dim(scaffolds); dim(bases)

## Removing rows and columns that do not contain data. 
scaffolds_clean <- scaffolds[scaffolds$ground.truth != 0,c(1:7)] 
dim(scaffolds_clean)

plot_perc_scaffolds <- function (inscaffolds){
  
  ## Calculating the percent of scaffolds identified as TP, FP and FN. 
  ## Note 100% is not the sum of the number of scaffolds in the test genomes, but the sum TP, FN and FP. 
  denom <- inscaffolds$TP + inscaffolds$FP.contams + inscaffolds$FN.contams
  TP.perc <- (inscaffolds$TP/denom) * 100
  FP.contam.perc <- (inscaffolds$FP.contams/denom) * 100 
  FN.contam.perc <- (inscaffolds$FN.contams/denom) * 100
  
  ## Merging the original table with calculated percent values. 
  scaffolds_counts_perc <- cbind(inscaffolds,TP.perc, FP.contam.perc, FN.contam.perc)
  scaffolds_perc <- scaffolds_counts_perc[,c("assembly","TP.perc","FN.contam.perc","FP.contam.perc")]
  ## Rearranging table so the columns can be stacked in ggplot. 
  melt_scaffolds_perc <- melt(scaffolds_perc)
  
  ggplot(melt_scaffolds_perc, aes(x = factor(assembly, levels = scaffolds_perc$assembly), y = value, fill = variable, )) + 
    geom_col(position = "fill",  alpha = 0.8) + 
    coord_flip() +
    scale_fill_manual(values = c("darkcyan", "grey 65", "darkorange"),
                      labels = c("TP","FN contaminant","FP contaminant")) + 
    scale_y_reverse(labels = c("",25,50,75,100), expand = expansion(add=c(0.005,0.005)))+
    xlab(label = "Assembly") +
    ylab(label = "Percent of scaffolds") +
    labs (fill = "") +
    theme_bw() + 
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10, colour = c("violetred",
                                                       "dodgerblue1", "dodgerblue1", "dodgerblue1", "dodgerblue1",
                                                       "salmon", "salmon",
                                                       "slateblue3","slateblue3","slateblue3","slateblue3","slateblue3", "slateblue3",
                                                       "brown","brown",
                                                       "black")),
      legend.position = "none",
      plot.margin = margin(l=1, t=1,b=1,r=1, unit ="cm"),
      panel.grid.major.x = element_blank()
      ) 
  
}

bases_clean <- bases[,c("assembly","TP.perc", "FP.contam.perc","FN.contam.perc")]
dim(bases_clean)

plot_perc_bases <- function(inbases) { ## Generates a similar figure as above but with bases as opposed to scaffolds 
  
  maxTP <- max(inbases$TP.perc)
  
  melt_bases_clean <- melt(bases_clean)
  melt_bases_no_zero<- melt_bases_clean[melt_bases_clean$value != 0,]
  
  ggplot(melt_bases_no_zero, mapping = aes(x = assembly, y =value, fill = variable)) +
    geom_col(position = position_stack(reverse = TRUE), alpha = 0.8) +
    coord_flip () +
    scale_y_continuous(limits=c(0,maxTP),expand = expansion(add=c(0.005,0.005)), labels = c("",0.005,0.010,0.015,"")) + 
    ## The break in the axis is necessary because axis transformations were not working great, and I wanted to retain the outlier
    scale_y_break(breaks = c(0.02,0.17), scales = 0.2,space = 0, expand = FALSE, ticklabels = c(0.175,0.19)) + 
    scale_fill_manual(values = c("darkcyan", "darkorange", "grey 65"),
                    labels = c("TP","FP contaminant","FN contaminant")) +
    xlab(label = "Assembly") +
    ylab(label = "Percent of bases") + # this would be bases that were correctly identified, and bases misidentified (including false negatives)
    theme_bw() +
    labs(fill = "") +
    theme (
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.x.top =  element_text(colour = "white"),
      axis.text.x.bottom = element_text(size = 10),
      axis.ticks = element_blank(),
      plot.margin = margin(l=0, t=1,b=1,r=1, unit ="cm"),
      legend.text = element_text(size = 10),
      panel.grid.major.x =  element_blank()) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) 
  
  
}

scaffolds_perc_plot <- plot_perc_scaffolds(scaffolds_clean)
bases_perc_plot <- plot_perc_bases(bases_clean)

scaffolds_perc_plot
bases_perc_plot

ggsave(filename = "bases_perc_plot.svg", plot = bases_perc_plot, width = 6, height = 4)
ggsave(filename = "scaffolds_perc_plot.svg", plot = scaffolds_perc_plot, width = 6, height = 4)
