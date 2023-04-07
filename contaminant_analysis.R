library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(ggbreak)
library(ggsci)
library(svglite)

## Analyzing the blast results of contaminant sequences. 

comp_IDs <- read.csv(file = "/Users/nadolinabrajuka/Documents/ROCKU/decontam_analyses/compiled_IDs_MAR3.csv", 
                     header = FALSE)
head(comp_IDs)
colnames(comp_IDs) <- c("VGP_ID","scaffold","acc_num", "species", "coverage","size")
commonnames <- read.csv(file = "/Users/nadolinabrajuka/Documents/ROCKU/decontam_analyses/fullnames_tolids.csv", 
                        header = TRUE)

summary(comp_IDs)
summary(as.factor(comp_IDs$species))

## Identifying the top contaminant ID for each unique VGPID-scaffold pair (unfiltered dataframe)
top_IDs <- comp_IDs[which(!duplicated(comp_IDs[,c("VGP_ID","scaffold")])),] 
dim(top_IDs)

## Creating a dataframe where contaminants identified as 'almond' or 'no identification' have been removed 
filt_IDs <- comp_IDs %>% 
  filter(species != 'almond') %>% 
  filter(species != 'no identification')
dim(filt_IDs)
## Identifying top ID from filtered dataframe 
filt_topIDs<- filt_IDs[which(!duplicated(filt_IDs[,c("VGP_ID","scaffold")])),] 
dim(filt_topIDs)

## checking to make sure all the VGP-ID and scaffolds line up between the two tables 
top_IDs$scaffold == filt_topIDs$scaffold; top_IDs$VGP_ID == filt_topIDs$VGP_ID

## merging the original top IDs and the filtered top IDs for comparison 
merge_df <- full_join(top_IDs, filt_topIDs, by = c("VGP_ID","scaffold"))
dim(merge_df); head(merge_df);tail(merge_df)

former_almonds <- merge_df[merge_df$species.x == "almond",] ##isolating rows of the merged table where almond was the initial ID 
dim(former_almonds)
summary(as.factor(former_almonds$VGP_ID))
summary(as.factor(former_almonds$species.y))
summary(former_almonds$coverage.x); summary(as.numeric(former_almonds$coverage.y))

filt_topIDs$coverage<-round(as.numeric(filt_topIDs$coverage),0)


## Grouping taxa together for the figure + removing the killer whale IDs 
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Delftia")] <- "Delftia"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Escherichia")] <- "E. coli"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Staphylococcus")] <- "S. aureus"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Leptospira")] <- "L. borgpetersenii"
filt_topIDs$species[filt_topIDs$species == "Shigella flexneri"] <- "S. flexneri"
filt_topIDs_nowhale  <- filt_topIDs[which(filt_topIDs$species != "killer whale"),]

final_df<- merge(filt_topIDs_nowhale , commonnames, by.x = "VGP_ID", by.y = "TOLID", all.x = TRUE) 
dim(final_df)
  
contam_ID_size_cov <- ggplot(final_df, mapping = aes (x = common.name, 
                                            y = species, 
                                            size = size, 
                                            color = coverage
                                            )) +
  geom_point(position = "jitter", alpha = 0.8) + 
  theme_bw() +
  scale_size_continuous(breaks = c(2000, 10000, 50000, 100000, 1000000), 
                        range =c(0.01, 10),
                        trans = "log2", 
                        labels = c(2, 10, 50, 100, 1000), 
                        name = "Contaminant size (Kb)") +
  scale_color_viridis_c(trans = "reverse", 
                        option = "F",
                        begin = 0.2,
                        end = 0.9)  +
  geom_hline(yintercept  = seq(1.5,12,1), linetype = "dotted", colour = "grey20") +
  geom_vline(xintercept  = seq(1.5,12,1), linetype = "dotted", colour = "grey20") +
   labs(color = "Alignment coverage (%)", y = "Contaminant ID", x = "Assembly")  +
  theme(axis.text.x =  element_text(angle = 30,
                                    hjust = 1),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold", vjust = 4))

ggsave(filename = "/Users/nadolinabrajuka/Documents/ROCKU/decontam_analyses/contam_ID_size_cov.svg", plot = contam_ID_size_cov, width = 7, height = 5)
  


  
