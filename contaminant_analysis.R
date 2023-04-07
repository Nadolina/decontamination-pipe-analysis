library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(ggbreak)
library(ggsci)
library(svglite)

## contaminant analysis but with the sizes of all the contaminants 

comp_IDs <- read.csv(file = "/Users/nadolinabrajuka/Documents/ROCKU/decontam_analyses/compiled_IDs_MAR3.csv", 
                     header = FALSE)
head(comp_IDs)
colnames(comp_IDs) <- c("VGP_ID","scaffold","acc_num", "species", "coverage","size")

summary(comp_IDs)
summary(as.factor(comp_IDs$species))

## subsetting rows where the contaminant species is not labeled 'almond' or 'no identification'
not_almonds <- comp_IDs[which(comp_IDs$species != 'almond'),] 
not_noID <- not_almonds[which(not_almonds$species != 'no identification'),]
dim(not_noID)

filt_topID_ind <- which(!duplicated(not_noID[,c("VGP_ID","scaffold")])) ##isolating the index for the top non-almond species ID for each unique VGP_ID-scaffold pair 
filt_topIDs <- not_noID[filt_topID_ind,] ## top ID for each unique VGP_ID/scaffold pair that isn't almond or "no identification"
dim(filt_topID)

top_IDs_3 <- comp_IDs[which(!duplicated(comp_IDs[,c("VGP_ID","scaffold")])),] ##CORRECTLY ISOLATED THE TOP ID FOR EACH UNIQUE ID-SCAFFOLD PAIR, USE THIS ONE 
## this would include almonds 

## checking to make sure all the VGP-ID and scaffolds line up between the two tables 
top_IDs_3$scaffold == filt_topIDs$scaffold; top_IDs_3$VGP_ID == filt_topIDs$VGP_ID

## merging the original top IDs and the filtered top IDs for comparison 
merge_df <- full_join(top_IDs_3, filt_topIDs, by = c("VGP_ID","scaffold"))
dim(merge_df); head(merge_df);tail(merge_df)

former_almonds <- merge_df[merge_df$species.x == "almond",] ##isolating rows of the merged table where almond was the initial ID 
dim(former_almonds)
summary(as.factor(former_almonds$VGP_ID))
summary(as.factor(former_almonds$species.y))
summary(former_almonds$coverage.x); summary(as.numeric(former_almonds$coverage.y))

filt_topIDs$coverage<-round(as.numeric(filt_topIDs$coverage),0)

filt_topIDs_size <- filt_topIDs[which(filt_topIDs$size < 100000),] ## filtering for contaminants less than 100kb in size

ggplot(data = filt_topIDs, mapping = aes (x = VGP_ID, 
                                          y = size/1000, 
                                          alpha = 0.5,
                                          size = coverage, 
                                          colour = species)) + 
  geom_point(position = "jitter") +
  scale_y_break(breaks = c(170,320,380,3500), scales = c(0.2,2), expand = TRUE) + 
  ylim(0,3700) +
  scale_size(range = c(1,8))


# with the points "dodging"
ggplot(filt_topIDs_size, mapping = aes (x = VGP_ID, 
                                        y = coverage, 
                                        size = size, 
                                        colour = species,
                                        alpha = 0.7)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Alignment coverage")+ 
  labs(color = "Contaminant ID", size = "Contaminant size") +
  guides( alpha = "none")

# with the points "jittered" 

ggplot(filt_topIDs_size, mapping = aes (x = VGP_ID, 
                                        y = coverage, 
                                        size = size, 
                                        colour = species,
                                        alpha = 0.7)) + 
  geom_point(position = "jitter") +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(filt_topIDs_size, mapping = aes (x = VGP_ID, 
                                          y = size, 
                                          size = coverage, 
                                          colour = species,
                                          alpha = 0.7)) +
  geom_point() +
  scale_size(range =c (1,10)) + 
  theme_bw() 

## grouping the taxa together 

filt_topIDs_size$species[str_detect(filt_topIDs_size$species, pattern = "Delftia")] <- "Delftia"
filt_topIDs_size$species[str_detect(filt_topIDs_size$species, pattern = "Escherichia")] <- "E. coli"
filt_topIDs_size$species[str_detect(filt_topIDs_size$species, pattern = "Staphylococcus")] <- "S. aureus"
filt_topIDs_size_nowhale <- filt_topIDs_size[which(filt_topIDs_size$species != "killer whale"),]

ggplot(filt_topIDs_size_nowhale, mapping = aes (x = VGP_ID, 
                                                y = coverage, 
                                                color = species,
                                                size = size,
                                                alpha = 0.7))  +
  geom_point(position = "jitter") +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Coverage") +
  labs(size = "Contaminant size",
       colour = "Contaminant ID") + 
  guides(alpha = "none")


ggplot(filt_topIDs_size_nowhale, mapping = aes (x = VGP_ID, 
                                        y = coverage, 
                                        shape = species,
                                        alpha = 0.7))  +
  geom_point(position = "jitter", data = filt_topIDs_size_nowhale %>% filter(size < 10000), color = "black", size = 3) +
  geom_point(position = "jitter", data = filt_topIDs_size_nowhale %>% filter(size >= 10000), color = "red", size = 3) + 
  geom_point(position = "jitter", data = filt_topIDs_size_nowhale %>% filter(size >= 50000), color = "blue", size = 3) +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Coverage") +
  labs(size = "Contaminant size",
       colour = "Contaminant ID") + 
  guides(alpha = "none")


install.packages("treemap")
library(treemap)

## something is weird with the coverage scale 
datafr <- filt_topIDs_size_nowhale[,c(4,5,6)]
treemap(datafr, index = "species", vSize = "size", vColor = "coverage", type = "value", palette = "RdYlBu", aspRatio = 2)

filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Delftia")] <- "Delftia"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Escherichia")] <- "E. coli"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Staphylococcus")] <- "S. aureus"
filt_topIDs$species[str_detect(filt_topIDs$species, pattern = "Leptospira")] <- "L. borgpetersenii"
filt_topIDs$species[filt_topIDs$species == "Shigella flexneri"] <- "S. flexneri"
filt_topIDs_nowhale  <- filt_topIDs[which(filt_topIDs$species != "killer whale"),]

commonnames <- read.csv(file = "fullnames_tolids.csv", header = TRUE)

full_filt_FINAL<- merge(filt_topIDs_nowhale , commonnames, by.x = "VGP_ID", by.y = "TOLID", all.x = TRUE) ## SOMETHING WENT WRONG WITH THIS MERGE
dim(full_filt_FINAL)

ggplot(filt_topIDs_nowhale, mapping = aes (x = VGP_ID, 
                                                y = coverage, 
                                                shape = species,
                                                alpha = 0.7))  +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size < 10000), color = "orange", size = 3) +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 10000), color = "red", size = 3) + 
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 50000), color = "blue", size = 3) +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 100000), color = "black", size = 3) +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Coverage") +
  labs(size = "Contaminant size",
       colour = "Contaminant ID") + 
  guides(alpha = "none")

ggplot(filt_topIDs_size_nowhale, mapping = aes (x = VGP_ID, 
                                                y = size,
                                                size = coverage,
                                                color = species,
                                                alpha = 0.8)) +
  geom_point(position ="jitter") +
  theme_bw() +
  scale_size_continuous(trans = "reverse", range = c(2,7)) +
  scale_color_uchicago() +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Size") +
  labs(size = "Alignment coverage",
       colour = "Contaminant ID") + 
  guides(alpha = "none")


####FINAL OPTIONS ---------


# FIGURE 1 
ggplot(filt_topIDs_size_nowhale, mapping = aes (x = VGP_ID, 
                                        y = coverage, 
                                        size = size, 
                                        colour = species,
                                        alpha = 0.7)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Alignment coverage")+ 
  labs(color = "Contaminant ID", size = "Contaminant size") +
    scale_colour_uchicago() +
  guides( alpha = "none")

# FIGURE 2 - same as Figure 1 but with a jitter instead of dodge; more clarity on density, but lose resolution on S.flexneri (bStrHab1)
ggplot(filt_topIDs_nowhale, mapping = aes (x = VGP_ID, 
                                                y = coverage, 
                                                size = size, 
                                                colour = species,
                                                alpha = 0.7)) + 
  geom_point(position = "jitter") +
  scale_size_continuous(trans = "log2") +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Alignment coverage")+ 
  labs(color = "Contaminant ID", size = "Contaminant size") +
  scale_colour_uchicago() +
  guides( alpha = "none")

# FIGURE 3 - binned sizes, no colour legend as of present 
ggplot(filt_topIDs_nowhale, mapping = aes (x = VGP_ID, 
                                           y = coverage, 
                                           shape = species,
                                           alpha = 0.7))  +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size < 10000), aes(color = "orange"), size = 3) +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 10000), aes(color = "red"), size = 3) + 
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 50000), aes(color = "blue"), size = 3) +
  geom_point(position = "jitter", data = filt_topIDs_nowhale %>% filter(size >= 100000), aes(color = "black"), size = 3) +
  scale_size_binned(range = c(1,8), breaks = c(1000,10000,25000,75000)) +
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Alignemnt coverage") + 
  xlab("VGP ID")+ ## alternate name "Individual" "Specimen" etc. 
  labs(color = "Contaminant size", ) + 
  guides(alpha = "none") + 
  scale_color_identity(guide = "legend")

# FIGURE 4 - size as size. 

ggplot(filt_topIDs_size_nowhale, mapping = aes (x = VGP_ID, 
                                                y = species, 
                                                size = size, 
                                                color = coverage)) +
  geom_point(position = "jitter") + 
  scale_size_binned(range = c(1,10), breaks = c(1000,10000,25000,50000,75000)) + 
  geom_vline(xintercept = seq(1.5,12.5,1), linetype = "dotted", colour = "grey20") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Alignemnt coverage") + 
  xlab("VGP ID")+ ## alternate name "Individual" "Specimen" etc. 
  labs(size = "Contaminant size", color = "Contaminant size") + 
  guides(alpha = "none") 
  
contam_ID_size_cov <- ggplot(full_filt_FINAL, mapping = aes (x = common.name, 
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

ggsave(filename = "contam_ID_size_cov.svg", plot = contam_ID_size_cov, width = 7, height = 5)
  


  
