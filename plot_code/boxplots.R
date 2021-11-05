      ## ---------------------------
      ##
      ## Script name: boxplots.R
      ##
      ## Purpose of script: Generate boxplots for groups of genes from in-vivo 
      ## tn-seq data
      ##
      ## Author: Jennifer J. Stiens
      ##
      ## Date Created: 2021-11-04
      ##
      ## Copyright (c) Jennifer J. Stiens, 2021
      ## Email: j.j.stiens@gmail.co.uk
      ##
      ## ---------------------------
      ##
      ## Load Packages:
      ##
      require(here)
      require(tidyverse)
      ##   
      ## Load Functions:
      ##
      ## ---------------------------
      ##

  
## load libraries
library(here)
library(tidyverse)

# load in data
load(here("resamp_data.RData"))

rd1_genes <-scan(here("gene_lists", "RD1_bcg_microti_genes.txt"), what="", sep="\n")
length(rd1_genes)

# any with ES or sites <5
data_ortho %>% 
  filter(Orf %in% rd1_genes) %>% 
  #filter(innoc_call=="ES")
  filter(TAs <5) %>%
  group_by(Orf, TAs) %>%
  distinct(Orf, TAs) 
# MB3902 (1), MB3904 (4), MB3905 (1), and MB3908 (3) have less than 5 insertion sites

rd_1 <- data_ortho %>%
  filter(Orf %in% rd1_genes) %>%
  filter(call!="ES") %>%
  #filter(tissue == "Lung") %>%
  mutate(signif = `Adj. p-value` < 0.05)
rd_1$signif <- as.factor(rd_1$signif)

tik_seq <- c(rep("grey40", 8), "violetred", "grey40", "violetred", "violetred","grey40", "grey40", "violetred", "grey40")

p_tog <- ggplot(rd_1, aes(x=Orf, y=log2FC, fill=tissue)) +
  geom_boxplot(outlier.shape = NA) +   # hides outliers 
  scale_fill_manual(values=c("deepskyblue3", "darkseagreen3")) +
  guides(fill = FALSE) +    #removes legend for fill in boxplot
  geom_jitter(alpha=0.8, size=1, width = .1,
              aes(color=signif)) +
  scale_color_manual(name = "Adj. p-value <0.05", 
                     labels = c("Non-significant", "Significant"), 
                     values = c("grey60", "purple3")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"),
        axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust=1, size=10, face = "bold.italic", color = tik_seq),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label="RD1, BCG and M.Microti genes") +
  facet_wrap( ~tissue)
p_tog

## Boxplots Mtbvac genes and novel genes

mtbvac_genes <- c("MB0780", "MB0781", "MB2955")
novel_genes <- c("MB1033", "MB2991c", "MB3747" )

#are any of these genes ES or <5 insertions?
data_ortho %>%
  filter(Orf %in% mtbvac_genes | Orf %in% novel_genes) %>%
  #filter(Sites <5) NONE
  filter(call == "ES") #none

#both plots together (facet wrap)
mtbvac <- data_ortho %>%
  filter(Orf %in% mtbvac_genes) %>%
  filter(call!="ES") %>%
  #filter(tissue == "Lung") %>%
  mutate(signif = `Adj. p-value` < 0.05)
mtbvac$signif <- as.factor(mtbvac$signif)

p <- ggplot(mtbvac, aes(x=Orf, y=log2FC, fill=tissue)) +
  geom_boxplot(outlier.shape = NA) +   # hides outliers 
  scale_fill_manual(values=c("deepskyblue3", "darkseagreen3")) +
  guides(fill = FALSE) +    #removes legend for fill in boxplot
  geom_jitter(alpha=0.8, size=2, width = .1,
              aes(color=signif)) +
  scale_color_manual(labels = c("Non-significant", "Significant"), values = c("grey60", "purple3"), name=c("Adj p-value <0.05")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"),
        axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust=1, size=10, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label="MTBVAC genes") +
  facet_wrap( ~tissue)
p

#novel genes
novel <- data_ortho %>%
  filter(Orf %in% novel_genes) %>%
  mutate(signif = `Adj. p-value` < 0.05)
novel$signif <- as.factor(novel$signif)

# are any es or <5 sites?
novel %>% filter(TAs <5) #no
novel %>% filter(call=="ES") #no

novel_labels <-c("MB1033/Rv1006", "MB2991c/Rv2967c", "MB3747/Rv3720")

p <- ggplot(novel, aes(x=Orf, y=log2FC, fill=tissue)) +
  geom_boxplot(outlier.shape = NA) +   # hides outliers 
  scale_fill_manual(values=c("deepskyblue3", "darkseagreen3")) +
  guides(fill = "none") +    #removes legend for fill in boxplot
  geom_jitter(alpha=0.8, size=2, width = .1,
              aes(color=signif)) +
  scale_color_manual(labels = c("Non-significant", "Significant"), values = c("grey60", "purple3"), name=c("Adj p-value <0.05")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"),
        axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust=1, size=10, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label="novel genes") +
  scale_x_discrete(breaks = novel_genes, labels= novel_labels) +
  facet_wrap( ~tissue)
p

