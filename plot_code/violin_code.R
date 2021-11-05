      ## ---------------------------
      ##
      ## Script name: violin_code.R
      ##
      ## Purpose of script: Create violin plots for in-vivo tnseq data
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

# load libraries
      library(here)
      library(tidyverse)
      
# load in data
      load(here("resamp_data.RData"))
      

# read in gene list
myco_genes <- scan(here("gene_lists/mycobactin_mendum.txt"), what="", sep="\n")
pdim_genes <- scan(here("gene_lists/pdim_mendum.txt"), what="", sep="\n")
carb_genes <- scan(here("gene_lists/carb_trans_mendum.txt"), what = "", sep="\n")
chol_kstr2_genes <- scan(here("gene_lists/chol_kstr2_mendum.txt"), what = "", sep="\n")
mce_genes <- scan(here("gene_lists/mce4_operon.txt"), what = "", sep="\n")
chol_ring_genes <- scan(here("gene_lists/chol_ring_cat_mendum.txt"), what = "", sep="\n")

## create Node gene group violin plot

node_data <- data_ortho %>%
  filter(TAs > 5) %>%
  filter(call !="ES") %>%
  filter(tissue == "Node") %>%
  group_by(filename) %>%
  mutate(sample_median = median(log2FC)) %>%
  ungroup() %>%
  mutate(norm_lfc = log2FC-sample_median) %>%
  group_by(Orf) %>%
  mutate(mean_lfc = mean(log2FC))

myco_node <- node_data %>%
    filter(Orf %in% myco_genes) %>%
    mutate(group = "mycobactin synth")
 
pdim_node <- node_data %>%
  filter(Orf %in% pdim_genes) %>%
  mutate(group = "pdim synth")

carb_n <- node_data %>%
  filter(Orf %in% carb_genes) %>%
  mutate(group = "carb transport")
  
chol_kstr2_n <- node_data %>%
  filter(Orf %in% chol_kstr2_genes) %>%
  mutate(group = "kstr2 regulon")
  
mce_n <- node_data %>%
  filter(Orf %in% mce_genes) %>%
  mutate(group = "chol import/mce4")

chol_ring_n <- node_data %>%
  filter(Orf %in% chol_ring_genes) %>%
  mutate(group = "chol ring catab")

base_nodes <- data_ortho %>%
  filter(TAs >5 ) %>%
  mutate(group = "background") %>%
  filter(call !="ES") %>%
  filter(tissue == "Node") %>%
  group_by(filename) %>%
  mutate(sample_median = median(log2FC)) %>%
  ungroup() %>%
  mutate(norm_lfc = log2FC-sample_median) %>%
  ungroup() %>%
  mutate(mean_lfc = NA) #set as NA so segments won't show
  
median(base_nodes$log2FC)

median(base_nodes$norm_lfc)
#0

data_node <- NULL
data_node <- full_join(myco_node, pdim_node)
data_node <- full_join(data_node, carb_n)
data_node <- full_join(data_node, chol_kstr2_n)
data_node <- full_join(data_node, mce_n)
data_node <- full_join(data_node, chol_ring_n)
data_node <- full_join(data_node, base_nodes)
data_node$group <- as.factor(data_node$group)

node <-ggplot(data=data_node, group=group) +
        geom_violin( aes(x = group, y = norm_lfc, fill = group),
                     trim=F, draw_quantiles = 0.5) +
        geom_segment( aes(
              x    = match(group, levels(group)) - 0.1,
			        xend = match(group, levels(group)) + 0.1,
		          y    = mean_lfc,
			        yend = mean_lfc ),
		          col='white') +
        scale_fill_manual(values=c("grey50", "red3",
                          "#00FF00","dodgerblue3", "yellow2",
                          "violetred2","turquoise3",
                          "chocolate4")) +
        #scale_fill_brewer(type = "qual", palette = "Spectral") +
        ggtitle(label = "Gene groups: Nodes", subtitle = "white = mean LFC for Orf; black = median LFC for group, normalised for median of background logFC") +
        theme(panel.background = element_rect(
                  fill = "grey75", 
                  colour = "grey75"),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(
                        angle = 90, 
                        vjust = 0.5, 
                        hjust=1,
                        size=10, 
                        face = "bold.italic"),
              plot.title = element_text(hjust = 0.5), 
              legend.position="none" ) 
node

