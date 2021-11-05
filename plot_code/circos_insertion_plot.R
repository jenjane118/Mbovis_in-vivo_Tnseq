      ## ---------------------------
      ##
      ## Script name: circos_insertion_plot.R
      ##
      ## Purpose of script: To create a circos-style plot showing insertion sites 
      ## over the entire Mbovis genome
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
      require(circlize)
      ##   
      ## Load Functions:
      ##
      ## ---------------------------
      ##
      
      library(circlize)
      library(dplyr)
      
      # using wig with all insertion sites included, including non-permissive sites
      f.wig <- read.table(file=here("insertion_files", "tpp_MbA27.wig"), header = FALSE, sep=" ", skip =2)
      df <- cbind(f.wig[,1],f.wig) #adds extra column with dupl position for start and end
      head(df)
      colnames(df) <- c("start","end", "value")
      df$logvalue <- log10((df$value+0.01))
      quantile(df[,3])
      ymax.f <- quantile(df[,3], probs=c(0.9))
      df$value <- ifelse(df$value>ymax.f,ymax.f,df$value)
      quantile(df[,3])
      ymax.f
      head(df)
      
      #read in the genes from the bed file
      beddata<- read.table(file=here("data_files", "LT708304_updated.bed"), header=FALSE, sep="\t")
      #keep the gene id
      beddata$geneid <- gsub(pattern=".*_","",beddata$V5)
      #keep the gene name (where there isn't one, the locus tag that folllows will take the place of the gene name)
      beddata$genename <- gsub(pattern="gene=|;locus_tag.*", "", beddata$V5)
      beddata$genename <- gsub(pattern="locus_tag.*_", "", beddata$genename)
      head(beddata)
      #only one genome displayed around the whole circle
      bovis = c("Mbovis_LT708304.1")
      colnames(beddata)[1] <- "chr"
      #region.size is the number of bases shown in the circle
      #region.size <- 10000
      region.size <- 4349904
      #how many genes in bedplus and bedminus can be shown within the region.size?
      numgenes <- dim(subset(beddata, beddata[,3] < region.size))[1]
      #limit the plotting of the wig to the region size we are showing
      endregion <- dim(subset(df, df[,2]< region.size))[1]
      
      
      # have to run whole plot together in one chunk
      circos.clear()
      
      png("MbA27_circos.png", width = 960, height=960)
      #1) initialise plot
      circos.initialize(sectors=beddata[1:numgenes,"chr"], xlim=c(0,region.size))
      #2) plot genome
      circos.trackPlotRegion(factors=beddata[1:numgenes,"chr"], 
                             ylim=c(0,120), 
                             track.height = 0.3, panel.fun=function(x,y){circos.genomicLines(region=df[1:endregion,1:2],                                       value=df[1:endregion,], 
                                                                                             numeric.column=3, 
                                                                                             type = "h",
                                                                                             col="blue")
                               circos.yaxis(labels.cex=0.75)
                             })
      #3) draw axis with nt size markers around outside of circle
      circos.genomicAxis(
        h = "top",
        major.at = NULL,
        labels = NULL,
        major.by = NULL,
        tickLabelsStartFromZero = TRUE,
        #labels.cex = 0.5*par("cex"),
        labels.cex = par("cex"),
        sector.index = get.current.sector.index(),
        track.index = get.current.track.index())
      
      title(main="MBovis TnSeq Input Library Insertions")
      dev.off()
      