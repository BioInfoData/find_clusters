#setwd("/path/to/find_clusters")


options(stringsAsFactors = F )

library("rmarkdown")
library("raster") 
library("ggplot2") 
library(RColorBrewer)

source("clust_find_distance_function.r")
source("clust_select_cells_function.r")
source("find_distance_plot_all.r")


### params ###

input_file = "input_example/_merged_1h_Arc_Nr4a1_.csv"
outDir = "output_example/"
cond = "1h"
myGene = "Egr2"
pos_cut = 6
do_shuffle = FALSE
reps = c(1:3)


### create output dirs ####

outDir_tmp = paste0(outDir,"/tmp/")
dir.create(outDir)
dir.create(outDir_tmp)


### get probability for positive neighbors for each positive cell in the data ###

find_distance_cells(input_file,outDir_tmp,cond,myGene,pos_cut, do_shuffel = FALSE)
find_distance_cells(input_file,outDir_tmp,cond,myGene,pos_cut, do_shuffel = TRUE)

### plot results of positive fraction environment ###

one_comb_all_anlysis(outDir, cond , myGene)

### select cells in positive environment for cluster analysis ###

select_cells(input_file,outDir_tmp,cond, myGene)

### find clusters in the data and calculate expression levels in each group of clusters ###

find_cluster = "clust_find_clusters_function.rmd"

for (numRep in reps){
  outFile = paste0(outDir,"/",myGene,"_",cond,"_",numRep,".html")
  render(find_cluster, output_dir = outDir, output_file = outFile,
         params = list(outDir_tmp = outDir_tmp, cond = cond, 
                       myGene =myGene, numRep = numRep, do_shuffle =do_shuffle ))
}



