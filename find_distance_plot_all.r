
get_probability_table = function(prob_files,myGene,area,cond, inDir){
  patern = paste0(myGene,"_",cond,"_rep_._",area,"_",".")
  if(cond =="both"){
    patern = paste0(myGene,".+","_",area,"_.+")
  }
  myGene_prob_files = prob_files[grep(patern,prob_files, perl = T)]
  
  
  
  for(i in 1:length(myGene_prob_files)){
    f = myGene_prob_files[i]
    rep = strsplit(f,"_")[[1]][4]
    type = strsplit(f,"_")[[1]][6]
    mycond = strsplit(f,"_")[[1]][2]
    if(type == "allProbability.csv"){
      type= "real"
    }
    tmp = read.csv(paste0(inDir,f), row.names = 1)
    tmp$cond = mycond
    tmp$area = area
    tmp$rep = rep
    tmp$gene = myGene
    tmp$type = type
    if (i == 1){
      finalDF = tmp
    }else{
      finalDF = rbind(finalDF,tmp)
    }
  }
  return(finalDF)
}

get_maen_sd = function(prob_table){
  mean_tab = aggregate(prob_table$probability, by = list( prob_table$type, prob_table$distance), mean)
  sd_tab = aggregate(prob_table$probability, by = list( prob_table$type, prob_table$distance), 
                     function(x) sd(x)/sqrt(length(x)))
  merge_df = merge(mean_tab,sd_tab, by = c("Group.1", "Group.2"), suffixes= c(".mean",".sd") )
  names(merge_df) = c("type","distance","mean","sd")
  merge_df$info = paste(merge_df$type,merge_df$combination,merge_df$distance,sep = "_")
  return(merge_df)
}

plot_grap = function(maen_sd_table,myGene,cond,area){
  p = ggplot(maen_sd_table, aes(y= mean, x= distance, color = type))+
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=4,
                  position=position_dodge(0.05)) +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    ggtitle(paste(myGene,cond,area)) #+
    #ylim(0.35,0.65)
  plot(p)
}




one_comb_all_anlysis = function(outDir,cond, myGene, area= "all", maxDist = 1000, myDist = F){
  

  inDir = paste0(outDir,"/tmp/")
  
  
  prob_files = list.files(inDir,pattern = "allProbability.csv")
  prob_table = get_probability_table(prob_files,myGene = "Egr2",area = "all",cond = "1h", inDir)
  
  maen_sd_table = get_maen_sd(prob_table)
  

  maen_sd_table_forPlot = maen_sd_table[maen_sd_table$distance <=  maxDist,] 

  pdf(paste0(outDir,"/",myGene,"_distance_probability_plot.pdf"), width = 12, height = 6)
  plot_grap(maen_sd_table_forPlot,myGene,cond,area)
  dev.off()
}

#################################


