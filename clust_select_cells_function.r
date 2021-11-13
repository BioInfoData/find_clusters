




############################## functions ################


get_data_rep = function(all_data,myGene,cond,numRep,myRadius, do_shuffle = FALSE,outDir_tmp){
  filePositive = paste0(myGene,"_",cond,"_rep_",numRep,"_all_ratio.csv")
  if(do_shuffle){
    filePositive = paste0(myGene,"_",cond,"_rep_",numRep,"_all_SHUFFLE_ratio.csv") 
  }
  
  positive_data = read.csv(paste0(outDir_tmp,filePositive))
  names(positive_data)[1] = "cellID"
  all_data$cellID = row.names(all_data)
  
  positive_with_info = merge(positive_data,all_data, by = "cellID")
  
  filt_pos = positive_with_info[!is.na(positive_with_info[,myRadius]),]
  filt_pos = filt_pos[filt_pos[,myRadius] != "not_pos",]
  filt_pos[,myRadius] = as.numeric(filt_pos[,myRadius])
  
  pos_cut = min(filt_pos[,myGene])
  pos_fraction = sum(positive_with_info[,myGene] >= pos_cut)/ nrow(positive_with_info)
  
  if(do_shuffle){
    pos_fraction = sum(positive_with_info[,myGene] >= gene_cutoff)/ nrow(positive_with_info)
  }
  
  increas_cut =0.8
  
  ind = filt_pos[,myRadius] >= increas_cut
  ind2 = filt_pos[,myRadius] <= 0
  selected_cells = filt_pos[ind,]
  not_selected_cells = filt_pos[!ind,]
  
  
  
  return(list(selected_cells,not_selected_cells, pos_fraction, increas_cut ))
}



################ main analysis function #################

select_cells = function(myFile,outDir_tmp,cond, myGene,do_shuffle = FALSE){
  
  # 3 means radius 20
  myRadius = 3
  
  all_data = read.csv(myFile)
  all_data = rotate_data(all_data)
  
  myShuffle= "real"
  if(do_shuffle){
    myShuffle = "SHUFFLE"
  }
  
  for(numRep in 1:3){
    data_rep = get_data_rep(all_data,myGene,cond,numRep,myRadius, do_shuffle = do_shuffle, outDir_tmp)
    selected_cells = data_rep[[1]]
    not_selected_cells  = data_rep[[2]]
    pos_fraction = data_rep[[3]]
    increas_cut = data_rep[[4]]
    
    radius_info = names(selected_cells)[myRadius]
    radius_info = gsub("X","",radius_info)
    radius_info = as.numeric(radius_info)
    radius_info = (radius_info/95) *10
    
    
    if(!do_shuffle){
      outFileData = paste0(outDir_tmp,"distance_plot_cells_rep_",numRep,"_",myGene, cond,"_",myShuffle,".csv")
      write.csv(selected_cells, outFileData, row.names = F)
      outFileData = paste0(outDir_tmp,"distance_plot_cells_rep_",numRep,"_",myGene, cond,"_",myShuffle,"NotSelected.csv")
      write.csv(not_selected_cells, outFileData, row.names = F)
      
      
      # print number selected
      num_sel = nrow(selected_cells)
      not_sel = nrow(not_selected_cells)
      total = num_sel+not_sel
      tmp_num_selected = data.frame(num_selected = num_sel, num_not = not_sel, 
                                    rep = numRep, gene= myGene, 
                                    fraction_selected = num_sel/total, fraction_not = not_sel/total)
      if(numRep ==1){
        final_num_selected = tmp_num_selected 
      }else{
        final_num_selected = rbind(final_num_selected,tmp_num_selected)
      }
    }
  }
  
  
  if(!do_shuffle){
    outFileNumSelected = paste0(outDir_tmp,"distance_numSelected_",myGene, cond,"_",myShuffle,".csv")
    write.csv(final_num_selected, outFileNumSelected, row.names = F)
  }

}
