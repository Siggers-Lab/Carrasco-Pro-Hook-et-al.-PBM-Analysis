library(tidyverse)
#library(magrittr)
library(dplyr)
library(gmodels)
library(pheatmap)
library(Biostrings)
library(philentropy)  #For the vector distances

##################
#  Filter the PBM data based on significance; it also enforces REPLICATE significance
filter_pbm_data <- function(df, ARG_zscore_cut, ARG_delta_cut, ARG_pval_cut){

  #Ensure these are numeric
  df$neg_log10_pval_fisher_ttest_3 <- as.numeric(as.character(df$neg_log10_pval_fisher_ttest_3))
  df$mean_REF_3 <- as.numeric(as.character(df$mean_REF_3))
  df$mean_SNP_3 <- as.numeric(as.character(df$mean_SNP_3))
  df$delta_3 <- as.numeric(as.character(df$delta_3))

  sig_conditions = which(
    df$neg_log10_pval_fisher_ttest_3 > ARG_pval_cut &
      (df$mean_REF_3 > ARG_zscore_cut | df$mean_SNP_3 > ARG_zscore_cut ) &
      ( df$delta_3 > ARG_delta_cut | df$delta_3 < -1*ARG_delta_cut)
  )

  dfs = droplevels(df[sig_conditions,])

  # REPLICATE FILTER - ensure that BOTH replicates are in the filtered matrix
  #
  # Each probe is identified by the columns: SNP_id_allele, PBM_exp, rep
  for(i in 1:nrow(dfs)){

    #Determine the ID and PBM experiment
    SNP_id_allele = dfs$SNP_id_allele[i]
    PBM_exp = dfs$PBM_exp[i]

    filter = which( dfs$SNP_id_allele == SNP_id_allele & dfs$PBM_exp == PBM_exp ) #Set ID and PBM filter criteria

    #If there are MORE than one row with this ID and PBM_exp - you know at least 2 replicates passed criteria
    #Set their sig_rep flag to 1
    td = dfs[filter,]
    if(nrow(td) > 1) { dfs$sig_rep[i] = 1}
    else{ dfs$sig_rep[i]=0}

  }

  return(dfs[which(dfs$sig_rep==1),]) #Get only thse that are sign for both replicates

}#end function
####################

##################
#  Filter the PBM data based on NO significance; it also enforces REPLICATE significance
filter_pbm_data_nodiff <- function(df, ARG_zscore_cut, ARG_delta_cut, ARG_pval_cut){

  #Ensure these are numeric
  df$neg_log10_pval_fisher_ttest_3 <- as.numeric(as.character(df$neg_log10_pval_fisher_ttest_3))
  df$mean_REF_3 <- as.numeric(as.character(df$mean_REF_3))
  df$mean_SNP_3 <- as.numeric(as.character(df$mean_SNP_3))
  df$delta_3 <- as.numeric(as.character(df$delta_3))

  #This is a filter for NOT significant
  sig_conditions = which(
    df$neg_log10_pval_fisher_ttest_3 < ARG_pval_cut &
      (df$mean_REF_3 > ARG_zscore_cut | df$mean_SNP_3 > ARG_zscore_cut ) &
      ( abs(df$delta_3) < abs(ARG_delta_cut))
  )

  dfs = droplevels(df[sig_conditions,])

  # REPLICATE FILTER - ensure that BOTH replicates are in the filtered matrix
  #
  # Each probe is identified by the columns: SNP_id_allele, PBM_exp, rep
  for(i in 1:nrow(dfs)){

    #Determine the ID and PBM experiment
    SNP_id_allele = dfs$SNP_id_allele[i]
    PBM_exp = dfs$PBM_exp[i]

    filter = which( dfs$SNP_id_allele == SNP_id_allele & dfs$PBM_exp == PBM_exp ) #Set ID and PBM filter criteria

    #If there are MORE than one row with this ID and PBM_exp - you know at least 2 replicates passed criteria
    #Set their sig_rep flag to 1
    td = dfs[filter,]
    if(nrow(td) > 1) { dfs$sig_rep[i] = 1}
    else{ dfs$sig_rep[i]=0}

  }

  return(dfs[which(dfs$sig_rep==1),]) #Get only thse that are sign for both replicates

}#end function
####################

######################
# Filter the MERGED data by PBM and MPRA, it also enforces REPLICATE significance
filter_pbm_mpra_data <- function(df, ARG_zscore_cut, ARG_delta_cut, ARG_pval_cut, ARG_skewFDR_cut){

  #Ensure these are numeric
  df$neg_log10_pval_fisher_ttest_3 <- as.numeric(as.character(df$neg_log10_pval_fisher_ttest_3))
  df$mean_REF_3 <- as.numeric(as.character(df$mean_REF_3))
  df$mean_SNP_3 <- as.numeric(as.character(df$mean_SNP_3))
  df$delta_3 <- as.numeric(as.character(df$delta_3))
  df$Skew.logFDR <- as.numeric(as.character(df$Skew.logFDR))

  sig_conditions = which(
    df$neg_log10_pval_fisher_ttest_3 > ARG_pval_cut &
      (df$mean_REF_3 > ARG_zscore_cut | df$mean_SNP_3 > ARG_zscore_cut ) &
      ( df$delta_3 > ARG_delta_cut | df$delta_3 < -1*ARG_delta_cut) &
      (df$Skew.logFDR > ARG_skewFDR_cut)
  )

  dfs = droplevels(df[sig_conditions,])

  # REPLICATE FILTER - ensure that BOTH replicates are in the filtered matrix
  #
  # Each probe is identified by the columns: SNP_id_allele, PBM_exp, rep
  for(i in 1:nrow(dfs)){

    #Determine the ID and PBM experiment
    SNP_id_allele = dfs$SNP_id_allele[i]
    PBM_exp = dfs$PBM_exp[i]

    filter = which( dfs$SNP_id_allele == SNP_id_allele & dfs$PBM_exp == PBM_exp ) #Set ID and PBM filter criteria

    #If there are MORE than one row with this ID and PBM_exp - you know at least 2 replicates passed criteria
    #Set their sig_rep flag to 1
    td = dfs[filter,]
    if(nrow(td) > 1) { dfs$sig_rep[i] = 1}
    else{ dfs$sig_rep[i]=0}

  }

  return(dfs[which(dfs$sig_rep==1),]) #Get only thse that are sign for both replicates

}#end function
#####################

######################
# Filter original MPRA datafile
filter_only_mpra_data <- function(df, ARG_skewFDR_cut){

  #Ensure these are numeric
  df$Skew.logFDR <- as.numeric(as.character(df$Skew.logFDR))

  sig_conditions = which(
    df$Skew.logFDR > ARG_skewFDR_cut
  )

  dfs = droplevels(df[sig_conditions,])

  return(dfs) #Get only thse that are sign for both replicates

}#end function
#####################



##################
#  Filter the PBM data based on significance; it also enforces REPLICATE significance
filter_pbm_data_zscore_only <- function(df, ARG_zscore_cut){

  #Ensure these are numeric
  df$mean_REF_3 <- as.numeric(as.character(df$mean_REF_3))
  df$mean_SNP_3 <- as.numeric(as.character(df$mean_SNP_3))

  sig_conditions = which(
      (df$mean_REF_3 > ARG_zscore_cut | df$mean_SNP_3 > ARG_zscore_cut )
  )

  dfs = droplevels(df[sig_conditions,])

  # REPLICATE FILTER - ensure that BOTH replicates are in the filtered matrix
  #
  # Each probe is identified by the columns: SNP_id_allele, PBM_exp, rep
  for(i in 1:nrow(dfs)){

    #Determine the ID and PBM experiment
    SNP_id_allele = dfs$SNP_id_allele[i]
    PBM_exp = dfs$PBM_exp[i]

    filter = which( dfs$SNP_id_allele == SNP_id_allele & dfs$PBM_exp == PBM_exp ) #Set ID and PBM filter criteria

    #If there are MORE than one row with this ID and PBM_exp - you know at least 2 replicates passed criteria
    #Set their sig_rep flag to 1
    td = dfs[filter,]
    if(nrow(td) > 1) { dfs$sig_rep[i] = 1}
    else{ dfs$sig_rep[i]=0}

  }

  return(dfs[which(dfs$sig_rep==1),]) #Get only thse that are sign for both replicates

}#end function
####################

#################
# Add COF support to each SNP_id_allele
#   For each SNP_id_allel  that is significant, determine how many experiments support this
#   Save this data in the $SNP_support column
summarize_COF_support = function(d){

  SNP_id_allele = unique(d$SNP_id_allele)

  # Each ROW for a unique SNP_id_allele is an indepdnent PBM experiment that supports it -- count them up
  for(i in SNP_id_allele){
    t = which(d$SNP_id_allele == i) #Get the rows that match this SNP_id_allele
    d[t,"SNP_support"] = length(t)/2 #Set the SNP_support column with the number that match (div by 2 so don't include reps)
    d[t,"SNP_support_code"]= ""

    #Determine if all the support is in the same 'direction'
    #Create a SNP_support_code string =  1 postivei delta, 0 negative delta
    code = vector()
    for(r in t){

      if (d[r,"delta_3"] < 0 ) { code = paste(code,"0",sep="") } else { code = paste(code,"1",sep="")}
    }

    d[t,"SNP_support_code"]= code
  }

  support = sort(unique(d$SNP_support))#Determine the different levels of support for SNP_alleles

  #Report the numbers for different levels of support
  sum = 0
  for(i in support) {
    t = which(d$SNP_support == i)
    tt = d[(d$SNP_support== i & d$cell == "SKMEL"),]
    print(paste("Num SNP_alleles with support",i," = ",length(t)/i/2,"  #SKMEL=",nrow(tt)/2)) #Divide by support # and 2 replicates
    sum = sum + length(t)/i/2
  }
  print(paste0("sum=",sum))
  #
  #######

  return(d)
} #end function
#################

#####################
# Read in Sebastian's SNP PWM predictions and TF family info
get_pwm_info = function(){

  dir = "/projectnb2/siggers/work/tsiggers/projects/2020_Cancer_SNVs/SNV_analysis/version2_redo_data_2021/"

    filename = paste0(dir,"predicted_drivers.tsv")
  d.pwm = read.table(file=filename,header=T,sep="\t") #This is from Sebastian's supplemental files

  #Read TF information (just the first few columns from Seabastian's tf_info.tsv)

  filename = paste0(dir,"tf_info_trimmed.txt")
  d.families = read.table(filename,header=T,sep="\t") #This is adapted from Sebastian's supplemental files tf_info.tsv

  #Process PWM data
  d.pwm = d.pwm %>%

    #Make a SNP_id_allele column (for joining with PBM data later)
    dplyr::mutate(SNP_id_allele=paste0("chr",chr,"-",start,"_",REF,"/",ALT)) %>%

    #Grab just some columns
    #Might consider adding some more (e.g., promoter/gene and the TF scores etc)
    dplyr::select(
      SNP_id_allele,
      promoter_name,
      promoter_id,
      tf_name,
      pwm_cisbp_id,
      skin_mut_freq,
      skin_significant
    ) %>%

    #Merge with the TF annotation data (allows us to see the TF families for each PWM)
    dplyr::inner_join(d.families,by= c("pwm_cisbp_id" = "Motif_ID")) %>%

    return()
} #end function
######################

####################
# Plot scores of the individual REF/SNP tiles probes (for all positions)
plot_REF_SNP_tiles = function(d,SNP_id_allele){


  xval =c(1,2,3,4,5,6)   #For plotting if you select that
  cv = c("black","red","blue")  #For plotting if you select that

    cvi = 1 %% length(cv) + 1  #For plotting (can drop, was used before)

    #Determine all experimentas for this SNP_id_allele
    pd = droplevels(d[which(d$SNP_id_allele == SNP_id_allele),])
#    message(paste(" PROBE=",probes[i]))

    for (j in 1:nrow(pd)){
      ## Visualize

      #    message(paste("      experiment=",pd[j,"PBM_exp"])," ",pd[j,"cell_line"],"  ",pd[j,"rep"],"    MEAN DELTA=",pd[j,"mean_REF_3"]-pd[j,"mean_SNP_3"])
      message(paste("      experiment=",pd$PBM_exp[j])," ",pd$cell_line[j],"  ",pd$rep[j],"    MEAN DELTA=",pd$mean_REF_3[j]-pd$mean_SNP_3[j])

      par(mfrow=c(1,3))
      yval = pd[j,c("REF_o1_.5","REF_o1_0","REF_o1_5","REF_o2_.5","REF_o2_0","REF_o2_5")]
      plot(xval,yval,ylim=c(-2,20),col=cv[cvi],main=paste(SNP_id_allele," ",pd[j,"PBM_exp"]," ",pd[j,"cell_line"],"  ",pd[j,"rep"],"   "))
      lines(xval,yval,col=cv[cvi])
      abline(v=3.5,col="red")
      abline(h=1,col="green")
      tt = paste("      ",probes[j], "REF")
      text(x=2,y=20,labels=tt)
      #    tt = paste("      ",pd[j,"PBM_exp"])
      tt = paste(pd$PBM_exp[j]," ",pd$cell_line[j],"  ",pd$rep[j])
      text(x=2,y=18,labels=tt)
      tt = paste("MEAN DELTA=",pd$mean_REF_3[j]-pd$mean_SNP_3[j])
      text(x=2,y=16,labels=tt)

      yval = pd[j,c("SNP_o1_.5","SNP_o1_0","SNP_o1_5","SNP_o2_.5","SNP_o2_0","SNP_o2_5")]
      plot(xval,yval,ylim=c(-2,20),col=cv[cvi])
      lines(xval,yval,col=cv[cvi])
      abline(v=3.5,col="red")
      abline(h=1,col="green")
      tt = paste(" SNP ")
      text(x=2,y=20,labels=tt)

      nval = pd[j,c("REF_o1_.5","REF_o1_0","REF_o1_5","REF_o2_.5","REF_o2_0","REF_o2_5")]
      dval = pd[j,c("SNP_o1_.5","SNP_o1_0","SNP_o1_5","SNP_o2_.5","SNP_o2_0","SNP_o2_5")]
      pval = abs(nval - dval)
      ymax = max(pval)
      plot(xval,pval,ylim=c(0,ymax+2),col=cv[cvi])
      lines(xval,pval,col=cv[cvi])
      abline(v=3.5,col="red")
      abline(h=2,col="green")
      # tt = paste(" REF-SNP")
      #  text(x=2,y=30,labels=tt)

      readline(paste("go"))

    } #end for j

}# end function

#####################

#####################
# Plot BARPLOT of Probes and TF/PWM prediction associations
make_barplot_SNPs_by_TF_family = function(dm,family_names){

  counts = vector()
  for(i in 1:length(family_names))
  {
    counts[i] = dm %>%
      #     filter(Family_Name == "Ets") %>%
      filter(Family_Name == family_names[i]) %>%
      select(SNP_id_allele) %>%  unique() %>% nrow()
  }

  #Make this as a barplot - or just calculate the number for all TFs
  db = dplyr::bind_cols(list(family_names,counts))
  names(db) = c("Family_Names","Counts")

  db =  dplyr::arrange(db,-Counts)

  # Basic barplot
  p<-ggplot(data=db, aes(x=reorder(Family_Names,Counts), y=Counts)) +
    geom_bar(stat="identity") +
    coord_flip()  #Make horizontal

  return(p)
} #end function
###################

##############
# Plot BARPLOT of COF support
#
# Categories =  UT_JURKAT_BRD4, UT_JURKAT_SRC1
# Frequence = Numbers of SNPs
#
make_barplot_COF_support = function(d){

  exp = unique(d$PBM_exp) #Get the Different PBM expers (e.g., UT_JURKAT_BRD4, etc)
  num = vector(mode="numeric") #Data Vector

   for(i in 1:length(exp)){
    xx = d[which(d$PBM_exp == exp[i]),] #Filter Data Frame for each PBM_exp
    num[i] = length(unique(xx$SNP_id_allele)) #Get Number SNP_id_alleles for this PBM_exp
   } #for i in 1:exp

  n = length(exp)+1
  exp[n] = "Aggregate"
  num[n] = length(unique(d$SNP_id_allele))

  ord = 1:length(num)
  Ad = data.frame(exp,num,ord)
  names(Ad) = c("PBM_exp","num","ord")
  Ads = Ad[order(-Ad$num),] #Re-order based on SUPPORT
  Ads$ord = 1:length(num)

  #    Ad = factor(PBM_exp,levels=Ads$PBM_exp)
   p = ggplot(data = Ads,aes(x=reorder(PBM_exp,-ord),y=num, fill=num)) +
      geom_bar(stat = "identity" ) +
     theme_minimal() + #Will make background white
      coord_flip()
  #p
  return(p)
} #end function
###############

##############
# Plot BARPLOT of COF support
# This will indicate both GAINS adn LOSSES
#
make_barplot_COF_support_sense = function(d){

  #Make a COF column with teh SRC1,TBL1XR1 etc
  d <- d %>%
    mutate(PBM_temp_exp = PBM_exp) %>%  #Duplicate column
    separate(PBM_temp_exp,into=c("stim","celltype","cof"),sep="_") %>%
    select(-stim,-celltype) #Drop these columns, not needed

  #
  #Extract all probes with a GAIN or a LOSS of binding
  d.loss <- d %>%
    filter(delta_3 < 0) %>%    #Filter out LOSSES
    mutate(cof_sense = paste0(cof,"_loss")) #Make new annotation column

  d.gain <- d %>%
    filter(delta_3 > 0) %>%    #Filter out LOSSES
    mutate(cof_sense = paste0(cof,"_gain")) #Make new annotation column

  d = rbind(d.loss,d.gain)

  exp = unique(d$cof_sense) #Get the Different COFs (e.g., BRD4_gain, BRD4_loss, .. )
  num = vector(mode="numeric") #Data Vector

  for(i in 1:length(exp)){
    xx = d[which(d$cof_sense == exp[i]),] #Filter Data Frame for each PBM_exp
    num[i] = length(unique(xx$SNP_id_allele)) #Get Number SNP_id_alleles for this PBM_exp
  } #for i in 1:exp

  n = length(exp)+1
  exp[n] = "Aggregate"
  num[n] = length(unique(d$SNP_id_allele))

  ord = 1:length(num)
  Ad = data.frame(exp,num,ord)
  names(Ad) = c("PBM_exp","num","ord")
  Ads = Ad[order(-Ad$num),] #Re-order based on SUPPORT
  Ads$ord = 1:length(num)

  #    Ad = factor(PBM_exp,levels=Ads$PBM_exp)
#  p = ggplot(data = Ads,aes(x=reorder(PBM_exp,-ord),y=num, fill=num)) +
   p = ggplot(data = Ads,aes(x=reorder(PBM_exp,-ord),y=num)) +
    geom_bar(stat = "identity" ) +
    theme_minimal() + #Will make background white
    coord_flip()
  #p
  return(p)
} #end function
###############


##############
# Plot BARPLOT of COF support
#
# Categories =  Support Number (Strength)  [1,2,3,4,5,6,7]
# Frequence = Numbers of SNPs
#
make_barplot_COF_support_v2 = function(d){

  support = sort(unique(d$SNP_support)) #Get the Different COF Support Numbers
  num = vector(mode="numeric") #Data Vector

  for(i in 1:length(support)){
    xx = d[which(d$SNP_support == support[i]),] #Filter Data Frame for each PBM_exp
    num[support[i]] = length(unique(xx$SNP_id_allele)) #Get Number SNP_id_alleles for this PBM_exp
    print(paste0(" i=",i," support #=",support[i],"  num=",num[support[i]]))
  } #for i in 1:exp


  ord = 1:length(num)
  Ad = data.frame(support,num,ord)
  names(Ad) = c("COF_support","num","ord")
  Ads = Ad[order(-Ad$COF_support),] #Re-order based on SUPPORT
  Ads$ord = 1:length(num)

  #    Ad = factor(PBM_exp,levels=Ads$PBM_exp)
  p = ggplot(data = Ads,aes(x=reorder(COF_support,-ord),y=num)) +
    geom_bar(stat = "identity" ) +
    theme_minimal()
  #p
  return(p)
} #end function
###############


##############
# Plot BARPLOT of COF support
make_heatmap_COF_support = function(d,filename){

  #Grab the COF part of the PBM_experiment name
  extract_cof_info = function(s){
    #PBM exp anme : UT_JURKAT_SRC4 -- this grabs the SRC4 part
     t = strsplit(as.character(s),"_")
     return(t[[1]][3])
  }

  #Make a COF column with teh SRC1,TBL1XR1 etc
  d = d %>%
      dplyr::mutate(cof = map_chr(PBM_exp,extract_cof_info))

  cell_cof_table = gmodels::CrossTable(d$cell,d$cof,prop.t=TRUE,prop.r=TRUE,prop.c=TRUE)

  write.table(cell_cof_table$t,file=filename,row.names=TRUE, col.names=TRUE,sep="\t",quote=FALSE) #Write data to file

#  library(gplots)
#  heatmap.2(cell_cof_table$t,scale="none")
#  heatmap.2(cell_cof_table$t, scale = "none", col = bluered(100),
#            trace = "none", density.info = "none")

  # For some details
  # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/#data-preparation
  p = pheatmap(cell_cof_table$t, cutree_rows = 1.5)
  return(p)
} #end function
###############

##############
# Plot BARPLOT of COF support
#
# This version separates GAIN and LOSS of binding
#
make_heatmap_COF_support_sense = function(d,filename){


  #Grab the COF part of the PBM_experiment name
  extract_cof_info = function(s){
    #PBM exp anme : UT_JURKAT_SRC4 -- this grabs the SRC4 part
    t = strsplit(as.character(s),"_")
    return(t[[1]][3])
  }

  #Make a COF column with teh SRC1,TBL1XR1 etc
  d = d %>%
    dplyr::mutate(cof = map_chr(PBM_exp,extract_cof_info))

  #Extract all probes with a GAIN or a LOSS of binding
  d.loss <- d %>%
    filter(delta_3 < 0) %>%    #Filter out LOSSES
    mutate(cof_sense = "loss") #Make new annotation column

  d.gain <- d %>%
    filter(delta_3 > 0) %>%    #Filter out LOSSES
    mutate(cof_sense = "gain") #Make new annotation column

  d.sense = rbind(d.loss,d.gain)

  cell_cof_table = gmodels::CrossTable(d.sense$cof,d.sense$cof_sense,prop.t=TRUE,prop.r=TRUE,prop.c=TRUE)

  write.table(cell_cof_table$t,file=filename,row.names=TRUE, col.names=TRUE,sep="\t",quote=FALSE) #Write data to file

  #  library(gplots)
  #  heatmap.2(cell_cof_table$t,scale="none")
  #  heatmap.2(cell_cof_table$t, scale = "none", col = bluered(100),
  #            trace = "none", density.info = "none")

  # For some details
  # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/#data-preparation
  p = pheatmap(cell_cof_table$t, cutree_rows = 1.5)
  return(p)
} #end function
###############

##############
# Plot BARPLOT of COF support for Diff and NoDiff probes
#
# This version separates GAIN and LOSS of binding, AND probes with NODIFF / NOCHANGE
#
make_heatmap_COF_support_sense_nodiff = function(d,d.nodiff,filename){



  #Grab the COF part of the PBM_experiment name
  extract_cof_info = function(s){
    #PBM exp anme : UT_JURKAT_SRC4 -- this grabs the SRC4 part
    t = strsplit(as.character(s),"_")
    return(t[[1]][3])
  }

  #Make a COF column with teh SRC1,TBL1XR1 etc
  d <- d %>%
    dplyr::mutate(cof = map_chr(PBM_exp,extract_cof_info))

  d.nodiff <- d.nodiff %>%
    dplyr::mutate(cof = map_chr(PBM_exp,extract_cof_info))

  #Extract all probes with a GAIN or a LOSS of binding
  d.loss <- d %>%
    filter(delta_3 < 0) %>%    #Filter out LOSSES
    mutate(cof_sense = "loss") %>% #Make new annotation column
    select(SNP_id_allele,cof,cof_sense) #Select just these cols, easier to bind with nodiff data below

  d.gain <- d %>%
    filter(delta_3 > 0) %>%    #Filter out LOSSES
    mutate(cof_sense = "gain") %>% #Make new annotation column
    select(SNP_id_allele,cof,cof_sense)

  d.nodiff <- d.nodiff %>%
    mutate(cof_sense = "no_change") %>%
    select(SNP_id_allele,cof,cof_sense)

  d.sense = rbind(d.loss,d.gain,d.nodiff)

  cell_cof_table = gmodels::CrossTable(d.sense$cof,d.sense$cof_sense,prop.t=TRUE,prop.r=TRUE,prop.c=TRUE)

  write.table(cell_cof_table$t,file=filename,row.names=TRUE, col.names=TRUE,sep="\t",quote=FALSE) #Write data to file

  #  library(gplots)
  #  heatmap.2(cell_cof_table$t,scale="none")
  #  heatmap.2(cell_cof_table$t, scale = "none", col = bluered(100),
  #            trace = "none", density.info = "none")

  # For some details
  # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/#data-preparation
  p = pheatmap(cell_cof_table$t, cutree_rows = 1.5)
  return(p)
} #end function
###############

############
# Plot BARPLOT of Probe Type Enrichment
# d = data frame
# df = Filtered data frame with only significant probes
make_barplot_analysis_by_probe_type = function(d,df,filename){

  types = unique(d$analysis_type)

  types = types[!grepl("literature_genes",types)] #Remove these probes (01/20/22) - decided to drop for paper

  total_values = vector()
  sig_values = vector()
  for(i in 1:length(types)){

    dt  = d %>%  filter(analysis_type == types[i])
    total = length(unique(dt$SNP_id_allele))  #Number of SNP_id_alleles in Total
    total_values = append(total_values,total)

    dft = df %>%  filter(analysis_type == types[i])
    sig = length(unique(dft$SNP_id_allele))   #Number of Significant SNP_id_alleles that passed cutoff
    sig_values = append(sig_values,sig)

    #   print(paste("i = ",i," type=",types[i]," total=",total,"  sig=",sig,"  ratio=",sig/total))
    print(paste("i = ",i," type=",types[i]," total=",total,"  sig=",sig))
  }

  v = as.data.frame(cbind(types,total_values,sig_values))
  colnames(v) = c("analysis_type","total","sig")
  v$total = as.numeric(v$total)
  v$sig = as.numeric(v$sig)

  v = v %>% dplyr::mutate(not_sig = total - sig)
  v = v %>% mutate(ratio = sig/total)
  v =  dplyr::arrange(v,-ratio) #Order the data frame based on 'ratio'

  write.table(v,file=filename,row.names=FALSE, col.names=TRUE,sep="\t") #Write data to file

  # Basic barplot
  p<-ggplot(data=v, aes(x=reorder(analysis_type,ratio), y=ratio, fill=ratio)) +
    #   geom_bar(stat="identity", color="blue",fill="steelblue") +
    #   geom_bar(stat="identity",fill="steelblue",width=0.9) +
    geom_bar(stat="identity",width=0.9) +
    theme_minimal() + #Will make background white
    coord_flip()  #Make horizontal

  contingency_table = as.matrix(v[,c("sig","not_sig")])
  x = chisq.test(x=contingency_table,correct=FALSE)
  print(paste0(" X-squared=",x$statistic," P-value < ",x$p.value))

  return(p)
}#end function
#
###########

############
# Plot BARPLOT of SNP register

make_barplot_SNP_registers = function(d.temp.seq){

  register_vals = unique(na.omit(as.numeric(d.temp.seq$register)))
  registers = c(-10:10)
  total_values = vector()

  for(i in 1:length(registers)){
    #Get all SNPs that match the target register
    dt  = d.temp.seq %>%  filter(register == registers[i])
    total_values[i] = length(unique(dt$SNP_id_allele))  #Number of SNP_id_alleles in Total
  }

  v = cbind.data.frame(register=registers,num=total_values)

#  v =  dplyr::arrange(v,-ratio) #Order the data frame based on 'ratio'

  # Basic barplot
#  p<-ggplot(data=v, aes(x=register, y=num, fill=num)) +
   p<-ggplot(data=v, aes(x=register, y=num)) +
    #   geom_bar(stat="identity", color="blue",fill="steelblue") +
    #   geom_bar(stat="identity",fill="steelblue",width=0.9) +
    geom_bar(stat="identity",width=0.9) +
    theme_minimal()
#  + #Will make background white
#    coord_flip()  #Make horizontal

  return(p)
}
#####

##########
# Plot Volcano Scatter Plot
make_volcano_scatter = function(d,d.filter,title){

  #  dx = "mean_REF_3"
  dx = "delta_3"
  dy = "neg_log10_pval_fisher_ttest_3"

  p = ggplot(data=d,aes_string(x=dx,y=dy)) +
    geom_point(col="grey",size=2) +
#    geom_smooth( method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
    geom_point(data=d.filter,aes_string(x=dx,y=dy),col="firebrick2", size=2) +
    labs(
      title = title,
      x = "Delta z-score",
      y = "Negative log10(p-value)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5,vjust=-1,size = 15),
      axis.title.x = element_text(size = 16),
      axis.text.x =  element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y =  element_text(size = 14))

  return (p)
}
#######

##########
# Plot Pairwise Scatter Plot (used for replicate comparisons)
make_replicate_scatter = function(d,title){


  #Make the data from REP1/REP2 as separatre columns to plot
  d.wide  = d %>%
    select(SNP_id_allele,PBM_exp,mean_REF_3,rep) %>%
    pivot_wider(
      names_from = "rep",
      values_from = c("mean_REF_3")
    )

  dx = "REP1"
  dy = "REP2"

#  cor = cor(d.wide$REP1,d.wide$REP2,method="spearman")
  cor = cor(d.wide$REP1,d.wide$REP2,method="pearson")
  cor = round(cor, digits=2)

  p = ggplot(data=d.wide,aes_string(x=dx,y=dy)) +
    geom_point(col="grey",size=2) +
    xlim(-3,15)+
    ylim(-3,15)+
#    geom_smooth( method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
#    geom_smooth(method='lm',formula=d.wide$REP2~d.wide$REP1) +
    geom_smooth(method='lm') +
    labs(
      title = title,
      x = "z-score",
      y = "z-score"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5,vjust=-1,size = 15),
      axis.title.x = element_text(size = 16),
      axis.text.x =  element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.y =  element_text(size = 14)) +
    annotate("text", x = 0, y = 12, size=6, label = paste0("r = ",cor))

  return (p)
}
#######


#########
#
# Get Sequence Data for all SNP probes (i.e., all SNP_id_alleles)
# (04/22) Also Get the TRIPLETs around the SNPs!! Used to look at mutation bias

get_SNP_sequence_data <- function(da){


  #GET all the SNP sequences
  da.SNP <- da %>%
    mutate(SNP_id_allele = paste0(da$SNP_id,"_",da$SNP_allele)) %>%
    #            filter(SNP_id_allele == "chr11-47448149_G/T") %>%
    #           filter(SNP_id == "chr11-47448149") %>%
    filter(probe_type == "SNP") %>%
    filter(probe_or == "o1") %>%
    filter(probe_repl == "r1") %>%
    filter(tile_order == "0") %>%
    select(all_of(c("SNP_id_allele","SNP_id","SNP_allele","SNV_pos_offset","target_seq","probe_type","analysis_type"))) %>%
    dplyr::rename(SNP_seq = target_seq)

  #GET all the REF sequence
  da.REF <- da %>%
    #      filter(SNP_id == "chr11-47448149") %>%
    filter(probe_type == "REF") %>%
    filter(probe_or == "o1") %>%
    filter(probe_repl == "r1") %>%
    filter(tile_order == "0") %>%
    select(all_of(c("SNP_id","target_seq"))) %>%
    dplyr::rename(REF_seq = target_seq)

  #Join so each line has one SNP (with sequence) and its REF sequence  (AND GET TRIPLETS!!)
  da.seq= dplyr::inner_join(da.SNP,da.REF,by='SNP_id') %>%
    mutate(SNP_triplet = substring(SNP_seq,SNV_pos_offset-1,SNV_pos_offset+1)) %>%  #Get Triplet
    mutate(REF_triplet = substring(REF_seq,SNV_pos_offset-1,SNV_pos_offset+1))      #Get Triplet

  return(da.seq)

}
#####



####
# Select specific MOTIFS from JASPAR/MEME list (this is a universalmotif object)
get_motifs_from_universalmotif_list <- function(db,name){

  #DEFINE FUNCTION
  select_motif = function(motif,motif_name){
  #       if(motif@altname == motif_name) { return(motif) }
  if(grepl(motif_name,motif@altname)) { return(motif)}
}

  m = db %>%
  map(select_motif,name) %>%
  discard(is.null)

  return(m)
}
###

#####
# Identify Differential binding of motifs in REF_seq and SNP_seq
#
# Return data frame with match data
identify_differential_binding <-function(motif,df,index){

#TEMP HARDCODE for DEBUG
#  REF.seq = d.temp.seq$REF_seq[i]
#  SNP.seq = d.temp.seq$SNP_seq[i]
#  SNP_id_allele = d.temp.seq$SNP_id_allele[i]
#  motif = m[[m.index]]
#  motif = m

   REF.seq = df$REF_seq[index]
   SNP.seq = df$SNP_seq[index]
   SNP_id_allele = df$SNP_id_allele[index]

   thresh = 0.005
   thresh.type = "pvalue"

   #Define DNAStrings
   REF.dna = DNAString(REF.seq)
   REF.dna.rc = reverseComplement(REF.dna)
   SNP.dna = DNAString(SNP.seq)
   SNP.dna.rc = reverseComplement(SNP.dna)

   #Scan DNA strings with motif
   results.REF.dna <- scan_sequences(motif,REF.dna,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
   if(nrow(results.REF.dna) > 0){
     results.REF.dna <- results.REF.dna %>%
       mutate(type = "REF") %>%
       mutate(tag = paste0(start,":+")) %>%
       mutate(strand = "+")
   }

   results.SNP.dna <- scan_sequences(motif,SNP.dna,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
   if(nrow(results.SNP.dna) > 0){
     results.SNP.dna <- results.SNP.dna %>%
       mutate(type = "SNP") %>%
       mutate(tag = paste0(start,":+")) %>%
       mutate(strand = "+")
   }

   results.REF.dna.rc = scan_sequences(motif,REF.dna.rc,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
   if(nrow(results.REF.dna.rc) > 0){
     results.REF.dna.rc <- results.REF.dna.rc %>%
       mutate(type = "REF") %>%
       mutate(tag= paste0(start,":-")) %>%
       mutate(strand = "-")
   }

   results.SNP.dna.rc = scan_sequences(motif,SNP.dna.rc,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
   if(nrow(results.SNP.dna.rc) > 0){
     results.SNP.dna.rc <- results.SNP.dna.rc %>%
       mutate(type = "SNP") %>%
       mutate(tag= paste0(start,":-")) %>%
       mutate(strand = "-")
   }

   results <- rbind(results.REF.dna,results.SNP.dna,results.REF.dna.rc,results.SNP.dna.rc)

   if(nrow(results)==0) {return(results)}

   results <- as_tibble(results) %>%
              mutate(processed = paste0("false")) %>%
              mutate(tag2 = paste0(motif,":",type,":",tag))


   SNP_id_alleles = vector()
   motif_names = vector()
   start_vals = vector()
   REF_scores = vector()
   SNP_scores = vector()
   strand_vals = vector()

   # Go through each line and determine differential matches
   for(i in 1:nrow(results)){

      if(results$processed[i] == "true"){next} #Skip if already processed

      if(results$type[i] == "REF"){

        find_tag = paste0(results$motif[i],":SNP:",results$tag[i]) #Find the data for the SNP at this position
        j = which(results$tag2 == find_tag) #Find match

        SNP_score = 0 #Set default
        if(length(j) > 0) { SNP_score = results$score[j] } #If match found
        results$processed[j] = "true"

        SNP_id_alleles = append(SNP_id_alleles,SNP_id_allele) #Get SNP_id_allele
        motif_names = append(motif_names,results$motif[i])    #Get Motif name
        start_vals = append(start_vals,results$start[i])  #Get match START
        REF_scores = append(REF_scores,results$score[i])  #Get REF score
        SNP_scores = append(SNP_scores,SNP_score)         #Get SNP score
        strand_vals = append(strand_vals,results$strand[i])
      }else{

        find_tag = paste0(results$motif[i],":REF:",results$tag[i]) #Find the data for the SNP at this position
        j = which(results$tag2 == find_tag) #Find match

        REF_score = 0 #Set default
        if(length(j) > 0) { REF_score = results$score[j] } #If match found
        results$processed[j] = "true"

        SNP_id_alleles = append(SNP_id_alleles,SNP_id_allele) #Get SNP_id_allele
        motif_names = append(motif_names,results$motif[i])    #Get Motif name
        start_vals = append(start_vals,results$start[i])  #Get match START
        REF_scores = append(REF_scores,REF_score)         #Get REF score
        SNP_scores = append(SNP_scores,results$score[i])  #Get SNP score
        strand_vals = append(strand_vals,results$strand[i])
      }
   }

   results.summary = cbind(SNP_id_allele=SNP_id_alleles,motif=motif_names,start=start_vals,REF_scores,SNP_scores,strand=strand_vals) %>%
                     as.data.frame() %>%
                     dplyr::mutate(diff_scores = as.numeric(REF_scores) - as.numeric(SNP_scores))

   #%>%
  #                   filter(diff_scores != 0)   #Remove when REF and SNP matches are the same (not differential)

   return(results.summary)
}
#
####

###
# Determine ETS register for the matches
determine_ETS_register <-function(results,m.data){

  #I have select both ETS motifs (ETV3, ETV6) so that they start at -1 ETS position
  #
  #     1   5    0   4               SNP is position 14 in + sense
  #REF "CCTTGCGTCATTTCCTGTAGTGTGCT
  #SNP "CCTTGCGTCATTTTCTGTAGTGTGCT
  #                  3  0    5   1   SNP is position 13 in - sense
  #Calculate the ETS position


  #Set Register for matches in FORWARD strand
  results.pos = results %>% filter(strand == "+")
    if(nrow(results.pos)>0){
      results.pos = results.pos %>% mutate(SNP_pos = 14)
     }

  #Set Register for matches in REVERSE COMPLEMENT strand
  results.neg = results %>% filter(strand == "-")
   if(nrow(results.neg)>0){
     results.neg = results.neg %>% mutate(SNP_pos = 13)
   }

  #Recombine results
  results <- rbind(results.pos,results.neg)

  if(nrow(results)==0) { return(results)} #If no matches - return empty frame


  #   14 - START.MATCH + START.REG
  #   14 - 14 + -2 = -2
  #   14 - 13 + -2 = -1
  #   14 - 9  +  0 = 5
  #   14 - 13 +  0 = 1
  #
  #   9    4
  #   CGGAATC    14 - 9 = 5 (reg)
  #        *
  #  13
  #   CGGAATC    14 -  13  = 1
  #    *
  #   4
  #   CGGAA      14 - 14 = 0 (reg)
  #  -1  1  2

  #Calculate the REGISTER values
  results <- dplyr::inner_join(results,m.data,by='motif')

  results <- results %>%
          mutate(register = as.numeric(SNP_pos)-as.numeric(start)+as.numeric(start_reg) )


  #Determine BEST register
  #
  # Issues to consider
  # 1. There might be multiple matches, so I need to identify which is best
  # 2. Generally it might be the match where REF-SNP is maximal (FIRST GUESS)
  # 3. But if I use a shortened ETS motif sometime is MISS the SNP, in which case
  #    the alternative might be the BEST HIT
  # 4. So, do these TWO filters,
  #  (1) BEST DIFFERENTIAL  -- this one isn't optimal. It biases towards weak REF hits were SNP is absent, making differential high
  #  (1) DIFFERENTIAL AND BEST SCORE -- this is more optimal filter
  # , (2) BEST HIT THAT DOESN"T CHANGE

  #Get sites where REF and SNP are different
  results.diff <- results %>%
          filter(diff_scores != 0)

    if(nrow(results.diff)>0){
        #Identify Position where the SCORE (REF or SNP) is maximal
        if(max(results.diff$REF_scores) > max(results.diff$SNP_score)){
              max.diff <- results.diff[which.max(results.diff$REF_scores),]
        }else{
              max.diff <- results.diff[which.max(results.diff$SNP_scores),]
        }
    }else { max.diff = results.diff} #Need to give it something to return

  #Get sites where REF and SNP are the same
  results.nodiff <- results %>%
          filter(diff_scores == 0)

    #Identify Position where the SCORE is maximal
    max.nodiff <- results.nodiff[which.max(results.nodiff$REF_scores),]

  #
  if(nrow(max.diff) >0) { return(max.diff)} #Preferentially return this one
  if(nrow(max.nodiff)>0) { return(max.nodiff)}

}

#
####


#####
# Identify Differential binding of motifs in REF_seq and SNP_seq
#
# Return data frame with match data
report_motif_binding_details <-function(motif,df,index){

  REF.seq = df$REF_seq[index]
  SNP.seq = df$SNP_seq[index]
  SNP_id_allele = df$SNP_id_allele[index]

  thresh = 0.005
  thresh.type = "pvalue"

  #Define DNAStrings
  REF.dna = DNAString(REF.seq)
  REF.dna.rc = reverseComplement(REF.dna)
  SNP.dna = DNAString(SNP.seq)
  SNP.dna.rc = reverseComplement(SNP.dna)

  #Scan DNA strings with motif
  results.REF.dna <- scan_sequences(motif,REF.dna,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
  if(nrow(results.REF.dna) > 0){
    results.REF.dna <- results.REF.dna %>%
      mutate(type = "REF") %>%
      mutate(tag = paste0(start,":+")) %>%
      mutate(strand = "+")
  }

  results.SNP.dna <- scan_sequences(motif,SNP.dna,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
  if(nrow(results.SNP.dna) > 0){
    results.SNP.dna <- results.SNP.dna %>%
      mutate(type = "SNP") %>%
      mutate(tag = paste0(start,":+")) %>%
      mutate(strand = "+")
  }

  results.REF.dna.rc = scan_sequences(motif,REF.dna.rc,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
  if(nrow(results.REF.dna.rc) > 0){
    results.REF.dna.rc <- results.REF.dna.rc %>%
      mutate(type = "REF") %>%
      mutate(tag= paste0(start,":-")) %>%
      mutate(strand = "-")
  }

  results.SNP.dna.rc = scan_sequences(motif,SNP.dna.rc,threshold = thresh, threshold.type = thresh.type) %>% as_tibble()
  if(nrow(results.SNP.dna.rc) > 0){
    results.SNP.dna.rc <- results.SNP.dna.rc %>%
      mutate(type = "SNP") %>%
      mutate(tag= paste0(start,":-")) %>%
      mutate(strand = "-")
  }

  results <- rbind(results.REF.dna,results.SNP.dna,results.REF.dna.rc,results.SNP.dna.rc)

  if(nrow(results)==0) {return(results)}

  results <- as_tibble(results) %>%
    mutate(processed = paste0("false")) %>%
    mutate(tag2 = paste0(motif,":",type,":",tag))

  return(results)
}
#
####


##############
# FILTER corec_motifs LIST for particular SNP
match_seed_in_corec_list <-
  function(
      corec_motifs,
      pattern = "NA"
             ) {

    match_function <- function(x,pattern){
      str_detect(x@seed_name,pattern)
    }

     return(
       corec_motifs %>%
         purrr::keep(match_function,pattern=pattern)
     )

     #Other way to code same thing with a loop
#    m = vector(mode="list") #Initialize a list for Matches
#     for(i in 1:length(corec_motifs)){
#      if(str_detect(corec_motifs[[i]]@seed_name,pattern)){
#       m = c(m,corec_motifs[[i]])
#      } #end if
#    } #end for
#    return(m)


  }#end function

#############
# MATCH pbm condition in corec_list
#
match_pbm_cond_in_corec_list <-
  function(
    corec_motifs,
    pattern = "NA"
  ) {

    match_function <- function(x,pattern){
      str_detect(x@pbm_condition,pattern)
    }

    return(
      corec_motifs %>%
        purrr::keep(match_function,pattern=pattern)
    )

    #Other way to code same thing with a loop
    #    m = vector(mode="list") #Initialize a list for Matches
    #     for(i in 1:length(corec_motifs)){
    #      if(str_detect(corec_motifs[[i]]@seed_name,pattern)){
    #       m = c(m,corec_motifs[[i]])
    #      } #end if
    #    } #end for
    #    return(m)


  }#end function

##########
#
# helper function to perform the reverse complentation of the sequences
get.rev.compl <- function(DNA_seq) {
  # returns the reverse complement of the input DNA sequence
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  temp_seq <- chartr('ACGTN', 'TGCAN', temp_seq)
  return(temp_seq)
}

#######
#
# This will determine the MUTATION SIGNATURE based on the COSMIC signature syntax
# https://cancer.sanger.ac.uk/signatures/sbs/
#
# This has all mutations in 6 groups
# C->A, C->G, C->T and T->A, T->G and T->C
# And then all triplets with the mutation in the middle.
#
make_mutation_signature <- function(seed_nuc,SNV_nuc,REF_triplet,SNP_triplet){

  # If seed_nuc = C or T, we are in the correct 'sense'
  # If not, all thngs need to be switched to reverse complement
  #
  if ( (seed_nuc == "G") | (seed_nuc == "A") ){
    seed_nuc = get.rev.compl(seed_nuc)
    SNV_nuc = get.rev.compl(SNV_nuc)
    REF_triplet = get.rev.compl(REF_triplet)
    SNP_triplet = get.rev.compl(SNP_triplet)
  }

  signature_tag = paste0(REF_triplet,"_",seed_nuc,"_",SNV_nuc)

  return(signature_tag)
}


###
#
# Get SNP comparison
#
# Massage data into a matrix witih
#                 UT_JURKAT_BRD4   UT_JURKAT_SCR1  etc etc.
# SNP_id_allele       4.5             13.4
#
# This format can be used to COMPARE/CORRELATE SNP similarity

get_SNP_comparison_format <- function(dp){


  cols_constant = c("SNP_id_allele","cell_line","PBM_exp")

  cols_to_average = c("mean_REF_3","mean_SNP_3","delta_3",
                      "rep")

  #Get all the data, and make a unique SNP x PBM_exp ID for binding later
  d = dp %>%
    filter( analysis_type == "cancer_snv_sig") %>%    #Select just the 2557 Cancer_snv_sig
    select(all_of(c(cols_constant,cols_to_average)))  %>%
    mutate(SNP_id_allele_exp = paste0(SNP_id_allele,"_",PBM_exp))  #Create a SNP_id_allele x PBM_exp ID

  # Get all the constant data that doesn't change across PBM replicates
  d_constant = d %>%
    select( all_of(c("SNP_id_allele_exp",cols_constant)) )  %>%
    distinct()  #Remove the duplicates that were present for all the rep1, rep2 for the data


  #Get average data for all the PBM columns
  # (this data frame is just these averaged values and the SNP_id_allele_exp column)
  d_averaged = d %>%
    dplyr::group_by(SNP_id_allele_exp) %>%
    dplyr::summarise(avg_mean_REF_3 = mean(mean_REF_3),
                     avg_mean_SNP_3 = mean(mean_SNP_3),
                     avg_delta_3 = mean(delta_3))


  #Join the constant and replicate averaged PBM data scores together.
  d.join = inner_join(d_constant,d_averaged,by="SNP_id_allele_exp") %>%
    filter(cell_line == "jurkat")

  #Pivot to Longer Format for PBM variables for Plotting
  d.wide  = d.join %>%
    #   select(!c("SNP_id_allele_exp")) %>%   #Remove this column - it messes up the pivot
    #   select( all_of(c("SNP_id_allele","PBM_exp","avg_mean_REF_3","avg_mean_SNP_3","avg_delta_3")) ) %>%
    #   select(all_of(c("SNP_id_allele","PBM_exp","avg_mean_REF_3"))) %>%
    select(all_of(c("SNP_id_allele","PBM_exp","avg_delta_3"))) %>%
    pivot_wider(
      names_from = "PBM_exp",
      #      values_from = c("avg_mean_REF_3","avg_mean_SNP_3","avg_delta_3")
      #        values_from = c("avg_mean_REF_3")
      values_from = c("avg_delta_3")
    ) %>%
    dplyr::mutate(ID=row_number())  #Add ID, will be used to join data below

  d.wide[is.na(d.wide)] <- 0

  #Get some meta data, will add back to the UMAP dataset

  d.wide.select = d.wide %>%
    #    select(SNP_id_allele,ID,UT_JURKAT_BRD4)
    #    select(SNP_id_allele,ID,UT_JURKAT_NCOR)
    select(SNP_id_allele,ID,
           UT_JURKAT_BRD4,
           UT_JURKAT_TBL1XR1,
           UT_JURKAT_NCOR,
           UT_JURKAT_SRC1,
           UT_JURKAT_MOF,
           UT_JURKAT_RBBP5)


  return(d.wide.select)
}

####
#
#

return_SNP_id_allele_to_Promoters <- function(d){

   #
   # Initialize data structure to return
   dr = d %>%
     select(SNP_id_allele,SNV_pos) %>%
     dplyr::mutate( promoter_name = "temp") %>%
     unique()

   for(i in 1:length(dr$SNP_id_allele)){
     #
     # Get all promoter information associated with each SNP_id_allele
     #
     dt = d %>%
       filter(SNP_id_allele == dr$SNP_id_allele[i]) %>%
       select(promoter_name) %>%
       unique()

     #
     # Make concatenated promoter name
     dr$promoter_name[i] = paste0(dt$promoter_name,collapse=":")
   }

   return(dr)
}

#######
#
# PAIRWISE - SNP COMPARISON.
#
# SNP values
#
compare_snps <- function(d,si,sj){


#d = dp.vals

  di = d %>%
 #   filter( SNP_id_allele %in% c(si,sj) )
    filter (SNP_id_allele == si) %>%
    select(3:8) %>%
    as.vector %>%
    as.numeric

  dj = d %>%
    filter( SNP_id_allele == sj ) %>%
    select(3:8) %>%
    as.vector() %>%
    as.numeric()

   vals = vector()
#  c = cor(x, y, method = c("pearson", "kendall", "spearman"))
  vals[1] = cor(di,dj, method = c("pearson"))
  vals[2] = cor(di,dj, method = c("spearman"))
  vals[3] = sqrt(sum((di - dj)^2))   #Euclidean distance

  return(vals)
}


