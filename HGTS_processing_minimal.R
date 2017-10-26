#! /usr/bin/Rscript

options(stringsAsFactors = F)
options(java.parameters = "- Xmx1024m")

# HTGS abstar processing - minimal input

# input_file="EK-gDNA-053017_S22_mergedpairedend_15overlap.csv"
# input_type="gDNA"
# dedup=T

library(optparse)

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, 
              help="Abstar output file name", metavar="character"),
  make_option(c("-e", "--exact_matches"), type="character", default=NULL, 
              help="Exact matches (gDNA) file name", metavar="character"),
  make_option(c("-c", "--cross_priming"), type="character", default=NULL, 
              help="Cross priming seqs (gDNA) file name", metavar="character"),
  make_option(c("-t", "--input_type"), type="character", default=NULL, 
              help="input type (gDNA or RNA)", metavar="character"),
  make_option(c("-d", "--dedup"), type="logical", default=NULL, 
              help="deduplication", metavar="character")),
  make_option(c("-l", "--lower_length"), type="integer", default=150, 
            help="lower cutoff for V gene length (bp)", metavar="number")),
  make_option(c("-u", "--upper_length"), type="integer", default=311, 
            help="upper cutoff for V gene length (bp)", metavar="number"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(stringr)
library(Biostrings)

HGTS_processing<-function(input_file,exact_matches,cross_priming,input_type,dedup,vlower,vupper){

  logfile<-file(paste0(sub(".csv","",input_file),"_",format(Sys.time(), format = "%Y-%m-%d-%H%M"),"_dedup_",dedup,"_report.txt"),open="a")
  
  cat(paste0("input file: ",input_file),file = logfile, sep="\n")
  cat(paste0("input type: ",input_type),file = logfile, sep="\n")
  cat(paste0("deduplication: ",dedup),file = logfile, sep="\n")
  cat("-----",file = logfile, sep="\n")
  
#infile<-read.csv("EK-gDNA-053017_S22_mergedpairedend_15overlap_sample.csv")  
infile<-read.csv(input_file)

cat(paste0("Input file read successfully. Number of reads: ",nrow(infile)),file = logfile, sep="\n")
cat("-----",file = logfile, sep="\n")

# filtering common to both gDNA and RNA

# filter I - kappa chain only
chaincol<-grep("chain",colnames(infile))
infile_I<-infile[which(infile[,chaincol]=="kappa"),]

cat("Filter I - kappa chain only",file = logfile, sep="\n")
cat(paste0("Reads removed: ",length(which(infile[,chaincol]!="kappa"))),file = logfile, sep="\n")
cat(paste0("Infile I rows: ",dim(infile_I)[1]),file = logfile, sep="\n")
cat("-----",file = logfile, sep="\n")

# filter 2 - “_nt_identity_v" – 95% and above
nt_identity_v_col<-grep("var_identity_nt",colnames(infile))
infile_II<-infile_I[which(infile_I[,nt_identity_v_col]>=95),]

cat("Filter II - nt_identity_v – 95% and above",file = logfile, sep="\n")
cat(paste0("Reads removed: ",length(which(infile_I[,nt_identity_v_col]<95))),file = logfile, sep="\n")
cat(paste0("Infile II rows: ",dim(infile_II)[1]),file = logfile, sep="\n")
cat("-----",file = logfile, sep="\n")

# filter III 
# Vk gene length >= 150bp and <=311bp

Vk_gene_length_col<-grep("v_gene_length",colnames(infile_II))

glength<-apply(infile_II,1,function(x) nchar(x[Vk_gene_length_col]))

infile_III<-infile_II[which(glength>=vlower&glength<=vupper),]

cat(paste0("Filter III - Vk gene length >= ",vlower," bp and <= ",vupper," bp"),file = logfile, sep="\n")
cat(paste0("Reads removed: ",length(which(glength<150))),file = logfile, sep="\n")
cat(paste0("Infile III rows: ",dim(infile_III)[1]),file = logfile, sep="\n")
cat("-----",file = logfile, sep="\n")


### Deduplication - gDNA
# I - all VJ combos

dedup_gdna<-function(infile_III){
None_raw_input_col<-grep("raw_input",colnames(infile_III))
None_seq_id_col<-grep("seq_id",colnames(infile_III))

dedup_out<-vector()
VJ<-unique(paste(infile_III$v_gene,infile_III$j_gene,sep="_"))

# for each combo, take pools of the same length

#getVJ_len_dups<-function(VJ_combo,infile){

flagger<-function(VJ_adapt,ad_pattern){
  NNs<-sapply(VJ_adapt$raw_input,function(x) substr(x,str_locate(x,ad_pattern)[1]-6,str_locate(x,ad_pattern)[1]-1))
  
  dedup<-VJ_adapt$seq_id[duplicated(NNs)]
  dedup
}

for(j in 1:length(VJ)){

VJ_combo<-VJ[j]
print(VJ_combo)

V<-strsplit(VJ_combo,"_")[[1]][1]
J<-strsplit(VJ_combo,"_")[[1]][2]
VJ_tab<-infile_III[which(infile_III$v_gene==V&infile_III$j_gene==J),c(None_seq_id_col,None_raw_input_col)]
#VJ_tab<-infile[which(infile$v_gene==V&infile$j_gene==J),]
VJ_len<-nchar(VJ_tab$raw_input)

n_occur<-data.frame(table(nchar(VJ_tab$raw_input)))

VJ_dups<-VJ_tab[VJ_len%in%n_occur$Var1[n_occur$Freq>1],]
VJ_len<-nchar(VJ_dups$raw_input)

for(i in 1:length(unique(VJ_len))){
  
  VJ_len_dups<-VJ_dups[VJ_len%in%unique(VJ_len)[i],]
  adapt1<-str_detect(VJ_len_dups$raw_input,"GACTCGT")
  adapt2<-str_detect(VJ_len_dups$raw_input,"CTGCTCCT")
  
  VJ_adapt1<-VJ_len_dups[adapt1,]
   if(dim(VJ_adapt1)[1]>1){
     dedup<-flagger(VJ_adapt1,"GACTCGT")
     if(length(dedup)>0){dedup_out=c(dedup_out,dedup)}
   }
  
  VJ_adapt2<-VJ_len_dups[adapt2,]
  if(dim(VJ_adapt2)[1]>1){
    dedup<-flagger(VJ_adapt2,"CTGCTCCT")
    if(length(dedup)>0){
      dedup_out=c(dedup_out,dedup)
      print(dedup)
      }
  }
  
}

}

infile_IV<-infile_III[!infile_III$seq_id%in%dedup_out,]
}

# infile_IV<-dedup_gdna(infile_III)

### Deduplication - RNA 

# I - all VJ combos

dedup_rna<-function(infile_rna_III){

None_raw_input_col<-grep("raw_input",colnames(infile_rna_III))
None_seq_id_col<-grep("seq_id",colnames(infile_rna_III))
  
dedup_out<-vector()
VJ<-unique(paste(infile_rna_III$v_gene,infile_rna_III$j_gene,sep="_"))

# for each combo, take pools of the same length

#getVJ_len_dups<-function(VJ_combo,infile){

for(j in 1:length(VJ)){
  
  VJ_combo<-VJ[j]
  print(VJ_combo)
  
  
  V<-strsplit(VJ_combo,"_")[[1]][1]
  J<-strsplit(VJ_combo,"_")[[1]][2]
  VJ_tab<-infile_rna_III[which(infile_rna_III$v_gene==V&infile_rna_III$j_gene==J),c(None_seq_id_col,None_raw_input_col)]
#  print(nrow(VJ_tab))
  
  vjset<-DNAStringSet(VJ_tab$raw_input)
  adapt1_Ir<-vmatchPattern("CCACGCGTGCCCTAT",vjset,min.mismatch = 0,max.mismatch = 1) # match first 15 bp of adaptor, one mismatch)
  adapt1<-unlist(lapply(adapt1_Ir,function(x) length(x)))
      
      
#      str_detect(VJ_tab$raw_input,"CCACGCGTGCCCTATAGT")
#    adapt2<-str_detect(VJ_len_dups$raw_input,"CTGCTCCT")
    
    VJ_adapt1<-VJ_tab[adapt1>0,]
    vjset<-DNAStringSet(VJ_adapt1$raw_input)
    adapt1f_Ir<-vmatchPattern("CCACGCGTGCCCTAT",vjset,min.mismatch = 0,max.mismatch = 1)
    adapt1f<-adapt1[adapt1>0]
    
    starts<-vector()
    NNs<-vector()
    
    if(nrow(VJ_adapt1)>0){       # handle nrow(VJ_adapt1)==0 exception
    for(i in 1:dim(VJ_adapt1)[1]){
      
      st<-adapt1f_Ir[[i]][[adapt1f[i]]][1]
      NNs=c(NNs,substr(VJ_adapt1$raw_input[i],st-6,st-1))
      starts=c(starts,st)
     
    } 
    }else NNs<-0
    
      dedup<-VJ_adapt1$seq_id[duplicated(NNs)]
      
      if(length(dedup)>0){dedup_out=c(dedup_out,dedup)}
    
}

infile_rna_IV<-infile_rna_III[!infile_rna_III$seq_id%in%dedup_out,]
}

#infile_rna<-read.csv("EK-RNA-sort041717-062717_S10_mergedpairedend_15overlap_sample.csv")

# filter I - kappa chain only
#chaincol<-grep("chain",colnames(infile_rna))
#infile_rna_I<-infile_rna[which(infile_rna[,chaincol]=="kappa"),]

# filter 2 - “_nt_identity_v" – 95% and above
#nt_identity_v_col<-grep("_nt_identity_v",colnames(infile_rna_I))
#infile_rna_II<-infile_rna_I[which(infile_rna_I[,nt_identity_v_col]>=95),]

# filter III 
# Vk gene length >= 150bp.

#Vk_gene_length_col<-grep("_germ_alignments_nt_var_query",colnames(infile_rna_II))

#glength<-apply(infile_rna_II,1,function(x) nchar(x[Vk_gene_length_col]))

#infile_rna_III<-infile_rna_II[which(glength>=150),]

if(dedup==T){
  
  if(input_type=="gDNA"){
    infile_IV<-dedup_gdna(infile_III)
   
  }else if(input_type=="RNA"){
    infile_IV<-dedup_rna(infile_III)
    
  }else{
    stop("please specify a valid input type ('gDNA' or 'RNA')")
  }
  
  cat(paste0("Deduplication - ",input_type),file = logfile, sep="\n")
  cat(paste0("Reads removed: ",nrow(infile_III)-nrow(infile_IV)),file = logfile, sep="\n")
  cat(paste0("Infile IV rows: ",dim(infile_IV)[1]),file = logfile, sep="\n")
  cat("-----",file = logfile, sep="\n")
  
  
}else infile_IV<-infile_III

## Output stats

# I. All Jks combined

# 1. Percentage of each V gene (Column F) from total number of reads

# v_genes<-unique(infile_IV$v_gene)

vgenes<-table(infile_IV$v_full)

vgenes_percent<-(vgenes/dim(infile_IV)[1])*100

# 2. Percentage of productive and non-productive (Column D) from the total reads

prod<-table(infile_IV$productive)
prod_percent<-(prod/dim(infile_IV)[1])*100

# 3. V gene percentage within productive and non-productive
# (e.g  % Vk1-88 non-productive within non-productive pool or within total pool, % Vk1-88 productive within productive pool or total pool)

v_genes<-unique(infile_IV$v_full)
v_genes_prod<-sapply(v_genes,function(x) length(which(infile_IV$v_full==x&infile_IV$productive=="yes")))
v_genes_nonprod<-sapply(v_genes,function(x) length(which(infile_IV$v_full==x&infile_IV$productive=="no")))

prod_total<-(v_genes_prod/dim(infile_IV)[1])*100
prod_prod<-(v_genes_prod/prod[2])*100

nonprod_total<-(v_genes_nonprod/dim(infile_IV)[1])*100
nonprod_nonprod<-(v_genes_nonprod/prod[1])*100

outtab_allJks<-cbind(vgenes[v_genes],vgenes_percent[v_genes],v_genes_prod,prod_total,prod_prod,v_genes_nonprod,nonprod_total,nonprod_nonprod)
#outtab_allJks=outtab_allJks[order(rownames(outtab_allJks)),]
colnames(outtab_allJks)[1:2]<-c("total occurrences","% total occurrences")

v_genes_ord<-sapply(v_genes,function(x) as.numeric(strsplit(x,"[-*]")[[1]][2]))
outtab_allJks=outtab_allJks[order(v_genes_ord),]

library(xlsx)
write.xlsx(outtab_allJks,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_output_stats_allJk.xlsx"))

# gDNA: keep only exact matches of Jk primer

if(input_type=="gDNA"){

  
###
  
exact_m<-read.delim(exact_matches)  

matchLRpatterns_bulk<-function(exact_m,seq){
 strseq<-BString(seq)
 out<-apply(exact_m,1,function(x)length(matchLRPatterns(x[2],x[3],0,strseq,max.Rmismatch =1))>0)
 out
}
    
infile_exact<-sapply(infile_IV$raw_input,function(x) any(matchLRpatterns_bulk(exact_m,x)))

######

#Jk_seqs<-c("CGTTCGGTGGAGGCACCAAGCTGGAAAT[ACTG]","ACGTTCGGAGGGGGGACCAAGCTGGAAATAAAACGTAAG[ACTG]","TTCACGTTCGGCTCGGGGACAAAGT[ACTG]","CGTTCGGTGCTGGGACCAAGCTGGAGCTGAAAC[ACTG]")


#infile_exact<-sapply(infile_IV$None_raw_input,function(x) any(str_detect(x,Jk_seqs)))

infile_IV_exact<-infile_IV[infile_exact,]

cat(paste0("Jk exact matches: ",dim(infile_IV_exact)[1]),file = logfile, sep="\n")

# reassign Jks
# Jk1 reassignment

cross_p<-read.delim(cross_priming)

Jk1_reassign<-subset(cross_p,Gene=="Jk1")
Jk2_reassign<-subset(cross_p,Gene=="Jk2")
Jk4_reassign<-subset(cross_p,Gene=="Jk4")
Jk5_reassign<-subset(cross_p,Gene=="Jk5")

Jk_reassign<-function(infile_IV,Jk_reassign_seqs,target_reassign){
  
None_raw_input_col<-grep("raw_input",colnames(infile_IV))
None_j_gene_gene_col<-grep("j_gene",colnames(infile_IV))

#Jk_reassign_count<-sapply(infile_IV$None_raw_input,function(x) any(str_detect(x,Jk_reassign_seqs)))
Jk_reassign_count<-sapply(infile_IV$raw_input,function(x) any(matchLRpatterns_bulk(Jk_reassign_seqs,x)))

block<-infile_IV[Jk_reassign_count,]
if(nrow(block)>0) block$j_gene<-target_reassign
#Jk_reassign_vec<-apply(infile_IV,1,function(x) ifelse(any(str_detect(x[None_raw_input_col],Jk_reassign_seqs)),target_reassign,x[None_j_gene_gene_col]))

cat(paste0(target_reassign," reassignments: ",length(which(Jk_reassign_count))),file = logfile, sep="\n")
    
return(list(block,Jk_reassign_count))
}

infile_IV_crossampl<-vector()

Jk1_reassign_vec<-Jk_reassign(infile_IV,Jk1_reassign,"IGKJ1")
infile_IV_crossampl=rbind(infile_IV_crossampl,Jk1_reassign_vec[[1]])

Jk2_reassign_vec<-Jk_reassign(infile_IV,Jk2_reassign,"IGKJ2")
infile_IV_crossampl=rbind(infile_IV_crossampl,Jk2_reassign_vec[[1]])

Jk4_reassign_vec<-Jk_reassign(infile_IV,Jk4_reassign,"IGKJ4")
infile_IV_crossampl=rbind(infile_IV_crossampl,Jk4_reassign_vec[[1]])

Jk5_reassign_vec<-Jk_reassign(infile_IV,Jk5_reassign,"IGKJ5")
infile_IV_crossampl=rbind(infile_IV_crossampl,Jk5_reassign_vec[[1]])

Jks_vec<-paste0("Jk",c(1,2,4,5),"_reassign_vec")
cat(paste0("Total: ", sum(sapply(Jks_vec,function(x) length(which(get(x)[[2]]))))),file = logfile, sep="\n")

infile_IV_combined<-rbind(infile_IV_exact,infile_IV_crossampl)
infile_IV_combined=infile_IV_combined[!duplicated(infile_IV_combined),]

}

individual_Jk_stats<-function(infile,targetJk){
  v_genes<-unique(infile$v_full)
  block<-vector()
  for(x in 1:length(v_genes)){
    total_Jk_gene<-infile[which(infile$v_full==v_genes[x]&infile$j_gene==targetJk),]
    Jk_total<-length(which(infile$j_gene==targetJk))
    total_Jk_gene_percent<-(nrow(total_Jk_gene)/Jk_total)*100
    Jk_prod<-length(which(total_Jk_gene$productive=="yes"))
    
    Jk_total_prod<-length(which(infile$j_gene==targetJk&infile$productive=="yes"))
    Jk_total_nonprod<-length(which(infile$j_gene==targetJk&infile$productive=="no"))
    
    Jk_prod_prod_percent<-(Jk_prod/Jk_total_prod)*100
    Jk_prod_total_percent<-(Jk_prod/Jk_total)*100
    
    Jk_nonprod<-length(which(total_Jk_gene$productive=="no"))
    Jk_nonprod_nonprod_percent<-(Jk_nonprod/Jk_total_nonprod)*100
    Jk_nonprod_total_percent<-(Jk_nonprod/Jk_total)*100
    
    block=rbind(block,c(nrow(total_Jk_gene),total_Jk_gene_percent,Jk_prod,Jk_prod_prod_percent,Jk_prod_total_percent,Jk_nonprod,Jk_nonprod_nonprod_percent,Jk_nonprod_total_percent))
  }
  colnames(block)<-c(targetJk,paste0(targetJk," %"),paste0(targetJk,"_prod"),paste0(targetJk,"_prod/total_",targetJk,"_prod"),paste0(targetJk,"_prod/total_",targetJk),paste0(targetJk,"_nonprod"),paste0(targetJk,"_nonprod/total_",targetJk,"_nonprod"),paste0(targetJk,"_nonprod/total_",targetJk))
  rownames(block)<-v_genes
  v_genes_ord<-sapply(v_genes,function(x) as.numeric(strsplit(x,"[-*]")[[1]][2]))
  block=block[order(v_genes_ord),]
  block
}

infile_IV_simple<-infile_IV[,c("v_full","v_gene","j_full","j_gene","seq_id","chain","cdr3_nt","cdr3_length","var_identity_nt","v_gene_length","productive","raw_input")]


if(input_type=="RNA"){
  Jk_individual_stats_tab<-do.call("cbind",lapply(c("IGKJ1","IGKJ2","IGKJ4","IGKJ5"),function(x) individual_Jk_stats(infile_IV,x)))
  write.xlsx(Jk_individual_stats_tab,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_output_stats_individual_Jks.xlsx"))
  write.xlsx2(infile_IV_simple,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_Abstar_input_simplified.xlsx"))

}else if(input_type=="gDNA"){
  Jk_stats_exact_tab<-do.call("cbind",lapply(c("IGKJ1","IGKJ2","IGKJ4","IGKJ5"),function(x) individual_Jk_stats(infile_IV_exact,x)))
  write.xlsx(Jk_stats_exact_tab,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_output_stats_individual_Jks_exact.xlsx"))

  Jk_stats_combined_tab<-do.call("cbind",lapply(c("IGKJ1","IGKJ2","IGKJ4","IGKJ5"),function(x) individual_Jk_stats(infile_IV_combined,x)))
  write.xlsx(Jk_stats_combined_tab,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_output_stats_individual_Jks_combined.xlsx"))
  
  exact_match<-sapply(infile_IV$seq_id,function(x) ifelse(x%in%infile_IV_exact$seq_id,"yes","no"))
  crossampl_match<-sapply(infile_IV$seq_id,function(x) ifelse(x%in%infile_IV_crossampl$seq_id,"yes","no"))
  exact_or_crossampl_match<-sapply(infile_IV$seq_id,function(x) ifelse(x%in%infile_IV_exact$seq_id|x%in%infile_IV_crossampl$seq_id,"yes","no"))
  
  infile_IV_simple=cbind(infile_IV_simple,exact_match,crossampl_match,exact_or_crossampl_match)
  write.xlsx2(infile_IV_simple,file=paste0(sub(".csv","",input_file),"_dedup_",dedup,"_Abstar_input_simplified.xlsx"))
  
}

close(logfile)

}

HGTS_processing(opt$input_file,opt$exact_matches,opt$cross_priming,opt$input_type,opt$dedup,opt$lower_length,opt$upper_length)


