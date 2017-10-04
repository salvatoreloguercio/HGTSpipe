# HTGS abstar processing

library(stringr)

infile<-read.csv("EK-gDNA-053017_S22_mergedpairedend_15overlap_sample.txt")

# filter I - kappa chain only

infile_I<-infile[which(infile$chain=="kappa"),]

# filter 2 - “var_identity_nt” column S – 95% and above

infile_II<-infile_I[which(infile_I$var_identity_nt>=95),]

# filter III - gDNA
# Vk gene length >= 150bp.
# For gDNA samples, Vk gene length is the length of the vdj_nt read length (column O) minus 34 for Jk1 reads, 39 for Jk2 reads, 27 for Jk4 reads and 38 for Jk5 reads (based on where nested primer falls on Jk exon) (Column J).  

jkadder_dna<-function(jgene){
  if(jgene=="IGKJ1"){jkadd<-34}
  else if(jgene=="IGKJ2"){jkadd<-39}
  else if(jgene=="IGKJ4"){jkadd<-27}
  else if(jgene=="IGKJ5"){jkadd<-38}
  jkadd
}

glength<-apply(infile_II,1,function(x) nchar(x[15])-jkadder_dna(x[10]))

infile_III<-infile_II[which(glength>=150),]

# filter III - RNA

infile_rna<-read.csv("EK-RNA-sort041717-062717_S10_mergedpairedend_15overlap_sample.txt")
infile_rna_I<-infile_rna[which(infile_rna$chain=="kappa"),]
infile_rna_II<-infile_rna_I[which(infile_rna_I$var_identity_nt>=95),]


jkadder_rna<-function(jgene){
  if(jgene=="IGKJ1"){jkadd<-38}
  else if(jgene=="IGKJ2"){jkadd<-39}
  else if(jgene=="IGKJ4"){jkadd<-38}
  else if(jgene=="IGKJ5"){jkadd<-38}
  jkadd
}

glength_rna<-apply(infile_rna_II,1,function(x) nchar(x[15])-jkadder_rna(x[10]))

infile_rna_III<-infile_rna_II[which(glength_rna>=150),]

### Deduplication - gDNA
# I - all VJ combos

dedup_out<-vector()
VJ<-unique(paste(infile_III$v_gene,infile_III$j_gene,sep="_"))

# for each combo, take pools of the same length

#getVJ_len_dups<-function(VJ_combo,infile){

flagger<-function(VJ_adapt,ad_pattern){
  NNs<-sapply(VJ_adapt$raw_input,function(x) substr(x,str_locate(x,ad_pattern)[1]-6,str_locate(x,ad_pattern)[1]-1))
  
  dedup<-VJ_adapt[duplicated(NNs),1]
  dedup
}

for(j in 1:length(VJ)){

VJ_combo<-VJ[j]
print(VJ_combo)

V<-strsplit(VJ_combo,"_")[[1]][1]
J<-strsplit(VJ_combo,"_")[[1]][2]
VJ_tab<-infile_III[which(infile_III$v_gene==V&infile_III$j_gene==J),]
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


### Deduplication - RNA 

# I - all VJ combos

dedup_out<-vector()
VJ<-unique(paste(infile_rna_III$v_gene,infile_rna_III$j_gene,sep="_"))

# for each combo, take pools of the same length

#getVJ_len_dups<-function(VJ_combo,infile){

for(j in 1:length(VJ)){
  
  VJ_combo<-VJ[j]
  print(VJ_combo)
  
  V<-strsplit(VJ_combo,"_")[[1]][1]
  J<-strsplit(VJ_combo,"_")[[1]][2]
  VJ_tab<-infile_rna_III[which(infile_rna_III$v_gene==V&infile_rna_III$j_gene==J),]
  #VJ_tab<-infile[which(infile$v_gene==V&infile$j_gene==J),]
 # VJ_len<-nchar(VJ_tab$raw_input)
  
#  n_occur<-data.frame(table(nchar(VJ_tab$raw_input)))
  
#  VJ_dups<-VJ_tab[VJ_len%in%n_occur$Var1[n_occur$Freq>1],]
#  VJ_len<-nchar(VJ_dups$raw_input)
  
 # for(i in 1:length(unique(VJ_len))){
    
#    VJ_len_dups<-VJ_dups[VJ_len%in%unique(VJ_len)[i],]
  library(Biostrings)
  
  vjset<-DNAStringSet(VJ_tab$raw_input)
  adapt1_Ir<-vmatchPattern("CCACGCGTGCCCTAT",vjset,min.mismatch = 0,max.mismatch = 1) # match first 15 bp of adaptor, one mismatch)
  adapt1<-unlist(lapply(a,function(x) length(x)))
      
      
#      str_detect(VJ_tab$raw_input,"CCACGCGTGCCCTATAGT")
#    adapt2<-str_detect(VJ_len_dups$raw_input,"CTGCTCCT")
    
    VJ_adapt1<-VJ_tab[adapt1>0,]
    vjset<-DNAStringSet(VJ_adapt1$raw_input)
    adapt1f_Ir<-vmatchPattern("CCACGCGTGCCCTAT",vjset,min.mismatch = 0,max.mismatch = 1)
    adapt1f<-adapt1[adapt1>0]
    
    starts<-vector()
    NNs<-vector()
    for(i in 1:dim(VJ_adapt1)[1]){
      
      st<-adapt1f_Ir[[i]][[adapt1f[i]]][1]
      NNs=c(NNs,substr(VJ_adapt1$raw_input[i],st-6,st-1))
      starts=c(starts,st)
     
    } 
      
      dedup<-VJ_adapt1[duplicated(NNs),1]
      
      if(length(dedup)>0){dedup_out=c(dedup_out,dedup)}
    
}

infile_rna_IV<-infile_rna_III[!infile_rna_III$seq_id%in%dedup_out,]


## Output stats

# I. All Jks combined

# 1. Percentage of each V gene (Column F) from total number of reads

# v_genes<-unique(infile_IV$v_gene)

vgenes<-table(infile_IV$v_gene)

vgenes_percent<-(vgenes/dim(infile_IV)[1])*100

# 2. Percentage of productive and non-productive (Column D) from the total reads

prod<-table(infile_IV$productive)
prod_percent<-(prod/dim(infile_IV)[1])*100

# 3. V gene percentage within productive and non-productive
# (e.g  % Vk1-88 non-productive within non-productive pool or within total pool, % Vk1-88 productive within productive pool or total pool)

v_genes<-unique(infile_IV$v_gene)
v_genes_prod<-sapply(v_genes,function(x) length(which(infile_IV$v_gene==x&infile_IV$productive=="yes")))
v_genes_nonprod<-sapply(v_genes,function(x) length(which(infile_IV$v_gene==x&infile_IV$productive=="no")))

prod_total<-(v_genes_prod/dim(infile_IV)[1])*100
prod_prod<-(v_genes_prod/prod[2])*100

nonprod_total<-(v_genes_nonprod/dim(infile_IV)[1])*100
nonprod_nonprod<-(v_genes_nonprod/prod[1])*100

vgenes=vgenes[v_genes]

# II. Jks 1,2,4,5 separated









######################## Scraps
# match adaptor

sapply(infile_rna_III$raw_input,admatcher)

admatcher<-function(x){
vjset<-DNAStringSet(infile_rna_III$raw_input)

a<-vmatchPattern("CCACGCGTGCCCTATAGTC",vjset)
b<-unlist(lapply(a,function(x) length(x)))









