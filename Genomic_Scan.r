#!/usr/bin/env Rscript

#========================================================================================#
#========================================================================================#
# SCANS GENOMIC : 

  # - ACP Global on all chromosome 
  # - ACP on sliding windows
  # - ACP on SNP range
  # - FST on sliding windows
  # - Global FST on sliding windows
  # - FST by SNP by chromosome
  # - PCadapt 
  # - Fixed genotype by group in SW
  # - Diversity (theta and TD) 
  # - GLOBAL SFS by Chromosome
  # - SFS in Sliding Windows: Euclidian distance of SFS in sliding windows vs global SFS on chromosome
  # - LD : to do
  # 

# Aim : Run all scan independently on each chr 

# Authors : Stefano Mona , Elise Gay, Romuald Laso-Jadart (MNHN - EPHE)
#           2023
#	Please inform the authors before sharing
#========================================================================================#
#========================================================================================#

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#

library(bedr)
library(qvalue)
library(vcfR)
library(rlist)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(reshape2)
library(pcadapt)
library(ggrepel)
library(gridExtra)
library(ggbreak)
library(IRanges)
library(stringr)
library(PopGenome)
library(basetheme)
library(network)

# Stefano Mona libraries housemade
#----------------------------------#
source("R_functions/libreria_filtri_VCF_e_SFS_unfolded.r")

# function to apply filters on genotype and depth
#--------------------------------------------------#
source("R_functions/fonction_filtre_cover.R")
source("R_functions/Filter_Ref.R")
source("R_functions/Filter_Het.R")
source("R_functions/Filter_Hom_ref.R")
source("R_functions/Filter_Na.R")

# get/set dir
#---------------#
# get directories where VCF needed for genomic scan are stored
dir_data="data/"
dir(dir_data)

# create folder for each scan (comment line if it is a re-run)
# check your current directory
getwd()

dir.create("Sliding_div")
dir.create("SlidingSFS")
dir.create("Sliding_PCA")
dir.create("SlidingFst")
dir.create("FstBySNP")
dir.create("Sliding_global_FST")
dir.create("pcadapt")
dir.create("freq_genotype/")

# Genomic scan run on each chromosome (1 vcf by chromosomes are needed)
# For each index : 
# - Needed input VCF, position vector and table of metadata about samples and populations are read in each part
# - Result table by windows are stored in corresponding folder
# - Plot are made for each chromosome and stored in PDF in corresponding folder

#========================================#
#========================================#
# ------ Load Metadata and variable ----
#========================================#
#========================================#

# Chr list : to run genomic scan
#---------#
chr_list=c("SUPER_10", "SUPER_14")


# Metadata pop
#--------------#

# the pop table has to be ordered in the same way as all the VCF header
#'''
#samples	pop	
#sample_1	pop1
#sample_2	pop1
#sample_3	pop1
#sample_4	pop2
#sample_5	pop2
#'''

table_pop=read.table("metadata/Samples_table.txt", header = TRUE, row.names = 1)
table_pop

samples=row.names(table_pop)
samples

pop=table_pop$pop
pop

sex=table_pop$sex
sex

pop1=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP1"),]))
pop1
pop2=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP2"),]))
pop2
pop3=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP3"),]))
pop3
pop4=row.names(as.data.frame(table_pop[which(table_pop$pop == "POP4"),]))
pop4

lista_pop<-list(pop1,pop2,pop3,pop4)
names(lista_pop)=c("pop1","pop2","pop3","pop4")

# for div computing 
lista_pop_all=list(pop1,pop2,pop3,pop4,c(pop1,pop2,pop3,pop4))
names(lista_pop_all)=c("pop1","pop2","pop3","pop4","pop1+pop2+pop3+pop4")

# get all genes coordinates in a bed file :
#--------------------------------------------#
gene=read.table("metadata/Genes.name.bed", header = TRUE)

# prefix of vcf to read depending on filters (according to the exact name of your VCF in yout local folder)
# Only put the 'base name' , the loop in each scan will fill the chr name 
#-----------------------------------------------------------------------------------------------------------#
VCF_0="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na0_60.vcf.gz"
VCF_20="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60.vcf.gz"
VCF_20_MAF005="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60_MAF005.vcf.gz"
VCF_20_MAF005_unzip="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60_MAF005.vcf"
position_20_MAF ="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60_MAF005.positions"
position_all="_TAG_Flowqual_Noindels_Norepeat_SNP_21ind.positions"
position_20="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60.positions"
position_0="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na0_60.positions"
BED_20_MAF005="_TAG_Flowqual_Noindels_Norepeat_SNP_21IND_DP10_50_10_200_Na20_60_MAF005.bed.bed"

# prepare basename for variable save in R env
base_name_20="VCF_F_20Na"
base_name_0="VCF_F_0Na"
base_name_position="VCF_F.position"
base_name_position_0="VCF_0.position"
base_name_position_20_MAF="VCF_20_MAF.position"
base_name_plink="Bed_Na20_MAF005"
base_name_MAF="VCF_F_20Na_MAF005"
base_name_position_20="VCF_20.position"

# read chr table (needed several time to plot index along chr)
# '''
# Chr		length
# SUPER_1	274107945
# SUPER_2	239769335
# SUPER_3	233268383
# SUPER_4	205643877
# ...
# '''
table_chr=read.table("metadata/table_chr_length.txt", header = TRUE)

#======================# 
#======================# 
# ------ ACP All chr----
#======================#
#======================#

# Compute PCA on the entire VCF given

# INPUT 
#-------#
# VCF with chosen filters : 20% Na and MAF filtered (MAF can also be put directly in the PCadapt function)
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# FUNCTION ARGS
#------------------#
# respca_10 = pcadapt(input = read.pcadapt(path_vcf, type = "vcf"),  min.maf=0.01, K = 10)
  # VCF : abslotue path to your VCF file (already filtered)
  # min.maf = minimum maf threshold
  # K = number of components to be computed

# OUTPUT
#---------#
# PCA plot according to  and population assignation (see metadata)

#----------#
# run ACP
#----------#
# To adapt before run the Loop :
# name_vcf = to save as variable, here "base_name_MAF"
# path_vcf = to read your vcf
# in ggsave : filename = and Path=
# scale_colour_manual : in ggplot, adapt the number of color and label to your pop

for (chr in chr_list){
  name_vcf=paste(base_name_MAF, chr, sep="_")
  path_vcf=paste(c(dir_data,chr,VCF_20_MAF005), collapse = "")
  respca_10 = pcadapt(input = read.pcadapt(path_vcf,
                                           type = "vcf"),
                      min.maf=0.05, K = 10)

  # add samples name in table result (same order as VCF)
  scores = data.frame(respca_10$scores)
  rownames(scores) = samples

  # Compute percent of variance explained by all axis (here 10 axis)
  PEV = paste(round((respca_10$singular.values^2)*100,digits=1),"%",sep="")

  # Plot PCA
  p=ggplot() +
    geom_point(data=scores,
               aes(x=scores$X1, y=scores$X2,
                   colour=pop),
               size=2)+

    labs(x=PEV[1], y = PEV[2]) +

    ggtitle(chr) +

    scale_colour_manual(name = "Sampling sites",
                         values=rainbow(length(unique(pop)))) +

    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5))
  
  ggsave(p,path="Sliding_PCA/", filename = paste(chr,"_PCA.pdf"), width = 10, height = 10, device = "pdf")
}

#===================================#
#===================================#
# ------ ACP Sliding windows-------
#===================================#
#===================================#

# info
#-------#
# INPUT
#------#
# VCF input : Noindels_NoRef_Na20.vcf.gz
# lista_pop = list of populations containing ID in each (list of list type) : created in metadata part

# FUNCTION : ARGS
#----------------#
# sliding_pca_adegenet_basepair<-function(dati, window, slide, min_n_snp)
  # dati : your good vcf (with filters for structure analyses) 
  # window : the size of the window (int)
  # slide : the size of the slide (int)
  # min_n_snp: minimum Nb of SNPs to compute PCA (100 should be enough, it is just to avoid error at the bounds of the scaffold)

# OUTPUT
#---------#
# Return the final table with windows and axis info for the whole VCF
# finale<-cbind(low_bound,upper_bound,pca_asse1, pca_asse2)	

# in "risultati" you have a matrix of four columns: 
# col1 = lower value of the window
# col2 = upper value of the window
# col3 = % of variance of the first axe
# col4 = % of variance of the second axe
# you may plot for example the average values between col1 & col2 (i.e.: more or less the midpoint of the window) against either col3 or col3+col4

#----------------------------#
# Run Stefano's functions
# & create plot
#----------------------------#
# To adapt to your data :
# name_vcf = to save as variable
# path_vcf = to read your vcf
# in ggsave : filename = and Path=

window=1000000
slide=500000
min_SNP=300

# Compute PC1 variance in SW
#-----------------------------#

for (chr in chr_list){
  # get VCF
  name_vcf=paste(base_name_MAF, chr, sep="_")
  path_vcf=paste(c(dir_data,chr,VCF_20_MAF005), collapse = "")
  name_vcf<-read.vcfR(path_vcf, verbose = TRUE) 

  # run pca
  risultati<-sliding_pca_adegenet_basepair(name_vcf,window, slide, min_n_snp=min_SNP)
  risultati=as.data.frame(risultati)
  
  # get number of SNP in interval saved in  'Nb_SNP_vec' vector
  vettore_snps_pos<-getPOS(name_vcf)
  low_bound<-seq(1, max(vettore_snps_pos)-slide, slide)
  upper_bound<-low_bound+window
  Nb_SNP_vec=c()
  for (i in 1:length(low_bound)) {
    Nb_SNP_vec_current=length(which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] ))
    if (Nb_SNP_vec_current > min_SNP){Nb_SNP_vec=c(Nb_SNP_vec, Nb_SNP_vec_current)}
    else{Nb_SNP_vec=c(Nb_SNP_vec, "Na")}
  }
  
  # Create and format table
  risultati$Nb_SNP_PCA=Nb_SNP_vec
  risultati$mid_wind=rowMeans(risultati[,c("low_bound", "upper_bound")])
  risultati=risultati[,c("low_bound", "upper_bound", "mid_wind", "Nb_SNP_PCA","pca_asse1", "pca_asse2")]
  colnames(risultati)=c("low_bound_PC", "upper_bound_PC", "mid_wind_PC", "Nb_SNP_PCA","pca_asse1", "pca_asse2")
  # save table
  write.table(risultati, paste("Sliding_PCA/", chr, "_",min_SNP,"_PCA_Na20_MAF.table", sep=""))

}

# Plot PC1 variance along chr
#-----------------------------#

for (chr in chr_list){
  table=read.table(paste("Sliding_PCA/", chr, "_",min_SNP,"_PCA_Na20_MAF.table", sep=""), header=TRUE)
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]
  
  # get current chr length
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  # plot 
  p=ggplot() +
    geom_point(aes(x=table$mid_wind_PC, y=table$pca_asse1, 
                   colour="black"),
               size=1)+
    
    labs(x="500K pb windows", y = "PC 1") +
    
    ggtitle(chr) +
    theme(text = element_text(size = 10), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
										  
    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000))) +
									
	# (optional, uncomment if needed) print start and end of gene on chr with line
    #geom_vline(xintercept = coord_gene$start, 
    #           color = "red") +
    #geom_vline(xintercept = coord_gene$end, 
    #           color = "blue") +
			   
    scale_y_continuous(breaks = seq(0, 1, by = 0.2))
  
  ggsave(p, path="Sliding_PCA/", width = 12, height = 12, device = "pdf", filename = paste(chr, "_",min_SNP,"_PCA_Na20_MAF.pdf", sep = ""))
}

#===================================#
#===================================#
# ------ ACP on SNPs range -------
#===================================#
#===================================#
# info
#-------#
# INPUT
# VCF input : Noindels_NoRef_Na20.vcf.gz

# Pipeline :
# Segment vcf file every 'n' SNPs rather than use sliding windows

# OUTPUT
# Return the final table with windows contained more than 'nb_SNP' chosen and their corresponding PCA values

# In 'table_PCA' you have a matrix of four col: 
# col1 = start position of the window
# col2 = end position of the window
# col3 = Mid size of the windows
# col3 = % of variance of the first axe
# col4 = % of variance of the second axe
# you may plot for example the average values between col3 (i.e.: more or less the midpoint of the window) against either col3 or col3+col4

#----------------------------#
# Run PCA on range SNP
#----------------------------#

# SET up variable
#-----------------#
# Set the value of n
nb_snp <- 1000

# initiate variables
pca_asse1<-c()
pca_asse2<-c()

table_PCA=c()

# Compute PC1 variance
#-----------------------------#
for (chr in chr_list){
  path_vcf=paste(c(dir_data,chr,VCF_20_MAF005), collapse = "")
  readvcf<-read.vcfR(path_vcf, verbose = TRUE) 
  
  # SPLIT VCF every N SNPs
  #-----------------------#
  # position vector
  vector_SNP<-getPOS(readvcf)
  
  # Cut the vector and store the segmented vectors and indices
  for (i in seq(1,length(vector_SNP)-nb_snp, by=nb_snp)) {
    start <- i
    end <- i + nb_snp

    # Extract the lines corresponding to the segment
    vcf_segment <- readvcf[c(start:end),]
    
    print("SNP start and end : ")
    print(start)
    print(end)
    print("vcf :")
    print(vcf_segment)
    
    # filter VCF with Na (avoid samples with only Na in the current region)
    vcf_segment_F=Filters_na_Ind_Pos(vcf_segment, 95.0, 95.0)
    
    # Get genotype from subseted vcfs
    genotype<-extract.gt(vcf_segment_F, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
    
    # Do not run pca if one sample is 100% Na, the programm will be kill
    if (dim(genotype)[2] < 21 & dim(genotype)[2] < 490){
      print(paste("no pca for the region ", start, end, "in chr ", chr, sep = " "))
      table_PCA<-rbind(table_PCA, c(vector_SNP[start],vector_SNP[end], vector_SNP[end]-vector_SNP[start], (vector_SNP[start]+vector_SNP[end])/2, "Na", "Na", length(getPOS(vcf_segment_F)), chr))
      }
    else{
      # transform matrix
      geno_temp<-t(genotype)
      for (j in 1:ncol(geno_temp)) {
        geno_temp[geno_temp[,j]=="./.",j]<-NA 
        }
      
      # run PCA
      mat_data<-df2genind(geno_temp, ploidy = 2, sep="/", pop=row.names(geno_temp))
      # get PC values
      x<-indpca(mat_data)
      new_eigen<-x$ipca$eig/sum(x$ipca$eig)
      pca_asse1<-new_eigen[1]
      pca_asse2<-new_eigen[2]
	  
      # make final table
      table_PCA_current<-cbind(vector_SNP[start],vector_SNP[end], vector_SNP[end]-vector_SNP[start], (vector_SNP[start]+vector_SNP[end])/2, pca_asse1, pca_asse2, length(getPOS(vcf_segment_F)), chr)
      table_PCA<-rbind(table_PCA, table_PCA_current)
    }
  }
  colnames(table_PCA)=c("low_bound_RSNP", "upper_bound_RSNP", "wind_size_RSNP", "mid_wind_RSNP", "pca_asse1_SNPrange", "pca_asse2_SNPrange", "Nb_SNP_range","chr")
  # save table
  write.table(table_PCA, paste("Sliding_PCA/", chr, "_", nb_snp, "_range_PCA_Na20_MAF.table", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
  # re-initiate variables for the next chromosome
  table_PCA=c()
  table_PCA_current=c()
  pca_asse1<-c()
  pca_asse2<-c()
  # remove vcfs of current chromosome
  list_to_remove=grep("vcf_segment_", ls(), value = TRUE)
  rm(list=list_to_remove)
}

# Plot PC1 variance along chr
#-----------------------------#

for (chr in chr_list){
  pca_result=read.table(paste("Sliding_PCA/", chr, "_", nb_snp, "_range_PCA_Na20_MAF.table", sep=""), header = TRUE, na.strings = "Na")
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  coord_gene=gene[which(gene$Chrom == chr),]
  
  # get current chr length
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  # plot 
  p=ggplot() +
    geom_point(aes(x=pca_result$mid_wind_RSNP, y=pca_result$pca_asse1_SNPrange, 
                   colour="black"),
               size=1)+
    
    labs(x="Pb", y = "PC 1") +
    
    ggtitle(paste(chr," PC1 
                  PCA on range SNPs windows of non limited sizes
                  start and end of gene coordinates (red, blue)", sep="")) +
    theme(text = element_text(size = 10), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
    
    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000))) +
    
	# (optional, uncomment if needed) print start and end of gene on chr with line
  geom_vline(xintercept = coord_gene$start, 
             color = "red") +
  geom_vline(xintercept = coord_gene$end, 
             color = "blue") +
			   
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),limits = c(0,1))
  p
  ggsave(p, path="Sliding_PCA/", width = 10, height = 10, device = "pdf", filename = paste(chr, "_", nb_snp, "_range_PCA_Na20_MAF.pdf", sep = ""))
}

#===============================================#
#===============================================#
# ------ FST on sliding windows on whole chr ----
#===============================================#
#===============================================#

#----------#
# Infos  : 
#----------#

# INPUT 
#-------#
# VCF with chosen filters : 20% Na
# Position vector : vector of total position called in the chr

# FUNCTION ARGS
#---------------------#
# risultati<-sliding_fst_spectre_stand(vcf_name, lista_pop, SW_Size, Overlap_size, 0.01, positions_vector)
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part
# SW_Size=100000 : Window size
# Overlap_size=50000 : slide size
# min_SNP=5 : nb min of SNP (default = 0)


# OUTPUT
#---------------------#
# table of pairwise FST values in each windows 
# "mid_wind","FST","low_bound", "upper_bound","Nb_SNP_FST"

#------------#
# Read VCF : 
#------------#
# Filters : Na max = 20% and basics filters (no indel, no repeat, no Lowqual)

# Read VCF INPUT
for (chr in chr_list){
  name_vcf=paste(c(chr,base_name_MAF), collapse = "")
  path_vcf=paste(dir_data,chr,VCF_20_MAF005, sep="")
  readvcf=read.vcfR(path_vcf)
  assign(name_vcf, readvcf, envir=parent.frame())
}

#-----------------------#
# Read position vector
#-----------------------#
for (chr in chr_list){
  name_pos=paste(c(chr,position_all), collapse = "")
  path_pos=paste(dir_data,name_pos,sep = "")
  readpos=read.table(path_pos, header = FALSE)
  pos_vector=readpos$V2
  var_pos=paste(chr,base_name_position, sep="_")
  assign(var_pos, pos_vector, envir=parent.frame())
}

#-----------------#
# Run the loop
#-----------------#
SW_Size = 1000000
Overlap_size = 500000
min_SNP=0

for (chr in chr_list){
  vcf_name=get(paste(chr,base_name_MAF, sep=""))
  positions_vector=get(paste(chr,base_name_position, sep="_"))
  
  # Run FST
  risultati<-sliding_fst_spectre_stand(vcf_name, lista_pop, SW_Size, Overlap_size, 0.05, positions_vector)
  # Get only FST result in first element of "resulati"
  name_file=paste("FST_SW", chr, sep="_")
  fst_table=as.data.frame(risultati[[1]])
  
  # get number of SNP tot in interval
  vettore_snps_pos<-getPOS(vcf_name)
  low_bound<-fst_table$per_plot-(Overlap_size+0.5)
  upper_bound<-fst_table$per_plot+(Overlap_size-0.5)
  Nb_SNP_vec=c()
  for (i in 1:length(low_bound)) {
    Nb_SNP_vec_current=length(which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] ))
    if (Nb_SNP_vec_current >= min_SNP){ # set nb SNP according to windows boundaries extracted from geno_table which contain only SW with SNP > min_SNP from FST function
      Nb_SNP_vec=c(Nb_SNP_vec, Nb_SNP_vec_current)
    }
    else{Nb_SNP_vec=c(Nb_SNP_vec, "NA")} # Print Na in if Nb SNP in SW < min_SNP in order to keep the rigth number of vector size
  }
  
  # Create and format table
  fst_table$low_bound=fst_table$per_plot-(Overlap_size+0.5)
  fst_table$upper_bound=fst_table$per_plot+(Overlap_size-0.5)
  fst_table$Nb_SNP_FST=Nb_SNP_vec
  colnames(fst_table)=c("mid_wind_FST",
                        "FST_POP1_POP2",
                        "FST_POP1_POP3",
                        "FST_POP1_POP4",
                        "FST_POP2_POP3",
                        "FST_POP2_POP4",
                        "FST_POP3_POP4",
                        "low_bound_FST", "upper_bound_FST","Nb_SNP_FST")
  
  fst_table=fst_table[,c("low_bound_FST", "upper_bound_FST", "mid_wind_FST", "Nb_SNP_FST",
                         "FST_POP1_POP2",
                         "FST_POP1_POP3",
                         "FST_POP1_POP4",
                         "FST_POP2_POP3",
                         "FST_POP2_POP4",
                         "FST_POP3_POP4")]
  
  assign(name_file,fst_table, envir=parent.frame())
  
  # Save FST table
  path_file=paste("SlidingFst/",name_file, sep="")
  path_file
  write.table(fst_table,path_file,quote=F)
  
}

#----------#
# Plot FST
#----------#

# Which FST you want to plot
# all pairs are = c("FST_POP1_POP2", "FST_POP1_POP3", "FST_POP1_POP4", "FST_POP2_POP3", "FST_POP2_POP4","FST_POP3_POP4")

pairs_FST="FST_POP1_POP2"

# Run the FST index
for (chr in chr_list){
  
  # read FST file from local folder
  name_file=paste("FST_SW", chr, sep="_")
  path_file=paste("SlidingFst/",name_file, sep="")
  sliding_fst=read.table(path_file, header = TRUE)
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  coord_gene=gene[which(gene$Chrom == chr),]
  
  # read chr length
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  # Format table : discard negative value
  sliding_fst_melted=melt(sliding_fst[,-c(1,2,4)],id.vars = "mid_wind_FST")
  for (fst in 1:nrow(sliding_fst_melted)) {
    if (sliding_fst_melted[fst,"value"] < 0) {
      sliding_fst_melted[fst,"value"] = 0
    }
  }
  sliding_fst_melted_pairs=subset(sliding_fst_melted, sliding_fst_melted$variable == pairs_FST)
  
  # plot
  plot_name=paste(chr, "_FST_",pairs_FST, sep="")
  p=ggplot() + 
    geom_point(data=sliding_fst_melted_pairs,aes(x=sliding_fst_melted_pairs$mid_wind,y=sliding_fst_melted_pairs$value),pch=1)+
    
    labs(y="FST", x = c(chr,"Windows of",SW_Size))+
    ggtitle("FST")+
    
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1, 
                                     size = 12),
          
          text = element_text(size=12),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
    
    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000))) +
    scale_y_continuous(breaks = seq(0, 0.8, 
                                    by = 0.1), limits = c(0,0.8)) +
    
    # (optional, uncomment if needed) print start and end of gene on chr with line
    geom_vline(xintercept = coord_gene$start, 
               color = "red") +
    geom_vline(xintercept = coord_gene$end, 
               color = "blue")
  
  assign(plot_name, p, envir=parent.frame())
  
  ggsave(path = "SlidingFst/", filename = paste(chr,pairs_FST ,"fst.pdf", sep = "_"), 
         p, device = "pdf", height = 10, width = 12)
}

# ------ FST global ------ 
#=========================================#
#=========================================#
# FST global in sliding windows on all pop
#=========================================#
#=========================================#

#----------#
# Infos  : 
#----------#

# INPUT 
#-------#
# VCF with chosen filters : here 20% Na + MAF 0.05
# Because VCF will be read by popgenome : VCF HAS TO BE GUNZIPED AND PUT IN ONE FILE by FOLDER

# FUNCTION ARGS
#---------------------#
# The function is not optimized to run on sliding windows : the sliding is made outside the function in the loop

# boot_popgenome<-function(dati, lista_pop, bootstrap)
# dati = vcf
# bootstrap (int)
# lista_pop = list of population 

# window=100000 : Window size
# slide=50000 : slide size
# min_SNP=5 : nb min of SNP (default = 0)

# How does it work :
#--------------------#
# run file :
# testdata=readData("data/popgenome/SUPER_1_vcf_popgenome/",format = "VCF")
# get vector of all position (= @biallelic.site), then get position indexes in each windows which passes the min nb of SNP filter. Position indexes are listed in splitting.data function :
# test_split=splitting.data(testdata, positions = list(which(unlist(testdata@region.data@biallelic.sites) > 200000 & unlist(testdata@region.data@biallelic.sites) < 800000)))
# check for results :
# split_vcf@region.data@biallelic.sites --> list of positions after selected windows
# testdata@region.data@biallelic.sites  --> list of all position in the vcf


# OUTPUT
#---------------------#
# from the stefano's function 'boot_popgenome' itself
# - global FST with bootstrap distribution 
# - pairwise FST with bootstrap distrib

# from the table results outputed by the loop :
# Table of results with Global FST and FST distrib value for 2.5% , 50%, and 97.5% of data
# upper, mid and lower bounds
# nb of SNP in the windows
# Na if no SNP in the windows

# STEP 1 
#---------#
# Gunzip each vcf and put each VCF in individuals folder name 'Chr_vcf_popgenome'

folder_basename="_vcf_popgenome" # create variable name of folder containing each VCF

#------------------------#
# STEP 2 : Run the loop
#------------------------#
# create variable used in loop
SW_Size = 1000000
Overlap_size = 500000
min_SNP=50 # only print windows : "if (length(sss)>min_SNP)" 
bootstrap=5 # 5 for testing but increase the nb of bootstrap for reliable analysis


# Loop on VCF
for (chr in chr_list) {
  
  # read data
  vcf_data=readData(paste("data/popgenome/", chr, folder_basename, sep = ""),format = "VCF")
  
  # get bounds of windows for each Chr (each VCF)
  dim_locus<-SW_Size
  interval<-Overlap_size
  vettore_snps_pos<-unlist(vcf_data@region.data@biallelic.sites)
  low_bound<-seq(1, max(vettore_snps_pos)+interval, interval) # create low bounds vector
  upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval) # create upper bound vector
  ris_fst_sliding=data.frame() # initiate a dataframe
  ris_fst_sliding=rbind(ris_fst_sliding, rep("Na",9)) # intiate 9 columns
  colnames(ris_fst_sliding)=c("global_FST", "2.5%", "50%", "97.5%", "nloci/snps", "mid_wind_GFst", "low_bound_GFst", "upper_bound_GFst", "Nb_SNP_GFst")
  
  # Loop on SW : Fill result table for each window
  for (i in 1:length(low_bound)) {
    sss=which(low_bound[i] < vettore_snps_pos &  upper_bound[i] > vettore_snps_pos) # get position vector includes in the current windows
    if (length(sss)>min_SNP) {
      # run the vcf split according to the windows
      split_vcf=splitting.data(vcf_data, positions = list(which(unlist(vcf_data@region.data@biallelic.sites) > low_bound[i] & unlist(vcf_data@region.data@biallelic.sites) < upper_bound[i])))
      # run stefano function of Global fst
      res=boot_popgenome(split_vcf, lista_pop, bootstrap)
      # add column to fill the table results
      ris_fst_sliding[i,c("global_FST", "2.5%", "50%", "97.5%", "nloci/snps")]=c(res[[1]])
      ris_fst_sliding[i,"mid_wind_GFst"] = (low_bound[i] + upper_bound[i])/2
      ris_fst_sliding[i,"low_bound_GFst"] = low_bound[i]
      ris_fst_sliding[i,"upper_bound_GFst"] = upper_bound[i]
      ris_fst_sliding[i,"Nb_SNP_GFst"] = length(sss)
    }
    
    else{
      # create line with "na" for the windows with nb_snp < threshold
      ris_fst_sliding[i,c("global_FST", "2.5%", "50%", "97.5%", "nloci/snps")] = "Na"
      ris_fst_sliding[i,"mid_wind_GFst"] = (low_bound[i] + upper_bound[i])/2
      ris_fst_sliding[i,"low_bound_GFst"] = low_bound[i]
      ris_fst_sliding[i,"upper_bound_GFst"] = upper_bound[i]
      ris_fst_sliding[i,"Nb_SNP_GFst"] = length(sss)
      
    }
    
  } # end loop on SW
  
  # write table results in local folder
  name_file=paste("Global_FST_SW", chr, sep="_")
  
  path_file=paste("Sliding_global_FST/",name_file, sep="")
  path_file
  write.table(ris_fst_sliding,path_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
} # end loop on chr

# plot golbal FST
#-----------------#

for (chr in chr_list){
  # read FST file from local folder
  name_file=paste("Global_FST_SW", chr, sep="_")
  path_file=paste("Sliding_global_FST/",name_file, sep="")
  global_FST_table=read.table(path_file, header = T, na.strings = "Na")
  
  # create plot name
  plot_name=paste(chr, "Global_FST_", sep="")

  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  coord_gene=gene[which(gene$Chrom == chr),]
  
  # read chr length
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  p=ggplot() + 
  geom_point(data=global_FST_table,aes(x=global_FST_table$mid_wind_GFst,y=global_FST_table$global_FST),pch=1)+
  
  labs(y="Global FST", x = c(chr,"Windows of",SW_Size))+
  ggtitle("FST")+
  
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1, 
                                   size = 12),
        
        text = element_text(size=12),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5)) +
  
  scale_x_continuous(breaks = seq(0, length_chr, 
                                  by = (5000000))) +
  scale_y_continuous(breaks = seq(0, 0.8, 
                                  by = 0.1), limits = c(0,0.8))
  
  # (optional, uncomment if needed) print start and end of gene on chr with line
  # geom_vline(xintercept = coord_gene$start, 
  #            color = "red") +
  # geom_vline(xintercept = coord_gene$end, 
  #            color = "blue")

assign(plot_name, p, envir=parent.frame())

ggsave(path = "Sliding_global_FST/", filename = paste(chr,"global_fst.pdf", sep = "_"), 
       p, device = "pdf", height = 10, width = 12)
}

#===============================================#

#===============================================#
#===============================================#
# ------ FST by SNP by chr ----
#===============================================#
#===============================================#

#----------#
# Infos  : 
#----------#

# INPUT 
#-------#
# VCF with chosen filters : here, 0% Na

# FUNCTION ARGS
#---------------------#
# risultati<-fst.hudson(vcf_name, lista_pop)
# vcf_name : vcf name are vcf read by vcfR save in local variable
# lista_pop = list of populations containing ID in each (list of list type) : created in metadata part

# OUTPUT
#---------------------#
# table of FST values = "pairwise_fst" / "pairwise_fst_nuc" are outputed by the function

# Formated table is outputed by the loop (pairwise FST at each position)
# "POP1.POP2"	"POP1.POP3" "POP1.POP4" "POP2.POP3" "POP2.POP4" "POP3.POP4" "position"


#------------#
# Read VCF : 
#------------#

# Read VCF INPUT
for (chr in chr_list){
  name_vcf=paste(chr, base_name_MAF, sep="")
  path_vcf=paste(c(dir_data,chr,VCF_20_MAF005), collapse = "")
  readvcf=read.vcfR(path_vcf)
  assign(name_vcf, readvcf, envir=parent.frame())
}

#-----------------#
# Run the loop
#-----------------#
for (chr in chr_list){
  
  # get vcf from R env (read by vcfR earlier)
  name_vcf=paste(chr, base_name_MAF, sep="")
  read_vcf=get(name_vcf)
  
  # Run FST + theta sliding windows
  risultati<-fst.hudson(read_vcf, lista_pop)
  
  # name file
  name_file=paste("FST_bySNP", chr, sep="_")
  
  # get positions :
  pos_vcf=getPOS(read_vcf)
  
  # Create and format table
  risultati_fst=as.data.frame(risultati$pairwise_fst_nuc)
  risultati_fst_pos=cbind(risultati_fst, pos_vcf)
  colnames(risultati_fst_pos)=c(colnames(risultati_fst), "position")
  
  assign(name_file,risultati_fst_pos, envir=parent.frame())
  
  # Save FST table
  path_file=paste("FstBySNP/",name_file, sep="")
  write.table(risultati_fst_pos,path_file,quote=F)
}

#----------#
# Plot FST
#----------#
head(fst_by_snp_melt)

for (chr in chr_list){
  # read FST file from local folder
  name_file=paste("FST_bySNP", chr, sep="_")
  path_file=paste("FstBySNP/",name_file, sep="")
  fst_by_snp=read.table(path_file, header = TRUE, na.strings = "NA")
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]
  
  # get size of the current chromosome
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  # Format table : discard negative value
  fst_by_snp[fst_by_snp < 0] = 0  
  
  # melt FST table to plot easly multiple pairwise comparaison in the same graph
  fst_by_snp_melt=melt(fst_by_snp, id.vars = c("position"))
  
  # loop to plot each pairwise comparison in indivisuals plot
  
  for (comp_pop in unique(fst_by_snp_melt$variable)) {

    # get only the current comparison to plot
    fst_by_snp_melt_current=fst_by_snp_melt[which(fst_by_snp_melt$variable == comp_pop),]
    
    # plot
    plot_name=paste(chr,"_",comp_pop ,"_FSTbySNP", sep="")
    
    p=ggplot() + 
      geom_point(data=fst_by_snp_melt_current,aes(x=fst_by_snp_melt_current$position, y=fst_by_snp_melt_current$value,
                                          colour=fst_by_snp_melt_current$variable, 
                                          group=fst_by_snp_melt_current$variable),
                 pch=1)+
      
      labs(y="FST", x = c(chr))+
      ggtitle("FST by snp")+
      
      theme(axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       hjust=1, 
                                       size = 12),
            
            text = element_text(size=12),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "white"),
            panel.background = element_rect(fill = "white",
                                            colour = "black",
                                            size = 0.5)) +
      scale_x_continuous(breaks = seq(0, length_chr, 
                                      by = (5000000))) +
      scale_y_continuous(breaks = seq(0, 1, 
                                      by = 0.1))
    
    # (optional, uncomment if needed) print start and end of gene on chr with line
    #geom_vline(xintercept = coord_gene$start, 
    #           color = "red") +
    #geom_vline(xintercept = coord_gene$end, 
    #           color = "blue") +
    
    assign(plot_name, p, envir=parent.frame())
    
    ggsave(path = "FstBySNP/", filename = paste(chr,comp_pop,"fst_by_snp_scan.pdf", sep = "_"), 
           p, device = "pdf", height = 10, width = 12)
  }
}

#====================#
#====================#
# ------ PCadapt ----
#====================#
#====================#

# INPUT 
#-------#
# See readme in the data/pcadapt folders to generat input for pcadapt

# bed files used for the pcadapt function are generated with plink :
#   - plink --vcf vcf_Na20_MAF005.vcf.gz" --double-id --allow-extra-chr --make-bed --out $vcf_Na20_MAF005.bed" in bash script

# positions vector from the input vcf
#   - bcftools view -H vcf_Na20_MAF005.vcf.gz | cut -f1,2 >> WGS_21_${chr}_GATK_TAG_Flowqual_Noindels_Norepeat_SNP_order_15Na20_MAF005.positions

# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# FUNCTION ARGS
#---------------------#
# PCadapt function  : pcadapt(input = read.pcadapt(path_bed,type="bed"),min.maf=0.01, K=2)
# input  = bed file
# min.maf 
# K

# OUTPUT
#---------------------#
#  write.table(Pval_all, pval_path,  quote=F) # pval files include "CHR","POS","pvalue","qvalue"
#  write.table(scores, scores_path ,quote=F) # scores file include PCA scores
#  write.table(outliers, outliers_path, quote=F) # outliers files include only pos which passe the alpha filter

# To adapt to your data in the loop
# name_pos / var_pos
# path_bed
# name_scores / name_pval / outliers_table


# Read position vector
#---------------------#
# position files are present in .txt format in "dir_data" with one column without header :

for (chr in chr_list){
  path_pos=paste(dir_data,chr,position_20_MAF, sep = "")
  readpos=read.table(path_pos, header = FALSE)
  var_pos=paste(chr,base_name_position_20_MAF, sep="_")
  assign(var_pos, readpos, envir=parent.frame())
}

# Set variables :
#---------------#
alpha <- 0.05 # put the alpha probability threshold

# Run PCAdapt
#-----------#
for (chr in chr_list){
  # get bed file
  path_bed=paste(c(dir_data, "pcadapt/",chr,BED_20_MAF005), collapse = "")
  
  # get position vector
  positions=get(paste(chr, base_name_position_20_MAF, sep="_"))

  # run pcadapt
  respca_10 = pcadapt(input = read.pcadapt(path_bed,type="bed"),min.maf=0.05, K=2)
  
  # format
  scores = data.frame(respca_10$scores)
  rownames(scores)= samples #vector of pop , sex , etc .. in same order as in the bed file
  popCol = num2col(as.color(pop))
  names(popCol) = rownames(scores)
  scores$col = popCol[rownames(scores)]
  
  # Get component variance explanation
  PEV = paste(round((respca_10$singular.values^2)*100,digits=1),"%",sep="")

  # Get qval for each pvalue
  respca_10$qval <- qvalue(respca_10$pvalues)$qvalues

  # get only pvalue (SNP) for which the Pvalue is above the threshold
  outliers <- which(respca_10$qval < alpha)
  print(length(outliers))

  # Create a table with all pvalue and add absolute position in each snp
  Pval_all = cbind(positions,respca_10$pvalues,respca_10$qval)
  colnames(Pval_all)= c("CHR","POS","pvalue","qvalue")
  # Save table in local env
  name_scores=paste(chr, ".scores", sep = "")
  scores_path=paste("pcadapt/", name_scores, sep = "")


  name_pval=paste(chr, ".pval", sep = "")
  pval_path=paste("pcadapt/", name_pval, sep = "")


  outliers_table=paste(chr, ".outliers", sep = "")
  outliers_path=paste("pcadapt/", outliers_table, sep = "")

  # Write table in local folder
  write.table(Pval_all, pval_path,  quote=F)
  write.table(scores, scores_path ,quote=F)
  write.table(outliers, outliers_path, quote=F)

}

# Plot Pval
#-----------#

# put qval threshold
qval=0.0001

# RUN Loop for plot
for (chr in chr_list){
  
  # read FST file from local folder
  pval_file=read.table(paste("pcadapt/", chr, ".pval", sep = ""), header = TRUE)
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]
  
  # get current chr length
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
    
  # plot
  p=ggplot() + 
    geom_point(data=pval_file,aes(x=pval_file$POS,
                                         y=-log10(pval_file$qvalue)),pch=19, color = "green")+
    
    labs(y="-log10(qvalue)", x = "Positions") +
    ggtitle(paste(chr, " PCAdapt", sep = ""))+
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1, 
                                     size = 12),
          
          text = element_text(size=12),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
    
    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000))) +
    scale_y_continuous(breaks = seq(0, 250, 
                                    by = 10), limits = c(0,250)) +

	# (optional, uncomment if needed) print start and end of gene on chr with line
    # geom_vline(xintercept = coord_gene$start, 
    #            color = "red") +
    # geom_vline(xintercept = coord_gene$end, 
    #            color = "blue") +
	
    geom_hline(yintercept = -log10(qval), 
               color = "black")
  
  ggsave(path = "pcadapt/", width = 12, height = 10, filename = paste(chr, ".pdf", sep = ""), p, device = "pdf")
}

#========================================#
#========================================#
# ------ Detect fixed genotype in pop ----
#========================================#
#========================================#

#-------#
# INFOS
#-------#

# INPUT 
#-------#
# VCF with chosen filters : 0% Na (because computing a proportion of 0/1 vs 0/0 no Na allowed)

# FUNCTION ARGS
#---------------------#
# geno_table=sliding_SNP_homo_fix(vcf_name, lista_pop, window=SW_Size, slide = Overlap_size)
# vcf_name = variable name of vcf saved in R env
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part
# SW_Size=100000 : Window size
# slide=50000 : slide size
  
# OUTPUT
#---------------------#
# table with the nb of fixed sites by windows
# "low_bound", "upper_bound", "mid_wind", "Nb_SNP", "Nb_Fixed_Sites"

#------------#
# Read VCF : 
#------------#
# read vcf with VCFR

for (chr in chr_list){
  name_vcf=paste(base_name_0, chr, sep="_")
  path_vcf=paste(c(dir_data,chr,VCF_0), collapse = "")
  readvcf=read.vcfR(path_vcf)
  assign(name_vcf, readvcf, envir=parent.frame())
}

#-----------------#
# Run the loop
#-----------------#
SW_Size=1000000
Overlap_size=500000

for (chr in chr_list){
  
  # read VCF, position and variable
  vcf_name=get(paste(base_name_0, chr, sep="_"))

  # get Ratio 
  geno_table=sliding_SNP_homo_fix(vcf_name, lista_pop, window=SW_Size, slide = Overlap_size)
  geno_table=as.data.frame(geno_table)
  
  # generate table in good format
  
  # get number of SNP tot in interval
  vettore_snps_pos<-getPOS(vcf_name)
  low_bound<-geno_table$per_plot-(Overlap_size+0.5)
  upper_bound<-geno_table$per_plot+(Overlap_size-0.5)
  Nb_SNP_vec=c()
  
  # create the final table : 
  for (i in 1:length(low_bound)) {
    Nb_SNP_vec_current=length(which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] ))
    if (Nb_SNP_vec_current > 0){Nb_SNP_vec=c(Nb_SNP_vec, Nb_SNP_vec_current)} # set if SNP >=0 because windows boundaries are extracted from geno_table which contain only SW with SNP > 0
    else{} # do not print windows in the SNP vector if there is no SNP in the windows
  }
  
  # Add windows and pos to final 
  geno_table$low_bound=geno_table$per_plot-(Overlap_size+0.5)
  geno_table$upper_bound=geno_table$per_plot+(Overlap_size-0.5)
  geno_table$Nb_SNP_tot=Nb_SNP_vec
  geno_table_range=geno_table[,c(3,4,1,5,2)]
  colnames(geno_table_range)=c("low_bound_Fix", "upper_bound_Fix", "mid_wind_Fix", "Nb_SNP", "Nb_Fixed_Sites")
  
  # write table
  write.table(geno_table_range,paste("freq_genotype/", chr,"Sliding_geno_count",sep=""),quote=F)
}

#--------#
# PLOT
#--------#
for (chr in chr_list){

  # read file from local folder
  geno_file=read.table(paste("freq_genotype/",  chr,"Sliding_geno_count",sep=""))

  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]

  # get size of the current chromosome
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]

  # plot
  plot_FixedSNP=ggplot() +
    geom_point(data=geno_file,
               aes(x=geno_file$mid_wind_Fix,y=geno_file$Nb_Fixed_Sites,
                   color = "red"),pch=19) +

    ylab("Nb Fixed Sites per Windows")  +
    ggtitle(paste(chr, "Fixed Sites Between pops", sep = " "))+

    # scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 0.5)) +

    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust=1,
                                     size = 12),

          text = element_text(size=12),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +

    scale_x_continuous(breaks = seq(0, length_chr,
                                    by = (5000000)))
									
	# (optional, uncomment if needed) print start and end of gene on chr with line
    #geom_vline(xintercept = coord_gene$start, 
    #           color = "red") +
    #geom_vline(xintercept = coord_gene$end, 
    #           color = "blue") +

  ggsave(path = "freq_genotype/", filename = paste(chr, "SW_FixedSNP.pdf", sep = "_"), plot_FixedSNP, device = "pdf", width = 20, height = 10)

}

#========================================#
#========================================#
# ------ Diversity on sliding windows ----
#========================================#
#========================================#

#------------#
# INFOS
#------------#

# INPUT 
#-------#
# VCF with chosen filters : 0% Na and heterozygous filters (SFS-based analyses)
# position vector : total position in the chr
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# FUNCTION ARGS
#------------------#
# SlidingSFS <- sliding_spectre_onepop_stand(VCf_pop, SW_Size, Overlap_size, positions_vector) :
# VCf_pop : VCF
# SW_Size=100000 : Window size
# slide=50000 : slide size
# positions_vector : numeric vector of position

# SlidingDiv <- calcola_TD_folded(folded_sfs)
# folded_sfs : folded SFS, here got from SlidingSFS results 

# Output :
#--------#
# 1. SFS files are created for each pop + all pop by windows : see Sliding_SFS_n files

# 2. Diversity files are created for each pop + all pop by windows : See Sliding_DiversityIndices_n files
# colnames : "Windows" "TajimasD" "Thetapi" "ThetaS" "S" "Thetapi_persite" "ThetaS_persite"

# 3. Sliding_TD_allpops : TD for all pop 
# colanmes : low_bound_theta upper_bound_theta mid_wind_theta Nmb_Positions_theta Nb_SNP_theta pop1 pop2 pop3 pop4 pop1+pop2+pop3+pop4

# 4. Sliding_theta_norm_allpops : theta for all pop 
# colnames : low_bound_TD upper_bound_TD mid_wind_TD Nmb_Positions_TD Nb_SNP_TD pop1 pop2 pop3 pop4 pop1+pop2+pop3+pop4

#------------#
# Read VCF : 
#------------#
# Filters : Na  = 0% and basics filters (no indel, no repeat, no Lowqual) + Rate het max = 80%

for (chr in chr_list){
  name_vcf=paste(chr, base_name_0, sep="_")
  path_vcf=paste(c(dir_data,chr,VCF_0), collapse = "")
  readvcf=read.vcfR(path_vcf)
  assign(name_vcf, readvcf, envir=parent.frame())
}

# Read position vector
#----------------------#
for (chr in chr_list){
  name_pos=paste(c(chr,position_0), collapse = "")
  path_pos=paste(dir_data,name_pos,sep = "")
  readpos=read.table(path_pos, header = FALSE)
  pos_vector=readpos$V2
  var_pos=paste(chr,base_name_position_0, sep="_")
  assign(var_pos, pos_vector, envir=parent.frame())
  
}

#-----------------#
# Run the loop
#-----------------#
SW_Size=1000000
Overlap_size=500000

# created in metadata part
lista_pop_all 

# initiate variable
Theta_norm_all_pop = c()
TD_all_pop = c()


# Run div by pop and by chr
for (chr in chr_list){
  
  # initiate variable
  Theta_norm_all_pop = c()
  TD_all_pop = c()
  
  # read VCF, position and variable
  vcf_name=get(paste( chr, base_name_0 ,sep="_"))
  positions_vector=get(paste(chr,base_name_position_0, sep="_"))
  dir(dir_data)
  # get div by pop defined in lista_pop
  for (pop in 1:length(lista_pop_all)) {
    print(pop)
    # initiate pop nb ind
    n_ind = length(lista_pop_all[[pop]])+2
    
    # subset vcf
    VCf_pop<-vcf_name[,c("FORMAT",lista_pop_all[[pop]])]
    
    # get SFS in SW
    SlidingSFS <-sliding_spectre_onepop_stand(VCf_pop, SW_Size, Overlap_size, positions_vector)
    SlidingSFS=as.data.frame(SlidingSFS)
    
    # Get diversity index from SFS with "calcola_TD_folded" function
    # Return Table with 4 column
    SlidingDiv = apply(SlidingSFS[,3:n_ind], 1,calcola_TD_folded)
    
    # Format diverity table 
    SlidingDiv_table = matrix(unlist(SlidingDiv), ncol=4, nrow=nrow(SlidingSFS), byrow=T)
    
    # Compute Thetapi and theta normalized 
    thetapi=SlidingDiv_table[,2]/SlidingSFS[,2]
    thetas=SlidingDiv_table[,3]/SlidingSFS[,2] # 
    
    # Add in SW SFS , all Div table , thetapi and theta
    SlidingDiv_table2 <- cbind(SlidingSFS[,1], SlidingDiv_table,thetapi,thetas)
    
    # Format table
    Theta_norm_all_pop = cbind(Theta_norm_all_pop,thetapi)
    TD_all_pop = cbind(TD_all_pop,SlidingDiv_table2[,2])
    colnames(SlidingDiv_table2) = c("Windows","TajimasD","Thetapi","ThetaS","S","Thetapi_persite","ThetaS_persite")
    SlidingDiv_table2=as.data.frame(SlidingDiv_table2)
    # write table
    write.table(SlidingDiv_table2, paste("Sliding_div/", chr,"Sliding_DiversityIndices_",pop,".txt",sep=""))
    write.table(SlidingSFS,paste("Sliding_div/", chr, "Sliding_SFS_",pop,".txt",sep=""),quote=F)
  }
  
  # get number of SNP tot in interval
  vettore_snps_pos<-getPOS(vcf_name)
  low_bound<-SlidingSFS$per_plot-(Overlap_size+0.5)
  upper_bound<-SlidingSFS$per_plot+(Overlap_size-0.5)
  Nb_SNP_vec=c()
  for (i in 1:length(low_bound)) {
    Nb_SNP_vec_current=length(which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] ))
    if (Nb_SNP_vec_current > 0){Nb_SNP_vec=c(Nb_SNP_vec, Nb_SNP_vec_current)} # In 'sliding_spectre_onepop_stand' function only print windows : 'if (length(sss)>0)'. Same for SNP_vec
    else{} # do not print windows in the SNP vector if there is no SNP in the windows
  }

  # Add windows and pos to final thetapi
  Theta_norm_all_pop_toplot = as.data.frame(cbind(SlidingSFS[,1:2],Theta_norm_all_pop))
  colnames(Theta_norm_all_pop_toplot) = c("mid_wind","Nmb_Positions",names(lista_pop_all))
  
  Theta_norm_all_pop_toplot$low_bound=SlidingSFS$per_plot-(Overlap_size+0.5)
  Theta_norm_all_pop_toplot$upper_bound=SlidingSFS$per_plot+(Overlap_size-0.5)
  Theta_norm_all_pop_toplot$Nb_SNP_theta=Nb_SNP_vec
  
  Theta_norm_all_pop_toplot = Theta_norm_all_pop_toplot[,c("low_bound","upper_bound","mid_wind","Nmb_Positions","Nb_SNP_theta", names(lista_pop_all))]
  colnames(Theta_norm_all_pop_toplot)=c("low_bound_theta","upper_bound_theta","mid_wind_theta","Nmb_Positions_theta","Nb_SNP_theta", names(lista_pop_all))
  
  # Add windows, nb pos to final TD table 
  TD_all_pop_toplot = as.data.frame(cbind(SlidingSFS[,1:2],TD_all_pop))
  colnames(TD_all_pop_toplot) = c("mid_wind","Nmb_Positions",names(lista_pop_all))
  
  TD_all_pop_toplot$low_bound=SlidingSFS$per_plot-(Overlap_size+0.5)
  TD_all_pop_toplot$upper_bound=SlidingSFS$per_plot+(Overlap_size-0.5)
  TD_all_pop_toplot$Nb_SNP_TD=Nb_SNP_vec


  TD_all_pop_toplot = TD_all_pop_toplot[,c("low_bound","upper_bound","mid_wind","Nmb_Positions","Nb_SNP_TD", names(lista_pop_all))]
  colnames(TD_all_pop_toplot)=c("low_bound_TD","upper_bound_TD","mid_wind_TD","Nmb_Positions_TD","Nb_SNP_TD", names(lista_pop_all))
  # write final thetapi and TD tables within each pop for the current chr
  write.table(Theta_norm_all_pop_toplot,paste("Sliding_div/", chr,"Sliding_theta_norm_allpops",".txt",sep=""),quote=F)
  write.table(TD_all_pop_toplot,paste("Sliding_div/", chr,"Sliding_TD_allpops",".txt",sep=""),quote=F)
  
}

#-----------------#
# Plot div index
#-----------------#
# input table
#--------------#
# head(Theta_norm_all_pop_toplot)
#           low_bound_TD upper_bound_TD mid_wind_TD Nmb_Positions_TD Nb_SNP_TD pop1   pop2  pop3   pop4   pop1+pop2+pop3+pop4
# temp_spectre 0          1e+07         5000000.5   719               719      0.508  0.130 0.344  0.217    0.551 
# temp_spectre.1 5e+06 1.5e+07 10000000.5 854 854 0.395  -0.0789  0.288  0.252 0.318 

# head(TD_all_pop_toplot)
#                   low_bound_TD upper_bound_TD mid_wind_TD Nmb_Positions_TD Nb_SNP_TD    pop1           pop2       pop3        pop4          pop1+pop2+pop3+pop4
# temp_spectre      0.00e+00       1.00e+07    5.00e+06     1363            1363        -0.215731867 -0.064029279  0.036227171 -0.30981514        -0.216637752
# temp_spectre.1    5.00e+06       1.50e+07    1.00e+07     1606            1606        -0.217699898  0.002703841 -0.039886432 -0.27042018        -0.193465700



# PLOT DIV
#--------------#
for (chr in chr_list){

  # read FST file from local folder
  theta_file=read.table(paste("Sliding_div/", chr,"Sliding_theta_norm_allpops",".txt",sep=""))
  theta_file_all=melt(theta_file[,-c(1,2,4,5)],id=c("mid_wind_theta"))

  TD_file=read.table(paste("Sliding_div/", chr,"Sliding_TD_allpops",".txt",sep=""))
  TD_melt=melt(TD_file[,-c(1,2,4,5)],id=c("mid_wind_TD"))

  # plot
  plot_theta=paste(chr, "_theta", sep="")
  plot_TD=paste(chr, "_TD", sep="")

  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]
  
  # current chr length 
  length_chr=table_chr[which(table_chr$Chr==chr),"length"]
  
  # plot theta
  plot_theta=ggplot() +
    geom_point(data=theta_file_all,
               aes(x=theta_file_all$mid_wind_theta,y=theta_file_all$value,
                   group=theta_file_all$variable,
                   color = variable),pch=1) +

    ylab("ThetaPi per site")  +
    ggtitle(paste(chr, "theta", sep = " "))+

   theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1, 
                                     size = 12),
          
          text = element_text(size=12),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
						
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +

    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000)))
    # (optional, uncomment if needed) print start and end of gene on chr with line
    # geom_vline(xintercept = coord_gene$start, 
    #            color = "red") +
    # geom_vline(xintercept = coord_gene$end, 
    #            color = "blue") +
    
  # plot TD
  plot_TD=ggplot() +
     geom_point(data=TD_melt,
                aes(x=TD_melt$mid_wind,y=TD_melt$value,
                group=theta_file_all$variable,color = variable),pch=1) +

     ylab("Tajima's D")  +
     ggtitle(paste(chr, "TD", sep = " "))+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
           text = element_text(size = 15)) +
    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000)))
    # (optional, uncomment if needed) print start and end of gene on chr with line
    # geom_vline(xintercept = coord_gene$start, 
    #            color = "red") +
    # geom_vline(xintercept = coord_gene$end, 
    #            color = "blue") +

  ggsave(path = "Sliding_div/", filename = paste(chr, "theta.pdf", sep = "_"), plot_theta, device = "pdf", width = 20, height = 10)
  ggsave(path = "Sliding_div/", filename = paste(chr, "TD.pdf", sep = "_"), plot_TD, device = "pdf", width = 20, height = 10)

}

#=======================================#
#=======================================#
# ------ GLOBAL SFS by Chromosome ----
#=======================================#
#=======================================#

# Global SFS by chr is needed to compute the SFS euclidian distancies between SFS in SW

# INPUT 
#-------#
# VCF with chosen filters : 0% Na and heterozygous filters (SFS-based analyses)
# position vector : total position in the chr
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part

# FUNCTION ARGS
#------------------#
# sfs_tot = SFS get with DNABin then site.spectrim function (R package)
# SlidingSFS <-calcola_normalized_foldedSFS(sfs_tot)

# Normalized by the spectrum by the number of sites 
# Stefano s function: calcola_normalized_foldedSFS
# ARGUMENT : Folded SFS. 
# Methode : - Divide the SFS vector by the rank of frequencies (sum of SFS). 
#           - asse_x = vector of number of ind / nb_tot_ind
# Plot the obtained SFS curve 
# Remember that the standardized of a constant size pop is parallel to the x-axis, 

# OUTPUT
#---------#
# norm_sfs : SFS normalized for the whole chromosome (write in local file not implemented yet)
# Plot of SFS by chromosome

#-----------------#
# read VCF input
#-----------------#

for (chr in chr_list){
  name_vcf=paste(chr, base_name_0, sep="_")
  path_vcf=paste(c(dir_data,chr,VCF_0), collapse = "")
  readvcf=read.vcfR(path_vcf)
  assign(name_vcf, readvcf, envir=parent.frame())
}

#----------------------#
# read position input
#----------------------#
for (chr in chr_list){
  name_pos=paste(c(chr,position_0), collapse = "")
  path_pos=paste(dir_data,name_pos,sep = "")
  readpos=read.table(path_pos, header = FALSE)
  pos_vector=readpos$V2
  var_pos=paste(chr,base_name_position_0, sep="_")
  assign(var_pos, pos_vector, envir=parent.frame())
  
}

#-------------------------------#
# Run SFS in total chromosome
#-------------------------------#
for (chr in chr_list){
  vcf_name=get(paste(chr, base_name_0, sep="_"))
  positions_vector=get(paste(chr,base_name_position_0, sep="_"))
  
  # Get DNA sequences
  DNA_bin <-vcfR2DNAbin(vcf_name,  extract.indels = TRUE,
                        consensus = FALSE,
                        extract.haps = F,
                        unphased_as_NA = F,
                        asterisk_as_del = FALSE,
                        ref.seq = NULL,
                        start.pos = NULL,
                        verbose = TRUE)
  
  # Compute folded site frequency spectrum with the function site.spectrum (library pegas required)
  sfs_tot<-site.spectrum(DNA_bin, folded = TRUE)
  sfs_tot<-as.vector(sfs_tot)
  
  # Run SFS normalization
  norm_sfs<-calcola_normalized_foldedSFS(sfs_tot)
  sfs_norm_name=paste(chr, "_norm_sfs", sep = "")
  assign(sfs_norm_name, norm_sfs, envir=parent.frame())
  
  # Plot each SFS
  plot_name=paste(chr, "_SFS_total", sep="")
  norm_sfs=as.data.frame(norm_sfs)
  write.table(norm_sfs, paste("SlidingSFS/", chr, "SFS_all.table", sep=""))
  
  # plot
  p=ggplot()+
    geom_line(data=norm_sfs, aes(x=norm_sfs$asse_x, y=norm_sfs$eta_2)) +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1, 
                                     size = 12),
          
          text = element_text(size=12),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5)) +
										  
    scale_y_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1))
  
  assign(plot_name, p, envir=parent.frame())
  ggsave(path = "SlidingSFS/", filename = paste(chr, "SFS_tot.pdf", sep = "_"), p, width = 10, height = 5, device = "pdf")
  
}
#======================================#
#======================================#
# ------ SFS in Sliding windows ----
#======================================#
#======================================#

# INPUT 
#-------#
# norm_sfs : SFS normalized for the whole chromosome  (See "GLOBAL SFS by Chromosome" section)
# VCF with chosen filters : 0% Na and heterozygous filters
# position vector : total position in the chr
# lista_pop = list of populations containing ID in each (list of list type) : create in metadata part


# FUNCTION ARGS
#------------------#
# risultati_spectre_slide = sliding_spectre_onepop_stand(vcf_name, window, slide, positions_vector)
# SW_Size=100000 : Window size
# slide=50000 : slide size
 
# SlidingSFS <-calcola_normalized_foldedSFS(sfs_tot)
# sfs_tot = Global folded SFS

# dist_euclid = calcola_dist_spettro<-function(data, vettore_average)
# data = SlidingSFS
# vettore_average = norm_sfs
# Compute the euclidian distancies between the total SFS computed on the whole chr and the local SFS computed in SW

# OUTPUT
#---------#
# table = "low_bound","upper_bound","mid_wind","Nmb_Positions", "Nb_SNP_SFS", "dist_eucli"
# plot of distancies along chr

# intiate variable 
window<-1000000
slide<-500000

#---------------------------#
# Run euclidian distancies
#---------------------------#
for (chr in chr_list){ 
  # read VCF already loaded
  vcf_name=get(paste(chr, base_name_0, sep="_"))
  sfs_norm=get(paste(chr, "_norm_sfs", sep = ""))
  positions_vector=get(paste(chr,base_name_position_0, sep="_"))
  
  # get SFS on sliding windows
  risultati_spectre_slide<-sliding_spectre_onepop_stand(vcf_name, window, slide, positions_vector)
  # normalize SFS   
  risultati_spectre_slide_norm<-t(apply(risultati_spectre_slide[,3:23],1,calcola_normalized_foldedSFS)) 
  # get only column of SFS normalized
  solo_norm<-risultati_spectre_slide_norm[,22:42]
  # compare with the average SFS 
  calcola_dist_spettro<-function(data, vettore_average){
    distanza<-sum((data-vettore_average)^2)
    return(distanza)
  }
  
  dist_euclid<-apply(solo_norm, 1, calcola_dist_spettro, sfs_norm)
  dist_euclid
  
  # get number of SNP tot in interval
  vettore_snps_pos<-getPOS(vcf_name)
  low_bound<-risultati_spectre_slide[,c(1)]-(slide+0.5)
  upper_bound<-risultati_spectre_slide[,c(1)]+(slide-0.5)
  Nb_SNP_vec=c()
  for (i in 1:length(low_bound)) {
    Nb_SNP_vec_current=length(which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] ))
    if (Nb_SNP_vec_current > 0){Nb_SNP_vec=c(Nb_SNP_vec, Nb_SNP_vec_current)} # In 'sliding_spectre_onepop_stand' function only print windows : 'if (length(sss)>0)'. Same for SNP_vec.
    else{} # do not print in the SNP vector if there is no SNP in the windows
  }
  
  # get table 
  table_graphe <- as.data.frame(risultati_spectre_slide[,c(1,2)])
  table_graphe$nb_SNP <- Nb_SNP_vec
  table_graphe$dist_eucli <- dist_euclid
  table_graphe$low_bound <- risultati_spectre_slide[,c(1)]-(slide+0.5)
  table_graphe$upper_bound <- risultati_spectre_slide[,c(1)]+(slide-0.5)
  table_graphe=table_graphe[,c(5,6,1,2,3,4)]
  colnames(table_graphe)<- c("low_bound_SFS","upper_bound_SFS","mid_wind_SFS","Nmb_Positions_SFS", "Nb_SNP_SFS", "dist_eucli")
  write.table(table_graphe, paste("SlidingSFS/", chr, "_dist_SFS.table", sep=""))
  
}

#-----------------#
# plot distancies
#-----------------#
for (chr in chr_list){ 
  # read table from local path
  table_graphe=read.table(paste("SlidingSFS/", chr, "_dist_SFS.table", sep=""))
  
  # (optional, uncomment if needed) read gene : Chrom,start,end,Gene_name
  # coord_gene=gene[which(gene$Chrom == chr),]
  
  # plot table
  graphe1 <- ggplot(data=table_graphe,
                    aes(x=mid_wind_SFS , y=dist_eucli),
                    pch=16) +
    geom_point()+
	
	theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1, 
                                     size = 12),
          
    text = element_text(size=12),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"),
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5)) +
									
	    scale_x_continuous(breaks = seq(0, length_chr, 
                                    by = (5000000)))
  
    # (optional, uncomment if needed) print start and end of gene on chr with line
    # geom_vline(xintercept = coord_gene$start, 
    #            color = "red") +
    # geom_vline(xintercept = coord_gene$end, 
    #            color = "blue") +
			   										  
  ggtitle(paste("SFS",chr))
  
  ggsave(paste("SFS_dist",chr, ".pdf"), graphe1, width = 11, height = 8, path = "SlidingSFS/")
}
