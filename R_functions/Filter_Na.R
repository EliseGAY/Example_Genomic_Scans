###############################################################
##### Filter Individues and Position by heterozygous rate #####
###############################################################

# Library
#----------#
library("vcfR")

# input
#----------#
# A VCF file

# arg
#-------#
# rate_na_max_POS (double) : Max percent of missing data site in one position
# rate_na_max_Ind (double) : Max percent of missing data site in one individus

# output
#-------#
# VCF file filtered for missing data rate by Individu and then by position

# Usage :
#--------#
# Filter_VCF=Filters_Het_Ind_Pos(vcf_file, rate_na_max_POS,rate_na_max_Ind)

# Function 
#-------------------------------------------------------------------------------------------------------------------------------------------#
# Get genotype from vfc with the vcfR package
# Count number of "./." per column and compute na rate (percent) by divided by the number of line (position)
# Create new vcf with only samples which pass the filter
# On the new VCF : Count number of "./." per line and compute heterozygous rate (percent) by divided by the number of column (individus)
# Create new vcf with only samples AND position which pass the filter
#-------------------------------------------------------------------------------------------------------------------------------------------#

Filters_na_Ind_Pos<-function(vfc,rate_na_max_POS,rate_na_max_Ind){
	# Filters_Het_Ind_Pos(vcf_file, rate_na_max_POS,rate_na_max_Ind)
	# Get genotypes data frame
	gt<-extract.gt(vfc,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	# get ind rate het (column)
	somme_na_ind = colSums(gt=="./.")
	rate_na_ind = (somme_na_ind/dim(gt)[1])*100.0
	samples_kept<-names(which(rate_na_ind<=rate_na_max_Ind))
	vcf_filtered_ind<-vfc[,c("FORMAT",samples_kept)]

	# get pos rate het (lines)
	gt_1<-extract.gt(vcf_filtered_ind,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	somme_na_Pos = rowSums(gt_1=="./.")
	rate_na_Pos = (somme_na_Pos/dim(gt_1)[2])*100
	Pos_kept<-which(rate_na_Pos<=rate_na_max_POS)
	vcf_filtered_na<-vcf_filtered_ind[c(Pos_kept),]
	
	return(vcf_filtered_na)
}

#------------#
# Example :
#------------#
