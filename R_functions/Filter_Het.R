# ##############################################################
# #### Filter Individues and Position by heterozygous rate #####
# ##############################################################

# Library
#----------#
library("vcfR")

# input
# ----------#
# A VCF file

# arg
# -------#
# rate_het_max_POS (double) : Max percent of heterozygous site in one position
# rate_het_max_Ind (double) : Max percent of heterozygous site in one individus
# rate_het_min_pos (int) : Min nb of heterozygous site in one position

# output
# -------#
# VCF file filtered for heterozygous rate by Individu and then by position

# Usage :
# --------#
# Filter_VCF=Filters_Het_Ind_Pos(vcf_file, rate_het_max_POS,rate_het_max_Ind)

# Function 
# -------------------------------------------------------------------------------------------------------------------------------------------#
# Get genotype from vfc with the vcfR package
# Count number of "0/1" per column and compute heterozygous rate (percent) by divided by the number of line (position)
# Create new vcf with only samples which pass the filter
# On the new VCF : Count number of "0/1" per line and compute heterozygous rate (percent) by divided by the number of column (individus)
# Create new vcf with only samples AND position which pass the filter
# -------------------------------------------------------------------------------------------------------------------------------------------#

Filters_Het_Ind_Pos<-function(vfc, rate_het_max_POS,rate_het_max_Ind, rate_het_min_pos){
	# Filters_Het_Ind_Pos(vcf_file, rate_het_max_POS,rate_het_max_Ind)
	# Get genotypes data frame
	gt<-extract.gt(vfc,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	# get ind rate het (column)
	somme_het_ind = colSums(gt=="0/1")
	rate_het_ind = (somme_het_ind/dim(gt)[1])*100
	samples_kept<-names(which(rate_het_ind<=rate_het_max_Ind))
	vcf_filtered_ind<-vfc[,c("FORMAT",samples_kept)]

	# get pos nb het min (lines)
	gt_1<-extract.gt(vcf_filtered_ind,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	somme_het_Pos = rowSums(gt_1=="0/1")
	Pos_kept<-which(somme_het_Pos>=rate_het_min_pos)
	vcf_filtered_Het_min<-vcf_filtered_ind[c(Pos_kept),]
	
	
	# get pos rate het (lines)
	gt_1<-extract.gt(vcf_filtered_Het_min,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	somme_het_Pos = rowSums(gt_1=="0/1")
	rate_het_Pos = (somme_het_Pos/dim(gt_1)[2])*100
	Pos_kept<-which(rate_het_Pos<=rate_het_max_POS)
	vcf_filtered_Het<-vcf_filtered_Het_min[c(Pos_kept),]
	
	return(vcf_filtered_Het)
}

# ------------#
# Example :
# ------------#
# vcf_F_DP_NoNa_GWS
#
# ***** Object of Class vcfR *****
# 219 samples
# 30 CHROMs
# 4,110 variants
# Object size: 14.2 Mb
# 0 percent missing data
# *****        *****         *****
# '

## Get genotypes data frame
# gt<-extract.gt(vcf_F_DP_NoNa_GWS,element="GT",mask=F,as.numeric=F,return.alleles = F,
               # convertNA = F,extract = T)

##get ind rate het (column)
# somme_het_ind = colSums(gt=="0/1")
# summary(somme_het_ind)
	"
	   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	   # 6.00   12.00   15.00   32.19   18.00 1479.00 
	"

# rate_het_ind = (somme_het_ind/dim(gt)[1])*100
# summary(rate_het_ind)
	"
	   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	 0.1460  0.2920  0.3650  0.7833  0.4380 35.9854 
	"
# samples_kept<-names(which(rate_het_ind<10.0))
# summary(samples_kept)
	"
	   # Length     Class      Mode 
		  # 217 character character 
	"
# vcf_filtered_ind<-vcf_F_DP_NoNa_GWS[,c("FORMAT",samples_kept)]
# vcf_filtered_ind
	"
	# ***** Object of Class vcfR *****
	# 217 samples
	# 30 CHROMs
	# 4,110 variants
	# Object size: 13.9 Mb
	# 0 percent missing data
	# *****        *****         *****
	"

# # get Pos rate het (lines) from vcf filtered by ind
# gt<-extract.gt(vcf_filtered_ind,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)

# somme_het_Pos = rowSums(gt=="0/1")
# summary(somme_het_Pos)
	"
	   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	   # 0.00    0.00    0.00    1.08    0.00  152.00 
	"
# rate_het_Pos = (somme_het_Pos/dim(gt)[2])*100
# summary(rate_het_Pos)
	"
	   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	 # 0.0000  0.0000  0.0000  0.4978  0.0000 70.0461 
	"
# Pos_kept<-which(rate_het_Pos<50)
# summary(Pos_kept)
# length(Pos_kept)
	"
	# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	# 1    1030    2059    2058    3084    4110

	# > length(Pos_kept)
	# [1] 4103
	"
# vcf_filtered_Het<-vcf_filtered_ind[c(Pos_kept),]
# vcf_filtered_Het

	"
	# ***** Object of Class vcfR *****
	# 217 samples
	# 30 CHROMs
	# 4,103 variants
	# Object size: 13.8 Mb
	# 0 percent missing data
	# *****        *****         *****
	"
