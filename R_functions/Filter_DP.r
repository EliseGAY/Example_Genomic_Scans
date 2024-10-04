###############################################################
##### Filter Individues and Position by sequencing depth #####
###############################################################

# Library
#----------#
library("vcfR")

# input
#----------#
# A VCF file

# arg
#-------#
# rate_min_cover_pos (double) : Max percent of missing data site in one position
# rate_max_cover_pos (double) : Max percent of missing data site in one position
# rate_min_cover_ind (double) : Max percent of missing data site by idividuals
# rate_max_cover_ind (double) : Max percent of missing data site by individuals

# output
#-------#
# VCF file filtered for sequencing depth rate by Individu and then by position

# Usage :
#--------#
# From a VCF file, filters individuals and position according to DP threshold
# Filter_VCF=fun_filtre_cover(vcf_file, rate_na_max_POS,rate_na_max_Ind)

# Function 
#-------------------------------------------------------------------------------------------------------------------------------------------#
# Get sequencing depth from vfc with the vcfR package
# apply mean sequencing depth by individuals (on column)
# Create new vcf with only samples which pass the filter
# On the new VCF : apply mean sequencing depth by position (on row)
# Create new vcf with only samples AND position which pass the DP filter
#-------------------------------------------------------------------------------------------------------------------------------------------#


#########################
# fonction filtre cover #
#########################

fun_filtre_cover<-function(vcf,rate_min_cover_pos,rate_max_cover_pos,rate_min_cover_ind,rate_max_cover_ind){

# get ind DP mean (column)
dpi<- extract.gt(vcf, element='DP', as.numeric = TRUE) 
moyenne_cover_par_indiv = apply(dpi,2,mean,na.rm=T)
# filter individuals according to the rate specified in arg
noms<-names(which(moyenne_cover_par_indiv>=rate_min_cover_ind & moyenne_cover_par_indiv<=rate_max_cover_ind))
# select final row (ind) to keep 
vcf_filtered_ind_posi<-vcf_filtered_posi[,c("FORMAT",noms)]

# get position DP mean (row)
dp<- extract.gt(vcf_filtered_ind_posi, element='DP', as.numeric = TRUE) 
moyenne_cover_par_position = apply(dp,1,mean,na.rm=T)
# filter individuals according to the rate specified in arg
vcf_filtered_posi<-subset(vcf, moyenne_cover_par_position>=rate_min_cover_pos & moyenne_cover_par_position<=rate_max_cover_pos)

return(vcf_filtered_posi)
}
