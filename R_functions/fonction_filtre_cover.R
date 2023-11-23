#########################
# fonction filtre cover #
#########################

#SEUIL SUR LA COUVERTURE (bas et haut): enlever les SNPs qui depassent les seuils. 
fun_filtre_cover<-function(vcf,taux_min_cover_pos,taux_max_cover_pos,taux_min_cover_ind,taux_max_cover_ind){

# Filters DP on position
dp<- extract.gt(vcf, element='DP', as.numeric = TRUE) 
moyenne_cover_par_position = apply(dp,1,mean,na.rm=T)
vcf_filtered_posi<-subset(vcf, moyenne_cover_par_position>=taux_min_cover_pos & moyenne_cover_par_position<=taux_max_cover_pos)

# Filters DP on individu
dpi<- extract.gt(vcf_filtered_posi, element='DP', as.numeric = TRUE) #récupère la cover 
moyenne_cover_par_indiv = apply(dpi,2,mean,na.rm=T)
noms<-names(which(moyenne_cover_par_indiv>=taux_min_cover_ind & moyenne_cover_par_indiv<=taux_max_cover_ind))
vcf_filtered_ind_posi<-vcf_filtered_posi[,c("FORMAT",noms)]

return(vcf_filtered_ind_posi)
}


# dpi<- extract.gt(vcf_initial, element='DP', as.numeric = TRUE) #récupère la cover 
# moyenne_cover_par_indiv = apply(dpi,2,mean,na.rm=T)
# noms<-names(which(moyenne_cover_par_indiv>20))
# vcf_filtered_ind<-vcf_initial[,c("FORMAT",noms)]
# View(vcf_filtered_ind@gt)
