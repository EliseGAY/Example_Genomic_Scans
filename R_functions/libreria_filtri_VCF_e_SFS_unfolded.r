#######################################
# Functions to run genomic scan
# Author : Stefano Mona (EPHE - MNHN)
#######################################

#library(vcfR)
#library(rlist)
#library(PopGenome)
#library("devtools")
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install_github("whitlock/OutFLANK")
#library("OutFLANK")
#library(pegas)
#library(ggplot2)
#library(adegenet)
#install_github("jgx65/hierfstat")
#library(hierfstat)




######### funzioni contenute in questa libreria
##### 1) elimina missing data (vuole come input solo un VCF)
##### 2) elimina siti fissati 0/0 (vuole come input solo un VCF)
##### 3) elimina siti fissati 1/1 (vuole come input solo un VCF)
##### 4) elimina siti che cadono nella parte terminale del frammento (vuole come input VCF e la posizione da cui poi tagliare, per es. 115)
##### 5) elimina siti troppo heterozigoti (vuole come input VCF e la soglia di het, per esempio 0.8)
##### 5a) elimina loci con siti troppo heterozigoti (vuole come input VCF e la soglia di het, per esempio 0.8)
##### 6) elimna loci con troppi SNP (vuole come input VCF e il max numero di siti permessi per locus, per esempio 3)
##### 7) calcola SFS unfolded (vuole come input VCF e il nome dell'outgroup. non ci devono essere dati mancanti)
##### 7a) calcola SFS unfolded (vuole come input VCF e il nome dell'outgroup. non ci devono essere dati mancanti): corretta per il fatto che l'outgroup non abbia eterozigoti
##### 8) elimina siti che non son presenti in tutte le pop fornite (vuole come input VCF e una lista con i nomi degli individui per popolazione)
##### 9) calola Fst globale e pairwise e fa bootstrap per entrambi. inoltre restituisce distribuzione Fst per sito. Ha bisogno della library rlist e popgenome: il vcf deve essere letto con popgenome
##### 9a) calola Fst globale e pairwise e fa bootstrap per entrambi modificata per triplodi. inoltre restituisce distribuzione Fst per sito. Ha bisogno della library rlist e popgenome: il vcf deve essere letto con popgenome
##### 10) calola outflanck usando le librerie vcfR e outflanck. restituisce lista con oggetto outflanck su cui fare il plot e anche fst, pvalues etc. se ci sono outliers, a cercare negli output!
##### 11) scrive file arp a partire da un percorso verso un vcf e una lista con le pop. vuole libreria vcfR. scrive output nella cartella del vcf col suo nome con estensione cambiata in arp. 
##### 12) elimino SNP non in HW. Ho bisogno di vcfR e pegas per il test. Come input vuole un vcfR e una lista di individui (una pop). Se non specifico niente, considera tutti gli ind presenti come una pop. Restituisce vcf pulito e il valore del test per gli SNP eliminati
##### 12a) elimino SNP non in HW. Ho bisogno di vcfR e pegas per il test, le funzioni 2,3 e 12. Come input vuole un percorso per un vcf e una lista di pop). Restituisce vcf pulito da SNP non HW in ogni pop e il num di SNP eliminati per pop.
##### 13) calcolo bootstrap Fis sulle pop. Input: vcfR e lista di pop e fisso bootstrap=1000. Ha bisogno di hierfstat, pegas, vcfR. Restituisce matrice dei CI delle Fis per pop (definite nella lista pop)
##### 14) Input: un vcfR e MAF (par default a 0.05). Ha bisogno solo di vcfR. Restituisce una lista con: 1) vettore frequenze per SNP; 2) vcfR pulito dagli SNPs con frequenza minore di MAF.
##### 15) Input: semplice vettore di SFS folded. Mi restituisce matrice per fare il plot dell'SFS nomalizzato (come in Lapierre et al. Genetics, 2017).
##### 16) crea reads a partire da un file arlequin. ha bisogno di rlist. Solo argomento fisso, percorso file arp. Il resto sono settati par default (lungheza read, coverage, errore  illumina). L'ultimo argomento é la versione di fastsimcoal. Di default é la fsc226, fare attenzione!! IL coverage é PER ALLELE!!!!
##### 17) Calcola theta pi a partire d'un SFS folded o unfolded. biogna solo dirgli se folded=T/F
##### 17a) come la 17 ma senza singleton.
##### 17b) come la 17 ma senza doubleton.
##### 18) elimina_SNP_cov:togli SNP con troppo o troppo poco cov. vuole input vcf, min cov, max cov. restituisce vcf pulito
##### 19) downsampling:downsample un SFS folded di X classi in X-n classi (x avere uno spettro con meno individui). come argomento vuole lo spetto di partenza e il nuovo numero di individui
##### 19a) downsample_numero_SNP:downsample un SFS folded ad un numero ridotto di SNPs. come argomento vuole lo spetto di partenza e il nuovo numero di SNPs
##### 20) tajima D: ha bisogno solo di un folded SFS. Per calcolare theta_P perà usa la funzione 17!!!!!
##### 21) fst_reynolds_sfs2D: ha bisogno solo diuno spettro 2D e un valore di maf. Restituisce un vettore con Fst di Reynolds et al. 1983(unweigheted, weighted)
##### 22) 22a) GENO: transforma un vcf in file .geno (0,1,2,9 o NaN se serve per 2D-SFS e Fst). Funziona per riga (usare apply sulla matrice dei genotipi)
##### 23) Conta allele alternativo: ha bisogno di matrice .geno e una lista pop (funziona per riga, usare apply). Conta l'allele ALT in ogni pop e mette il totale in ultima colonna
##### 24) calcola 2S-SFS pairwise per fastsimcoal. Ha bisogno di funzione 23 e 24. vuole un vcf e una lista di pop. Restitiuisce lista di matrici pairwise
##### 25) calcola fst di reynolds pairwise. Ha busogno di matrice prodotta da FUNZIONE 23 come input. Restituisce vettore di Fst
##### 26) calcola bootstrap e fst di reynolds pairwise. Ha busogno FUNZIONE 25 etc. Restituisce lista con elemnto Fst osservate e elemento Fst bootstrap
##### 27) 27)a sliding windows di Fst pairwise e SFS. La versione a) fa solo SFS, la uso quando ho una sola pop.
##### 28) Vuole un vcfR e un vettore di pos totali chiamate. Parametri: lunghezza locus e distanza dal prossimo locus. Risultato: vcfR dei loci selezionati + numero tot basi chiamate.
##### 29 a) Funzione per ANtoine, calcola Num di SNP all het in pop1 et all fixed (either 0/0 or 1/1) in pop2
##### 29 b) Funzione per ANtoine, calcola Num di SNP all het in pop1 et all fixed (only 0/0) in pop2
##### 29 c) Funzione per Elise, calcola Num di SNP fixed between pop1 and pop2
##### 30 a) Funzione per ANtoine, sliding window PCA per SNP. Takes vcf, returns positions (up and low) of the window and the % of variabce of the 2 first axis
##### 30 b) Funzione per ANtoine, sliding window PCA per base pairs. Takes vcf, returns positions (up and low) of the window and the % of variabce of the 2 first axis
##### 31)   FST by SNP

################# FUNZIONE 19  #############################
########################### downsample un SFS folded di X classi in X-n classi. 

downsampling<-function(spettro_or, new_size){
l_spettro_or<-length(spettro_or)
risultato<-c()
for (i in 1:l_spettro_or) {
for (j in 1:spettro_or[i]) {
sequenza<-c(rep(1,i),rep(0,((2*l_spettro_or)-i)))
resampling<-sample(sequenza,2*new_size)
risultato<-c(risultato, (sum(resampling)))
}
}

finale<-table(risultato)
finale<-finale[-1]

##### controllo che ci sian tutte le classi, in caso contrario aggiungo una classe con valore 0.
##### mi serve perché altrimenti potrei perdere il buon numero di bin!!!!!
aaa=as.numeric(row.names(finale))
new_ris<-c()
a<-0
for (i in 1:new_size){

if (is.element(i, aaa)) {
new_ris<-c(new_ris,finale[i-a])
}
else {
new_ris<-c(new_ris,0)
a=a+1
}
}
#### dentro new_ris ho quindi il vettore con tutte le classi piu eventualmente quelle che sono passati da MAF a aver freq > 50%
#### ora devo aggiungere questi SNP alla loro classe corretta, secondo sempre la loro frequnza
#### che é in aaa. questi SNP son in finale ma non in new_size
for (i in 1: length(aaa)){
if (aaa[i]>new_size) {
new_ris[(2*new_size)-aaa[i]]<-new_ris[(2*new_size)-aaa[i]]+finale[aaa[i]]
}
}
new_ris<-as.numeric(new_ris)
##### adesso ho i valori corretti i new_size
return(new_ris)
}
#################################################

################# FUNZIONE 19a  #############################
########################### downsample un SFS folded per numero di SNP. 

downsample_numero_SNP<-function(sfs_originale, n_snps){ #### sfs_originale is the sfs to downsample; n_snp the number of snps we want in the new SFS
per_media<-c()
for (j in 1:1000){

totale<-c()
for (i in 1:length(sfs_originale)){
totale<-c(totale, rep(i,sfs_originale[i]))
}

tot_resampling<-sample(totale,n_snps)
qqq<-table(tot_resampling)

aaa=as.numeric(row.names(qqq))
new_ris<-c()
a<-0
for (i in 1:length(sfs_originale)){

if (is.element(i, aaa)) {
new_ris<-c(new_ris,qqq[i-a])
}
else {
new_ris<-c(new_ris,0)
a=a+1
}
}
per_media<-cbind(per_media,as.vector(new_ris))

}

finale_res_sfs<-apply(per_media,1, mean)

return(round(finale_res_sfs))

}
#################################################





################# FUNZIONE 18  #############################
########################### TOGLI SNPS CON COV < minco E > maxco
elimina_SNP_cov<-function(dati, minco, maxco){

coverage<-extract.gt(dati, element = "DP", mask = FALSE, as.numeric=T) ### prendo il coverage
coverage[coverage == 0] <- NA
cov_per_snp<-apply(coverage, 1, mean,na.rm=TRUE)

elimin_cov<-c()
for (i in 1: length(cov_per_snp)) {
if (cov_per_snp[i]<minco | cov_per_snp[i]>maxco ) {elimin_cov<-c(elimin_cov,i)}
}

if (length(elimin_cov)>0) {
dati_clean<-dati[-elimin_cov,]
}
else {
dati_clean<-dati
}

return(dati_clean)
}

#################################################



################# FUNZIONE 1  #############################
##############################################################################################################################
elimina_missing_data<-function(dati){
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
vettore_snps<-getID(dati)
sss=(strsplit(vettore_snps, "_"))
aaa=matrix(unlist(sss),ncol=2,byrow=TRUE)
eee=as.numeric(aaa[,2]) ##### vettore contenente le posizioni di ogni SNP sul suo frammento

snps_tot<-length(eee)
elimina_snps<-c()  #### 

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

for (i in 1:snps_tot) {
temp<-0
for (j in 1:ncol(genotype)) {
if (genotype[i,j]=="./.") {temp<-temp+1}
}
if (temp > 0) {elimina_snps<-rbind(elimina_snps,i)}
}


if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}

return(dati_clean)
}
#############################################################################################################################



################# FUNZIONE 1a per STACKS 2.5  #############################
##############################################################################################################################
elimina_missing_data_stacks25<-function(dati){

elimina_snps<-c()  #### 

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

for (i in 1:nrow(genotype)) {
temp<-0
for (j in 1:ncol(genotype)) {
if (is.na(genotype[i,j])) {temp<-temp+1}
}
if (temp > 0) {elimina_snps<-rbind(elimina_snps,i)}
}


if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}

return(dati_clean)
}
#############################################################################################################################




################# FUNZIONE 2  #############################
################################# tolgo SNP fissati per 0/0 poter calcolare il numero di loci variabili ################
elimina_snp_fixed_0<-function(dati){ 

vettore_snps<-getID(dati)

snps_tot<-length(vettore_snps) #### numero snp
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
n_ind<-length(vettore_nomi) ### numero individui


elimina_snps_0<-c()  ####

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

for (i in 1:snps_tot) {
temp<-0
for (j in 1:ncol(genotype)) {
if (genotype[i,j]=="0/0") {temp<-temp+1}
}
if (temp == n_ind) {elimina_snps_0<-rbind(elimina_snps_0,i)}
}

if (length(elimina_snps_0)>0) {
dati_clean<-dati[-elimina_snps_0,]
}
else {
dati_clean<-dati
}

return(dati_clean)
}
######################################################################################################################

################# FUNZIONE 3  #############################
################################# tolgo SNP fissati per 1/1 poter calcolare il numero di loci variabili ################
elimina_snp_fixed_1<-function(dati){ 

vettore_snps<-getID(dati)

snps_tot<-length(vettore_snps) #### numero snp
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
n_ind<-length(vettore_nomi) ### numero individui


elimina_snps_1<-c()  ####

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

for (i in 1:snps_tot) {
temp<-0
for (j in 1:ncol(genotype)) {
if (genotype[i,j]=="1/1") {temp<-temp+1}
}
if (temp == n_ind) {elimina_snps_1<-rbind(elimina_snps_1,i)}
}

if (length(elimina_snps_1)>0) {
dati_clean<-dati[-elimina_snps_1,]
}
else {
dati_clean<-dati
}



return(dati_clean)
}
########################################################################

################# FUNZIONE 4  #############################
################################# FUNZIONE CHE TOGLIE SNPs che cadono nella parte finale del frammento
elimin_snps_fine_frammento<-function (dati, val_frammento) { ### dati é un oggetto vcf, val_frammento é la posizione del frammento dopo la quale elminiare tutti gli snps
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
vettore_snps<-getID(dati)
sss=(strsplit(vettore_snps, "_"))
aaa=matrix(unlist(sss),ncol=2,byrow=TRUE)
eee=as.numeric(aaa[,2]) ##### vettore contenente le posizioni di ogni SNP sul suo frammento
snps_tot<-length(eee)
elimina_snps<-c()  #### con questo ciclo creo un vettore dove metto tutte le posizioni corrispondenti a SNP che cade dall'87 alla fine del frammento
for (i in 1:snps_tot) {
if (eee[i] >= val_frammento) {elimina_snps<-rbind(elimina_snps,i)}
}
if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}
return(dati_clean)
}
#######################################################################




################# FUNZIONE 4a STACKS 2.5  #############################
################################# FUNZIONE CHE TOGLIE SNPs che cadono nella parte finale del frammento
elimin_snps_fine_frammento_stacks25<-function (dati, val_frammento) { ### dati é un oggetto vcf, val_frammento é la posizione del frammento dopo la quale elminiare tutti gli snps
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
vettore_snps<-getCHROM(dati)
vettore_snps_pos<-getPOS(dati)
snps_tot<-length(vettore_snps_pos)

elimina_snps<-c()  #### con questo ciclo creo un vettore dove metto tutte le posizioni corrispondenti a SNP che cade dall'87 alla fine del frammento
for (i in 1:snps_tot) {
if (vettore_snps_pos[i] >= val_frammento) {elimina_snps<-rbind(elimina_snps,i)}
}

if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}
return(dati_clean)
}
#######################################################################


################# FUNZIONE 5  #############################
################################# FUNZIONE togli  SNPs sono eterogizoti al val_hetero
togli_hetero<-function (dati, val_hetero) {

genotype_clean<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

lista_het<-c()

for (j in 1:nrow(genotype_clean)) {
temp<-0
for (i in 1:ncol(genotype_clean)) {
if (genotype_clean[j,i]=="0/1") {temp<-temp+1}
}
if (temp/(ncol(genotype_clean)) >= val_hetero) {lista_het<-rbind(lista_het,j)}
}

if (length(lista_het)>0) {
dati_clean_het<-dati[-lista_het,]
}
else {
dati_clean_het<-dati
}


return(dati_clean_het)
}

#######################################################################

#########################################################################################################
################################# FUNZIONE 5A togli loci con  SNPs sono eterogizoti al val_hetero
togli_loci_hetero<-function (dati, val_hetero) {

genotype_clean<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

lista_het<-c()

for (j in 1:nrow(genotype_clean)) {
temp<-0
for (i in 1:ncol(genotype_clean)) {
if (genotype_clean[j,i]=="0/1") {temp<-temp+1}
}
if (temp/(ncol(genotype_clean)) >= val_hetero) {lista_het<-rbind(lista_het,j)}
}


vettore_snps<-getID(dati)
sss=(strsplit(vettore_snps, "_"))
aaa=matrix(unlist(sss),ncol=2,byrow=TRUE)
eee=as.numeric(aaa[,2]) ##### vettore contenente le posizioni di ogni SNP sul suo frammento
n_loci<-length(table(as.numeric(aaa[,1])))   ###numero loci
n_snp<-length(vettore_snps)
snps_frag<-as.numeric(aaa[,1])### c


lista_loci<-snps_frag[lista_het]

if (length(lista_loci)>0) {

index_lista_loci<-c()
for(i in 1:length(lista_loci))
{
for(j in 1:length(snps_frag))
{
if(lista_loci[i]==snps_frag[j])
index_lista_loci<-c(index_lista_loci,j)
}
}
dati_clean<-dati[-index_lista_loci,]
} else {
dati_clean<-dati
}



return(dati_clean)
}

################################################################################################


################# FUNZIONE 6  #############################
################################ FUNZIONE elimino loci con troppi SNPs usando max_snp#######

togli_loci_troppi_SNP<-function(dati, max_snp) {


vettore_snps<-getID(dati)
sss=(strsplit(vettore_snps, "_"))
aaa=matrix(unlist(sss),ncol=2,byrow=TRUE)
eee=as.numeric(aaa[,2]) ##### vettore contenente le posizioni di ogni SNP sul suo frammento
n_loci<-length(table(as.numeric(aaa[,1])))   ###numero loci
n_snp<-length(vettore_snps)
snps_frag<-as.numeric(aaa[,1])### c

per_elim_loc<-table(snps_frag)
lista_loci<-c()

for (i in 1:dim(per_elim_loc)) {
if (as.matrix(per_elim_loc[i])>max_snp) {
aaa<-rownames(as.matrix(per_elim_loc[i]))
lista_loci<-cbind(lista_loci,aaa)
}
}
lista_loci<-as.numeric(lista_loci) #### contiene i nomi dei loci da elminare

if (length(lista_loci)>0) {

index_lista_loci<-c()
for(i in 1:length(lista_loci))
{
for(j in 1:length(snps_frag))
{
if(lista_loci[i]==snps_frag[j])
index_lista_loci<-c(index_lista_loci,j)
}
}
dati_clean<-dati[-index_lista_loci,]
} else {
dati_clean<-dati
}


return(dati_clean)
}
################################################################################################

################# FUNZIONE 6a STACKS 2.5  #############################
################################ FUNZIONE elimino loci con troppi SNPs usando max_snp#######

togli_loci_troppi_SNP_stacks25<-function(dati, max_snp) {

vettore_snps<-getCHROM(dati)
vettore_snps_pos<-getPOS(dati)

n_loci<-length(table(as.numeric(vettore_snps)))   ###numero loci
n_snp<-length(vettore_snps)

snps_frag<-as.numeric(vettore_snps)###

per_elim_loc<-table(snps_frag)
lista_loci<-c()

for (i in 1:dim(per_elim_loc)) {
if (as.matrix(per_elim_loc[i])>max_snp) {
aaa<-rownames(as.matrix(per_elim_loc[i]))
lista_loci<-cbind(lista_loci,aaa)
}
}
lista_loci<-as.numeric(lista_loci) #### contiene i nomi dei loci da elminare

if (length(lista_loci)>0) {

index_lista_loci<-c()
for(i in 1:length(lista_loci))
{
for(j in 1:length(snps_frag))
{
if(lista_loci[i]==snps_frag[j])
index_lista_loci<-c(index_lista_loci,j)
}
}
dati_clean<-dati[-index_lista_loci,]
} else {
dati_clean<-dati
}
return(dati_clean)
}
################################################################################################



################# FUNZIONE 7  #############################
################################ funziona senza MISSING DATA #########################
MySFS_unfolded<-function(x,outgroup_nom){
x_gt<-extract.gt(x,element = "GT",convertNA = F,return.alleles = F)
x_name<-rownames(x_gt)
#x_name<-str_split_fixed(x_name,pattern ="_",n = 2) #split le nom des locus
#x_name<- x_name[,1]
#print(length(x_name))
#print("Nombre de SNPs avant")

#print(length(unique(x_name)))
#print("Nombre de Loci avant")



########### on enlève les hétérozbgote
outgroup_nom<-as.character(outgroup_nom)
outgroup<-which(colnames(x_gt)==outgroup_nom)
print(outgroup)
print("index_outgroup")
a<-x_gt
torm<-NULL
for(i in 1:nrow(a)){
  if(a[i,outgroup]== "0/1")
    torm<-rbind(torm,i)}


print(length(torm))
print("nombre de SNPs heterozygote chez l'outgroup")

x_name_after <-x_name[-torm]
print(length(x_name_after))
print("Nombre de SNPs apres")

print(length(unique(x_name_after)))
print("Nombre de Loci apres")
a1<-a[-torm,]
a<-a[-torm,-outgroup]
outgroup_1<-outgroup+1
clean_vcf<-x[-torm,-outgroup_1]

SNPs_freq=NULL
for(i in 1:nrow(a)){
  n=0
  n1=ncol(a)*2
  for(b in 1:ncol(a)){
    if(a[i,b]==a1[i,outgroup]){
      n=n+2}
      if(a[i,b]=="0/1")
        n=n+1}
  SNPs_freq<-rbind(SNPs_freq,(n1-n))
}

rownames(SNPs_freq)<- x_name_after


akl<-table(SNPs_freq,exclude = NULL)
print(sum(akl))
mytable=NULL
for (i in 0:(ncol(a)*2)){
  n=0
  for (y in 1:nrow(SNPs_freq)){
    if(SNPs_freq[y,]==i){
      n=n+1}}
  mytable<-rbind(mytable,n)}
 

print("verification de la valeur de la SFS")
return(list(as.matrix(mytable),clean_vcf))
}

##### la prima e l'ultima classe sono i siti fissati, vanno tolti###############
################################################################################################################


################# FUNZIONE 7a  #############################
################################ funziona senza MISSING DATA #########################
MySFS_unfolded_7a<-function(x,outgroup_nom){
x_gt<-extract.gt(x,element = "GT",convertNA = F,return.alleles = F)
x_name<-rownames(x_gt)
#x_name<-str_split_fixed(x_name,pattern ="_",n = 2) #split le nom des locus
#x_name<- x_name[,1]
#print(length(x_name))
#print("Nombre de SNPs avant")

#print(length(unique(x_name)))
#print("Nombre de Loci avant")



########### on enlève les hétérozbgote
outgroup_nom<-as.character(outgroup_nom)
outgroup<-which(colnames(x_gt)==outgroup_nom)
print(outgroup)
print("index_outgroup")
a<-x_gt
torm<-NULL
for(i in 1:nrow(a)){
  if(a[i,outgroup]== "0/1")
    torm<-rbind(torm,i)}


print(length(torm))
print("nombre de SNPs heterozygote chez l'outgroup")

if (length(torm)>0) {
x_name_after <-x_name[-torm]
}
else {
x_name_after <-x_name
}

print(length(x_name_after))
print("Nombre de SNPs apres")

print(length(unique(x_name_after)))
print("Nombre de Loci apres")
if (length(torm)>0) {
a1<-a[-torm,]
a<-a[-torm,-outgroup]
outgroup_1<-outgroup+1
clean_vcf<-x[-torm,-outgroup_1]
}
else {
a1<-a
a<-a[,-outgroup]
outgroup_1<-outgroup+1
clean_vcf<-x[,-outgroup_1]
}




SNPs_freq=NULL
for(i in 1:nrow(a)){
  n=0
  n1=ncol(a)*2
  for(b in 1:ncol(a)){
    if(a[i,b]==a1[i,outgroup]){
      n=n+2}
      if(a[i,b]=="0/1")
        n=n+1}
  SNPs_freq<-rbind(SNPs_freq,(n1-n))
}

rownames(SNPs_freq)<- x_name_after


akl<-table(SNPs_freq,exclude = NULL)
print(sum(akl))
mytable=NULL
for (i in 0:(ncol(a)*2)){
  n=0
  for (y in 1:nrow(SNPs_freq)){
    if(SNPs_freq[y,]==i){
      n=n+1}}
  mytable<-rbind(mytable,n)}
 

print("verification de la valeur de la SFS")
return(list(as.matrix(mytable),clean_vcf))
}
##### la prima e l'ultima classe sono i siti fissati, vanno tolti###############
################################################################################################################


################# FUNZIONE 8  #############################
#################################################################################
elimina_snp_non_condivisi<-function(dati, lista_pop) {  ###2 parametri: 1 vcf e una lista contenente i vettori dei nomi di ciascun ind per pop. restituisce un vcf con snp che siano in tutte le pop della lista

format<-c("FORMAT")

dati1<-dati[,c(format,unlist(lista_pop))]
n_pop<-length(lista_pop)


genotype<-extract.gt(dati1, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
snps_tot<-nrow(genotype)

elimina_snps<-c()  #### 
for (i in 1:snps_tot) {

pop_missing<-0

for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'é in almeno un individuo
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
#### fine ciclo
temp<-0
for (j in born_inf:born_sup) {
if (genotype[i,j]=="./.") {temp<-temp+1}
}
if (temp == (born_sup-born_inf+1)) {pop_missing<-pop_missing+1}

}
if (pop_missing > 0) {elimina_snps<-rbind(elimina_snps,i)}
}

if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}


return(dati_clean)
}
#### resistituisce il nuovo vcf pulito dagli snps che non son presenti in tutte le popolazioni
##################################################################################################################



################# FUNZIONE 9  #############################
#################################################################################
boot_popgenome<-function(dati, lista_pop, bootstrap) {  ### dati é un vcf letto con popgenome!

snps_tot<-dati@n.sites

all<-unlist(lista_pop)
n_pop<-length(lista_pop)

##### bootstrap sulla global Fst#######################################
ris_boot_fst_all<-c()
for (i in 1:bootstrap){

qq<-sample(all)  #### faccio uno shaffle degli individui per allocarli random alle popolazioni
new_list<-list()
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-as.vector(qq[born_inf:born_sup])
new_list<-list.append(new_list,ss)
}
### fine ciclo


dati_PG2<-set.populations(dati,new_list)
ee<-F_ST.stats(dati_PG2)
ris_boot_fst_all<-cbind(ris_boot_fst_all,ee@nucleotide.F_ST)
#ris_boot_fst_pairwise<-cbind(ris_boot_fst_pairwise,ee@nuc.F_ST.pairwise)
}
###### fine bootstrap che mi ritorna vettore ris_boot_fst_all

### calcolo valori nei veri dati
ris_veri<-set.populations(dati,lista_pop)
mat_vera<-F_ST.stats(ris_veri)
####

#### scrivo file finali coi veri valori e i CI per bootstrap
ris_FST_all<-c(as.vector(mat_vera@nucleotide.F_ST),as.vector(quantile(ris_boot_fst_all, prob=c(0.025, 0.5, 0.975))),format(snps_tot, scientific=F))
names(ris_FST_all)<-c("Fst","2.5%","50%","97.5%","nloci/snps")

####################################### bootstrap sulle comparazioni pairwise
finale_pairwise<-c()
for (k in 1:(n_pop-1)){
for (h in (k+1):n_pop){

lista_pairwise<-lista_pop[c(k,h)[]] #### creo la lista di due popolazioni

all_pairwise<-unlist(lista_pairwise)
n_pop_new<-length(lista_pairwise)

##### bootstrap Fst#######################################
ris_boot_fst_pairwise<-c()
for (i in 1:bootstrap){

qq<-sample(all_pairwise)  #### faccio uno shaffle degli individui per allocarli random alle popolazioni
new_list<-list()
for (z in 1:n_pop_new){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pairwise[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pairwise[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pairwise[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-as.vector(qq[born_inf:born_sup])
new_list<-list.append(new_list,ss)
}
### fine ciclo


dati_PG2<-set.populations(dati,new_list)
ee<-F_ST.stats(dati_PG2)
ris_boot_fst_pairwise<-cbind(ris_boot_fst_pairwise,ee@nucleotide.F_ST)
}
######
ris_veri_pairwise<-set.populations(dati,lista_pairwise)
mat_vera_pairwise<-F_ST.stats(ris_veri_pairwise)
####
ris_FST_pairwise<-c(as.vector(mat_vera_pairwise@nucleotide.F_ST),as.vector(quantile(ris_boot_fst_pairwise, prob=c(0.025, 0.5, 0.975))))
finale_pairwise<-rbind(finale_pairwise,ris_FST_pairwise)
}
}
################################################################################

vett_nomi_righe<-c()
for (k in 1:(n_pop-1)){
for (h in (k+1):n_pop){
aaa<-paste("pop_",k,"/","pop_",h, sep="")
vett_nomi_righe<-rbind(vett_nomi_righe,aaa)
}
}
colnames(finale_pairwise)<-c("Fst","2.5%","50%","97.5%")
rownames(finale_pairwise)<-vett_nomi_righe




aaa<-detail.stats(ris_veri,site.FST=TRUE)
dist_snp<-aaa@region.stats@site.FST[[1]] ### contiene il vettore di fst per SNP

finale<-list(ris_FST_all,finale_pairwise,dist_snp)
return(finale)

}
##################################################################################################################


################# FUNZIONE 9a modficata triploidi  #############################
#################################################################################
boot_popgenome_trip<-function(dati, lista_pop, bootstrap) {  ### dati é un vcf letto con popgenome!

snps_tot<-dati@n.sites2

all<-unlist(lista_pop)
n_pop<-length(lista_pop)

##### bootstrap sulla global Fst#######################################
ris_boot_fst_all<-c()
for (i in 1:bootstrap){

qq<-sample(all)  #### faccio uno shaffle degli individui per allocarli random alle popolazioni
new_list<-list()
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-as.vector(qq[born_inf:born_sup])
new_list<-list.append(new_list,ss)
}
### fine ciclo


dati_PG2<-set.populations(dati,new_list, diploid=FALSE,triploid=TRUE)
ee<-F_ST.stats(dati_PG2)
ris_boot_fst_all<-cbind(ris_boot_fst_all,ee@nucleotide.F_ST)
#ris_boot_fst_pairwise<-cbind(ris_boot_fst_pairwise,ee@nuc.F_ST.pairwise)
}
###### fine bootstrap che mi ritorna vettore ris_boot_fst_all

### calcolo valori nei veri dati
ris_veri<-set.populations(dati,lista_pop, diploid=FALSE,triploid=TRUE)
mat_vera<-F_ST.stats(ris_veri)
####

#### scrivo file finali coi veri valori e i CI per bootstrap
ris_FST_all<-c(as.vector(mat_vera@nucleotide.F_ST),as.vector(quantile(ris_boot_fst_all, prob=c(0.025, 0.5, 0.975))),format(snps_tot, scientific=F))
names(ris_FST_all)<-c("Fst","2.5%","50%","97.5%","numloci/snps")

####################################### bootstrap sulle comparazioni pairwise
finale_pairwise<-c()
for (k in 1:(n_pop-1)){
for (h in (k+1):n_pop){

lista_pairwise<-lista_pop[c(k,h)[]] #### creo la lista di due popolazioni

all_pairwise<-unlist(lista_pairwise)
n_pop_new<-length(lista_pairwise)

##### bootstrap Fst#######################################
ris_boot_fst_pairwise<-c()
for (i in 1:bootstrap){

qq<-sample(all_pairwise)  #### faccio uno shaffle degli individui per allocarli random alle popolazioni
new_list<-list()
for (z in 1:n_pop_new){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pairwise[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pairwise[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pairwise[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-as.vector(qq[born_inf:born_sup])
new_list<-list.append(new_list,ss)
}
### fine ciclo


dati_PG2<-set.populations(dati,new_list, diploid=FALSE,triploid=TRUE)
ee<-F_ST.stats(dati_PG2)
ris_boot_fst_pairwise<-cbind(ris_boot_fst_pairwise,ee@nucleotide.F_ST)
}
######
ris_veri_pairwise<-set.populations(dati,lista_pairwise, diploid=FALSE,triploid=TRUE)
mat_vera_pairwise<-F_ST.stats(ris_veri_pairwise)
####
ris_FST_pairwise<-c(as.vector(mat_vera_pairwise@nucleotide.F_ST),as.vector(quantile(ris_boot_fst_pairwise, prob=c(0.025, 0.5, 0.975))))
finale_pairwise<-rbind(finale_pairwise,ris_FST_pairwise)
}
}
################################################################################

vett_nomi_righe<-c()
for (k in 1:(n_pop-1)){
for (h in (k+1):n_pop){
aaa<-paste("pop_",k,"/","pop_",h, sep="")
vett_nomi_righe<-rbind(vett_nomi_righe,aaa)
}
}
colnames(finale_pairwise)<-c("Fst","2.5%","50%","97.5%")
rownames(finale_pairwise)<-vett_nomi_righe




aaa<-detail.stats(ris_veri,site.FST=TRUE)
dist_snp<-aaa@region.stats@site.FST[[1]] ### contiene il vettore di fst per SNP

finale<-list(ris_FST_all,finale_pairwise,dist_snp)
return(finale)

}
##################################################################################################################


################# FUNZIONE 10  #############################
#################################################################################
calcola_outflanck<-function(dati, lista_pop) {  ###


format<-c("FORMAT")

dati1<-dati[,c(format,unlist(lista_pop))]

###### crea input per outflanck ###################
genotype_clean<-extract.gt(dati1, element = "GT", mask = FALSE, as.numeric = F ,return.alleles = F, IDtoRowNames = TRUE, extract = TRUE,convertNA = TRUE) ### prendo il coverage per fare heatmap
for (i in 1:ncol(genotype_clean)) #coding all 0/0 as 0, 0/1 as 1 and 1/1 as 2 (0=ref, 1=alt)
{
genotype_clean[genotype_clean[,i]=="0/0",i]<-0 #homo anc
genotype_clean[genotype_clean[,i]=="0/1",i]<-1 #het
genotype_clean[genotype_clean[,i]=="1/1",i]<-2#homo alt
genotype_clean[genotype_clean[,i]=="./.",i]<-9# aggiunto da me perché ho i missing data
genotype_clean[genotype_clean[,i]=="./0",i]<-9# aggiunto da me perché ho i missing data
genotype_clean[genotype_clean[,i]=="./1",i]<-9# aggiunto da me perché ho i missing data
print (i)
}
genotype_clean[is.na(genotype_clean)] <- 9
SNPmat<-t(genotype_clean) #### matrice dati che vuole OUTFLANCK

vettore_snps<-getPOS(dati1)
locusNames<-vettore_snps


n_pop<-length(lista_pop)
vett_pop<-c()
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-rep(z,(born_sup-born_inf+1))
vett_pop<-c(vett_pop,ss)
vett_pop<-as.factor(vett_pop)
}
### fine ciclo che mi da il vettore di fattori per identificare la pop d'appartenenza di ciascun individuo


FstDataFrame <- MakeDiploidFSTMat(SNPmat, locusNames, vett_pop)



plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, xlim = c(-0.01,1.0), ylim = c(-0.01, 1.0), pch = 20)
hist(FstDataFrame$FSTNoCorr) 


OF <- OutFLANK(FstDataFrame, NumberOfSamples=n_pop, qthreshold = 0.05, RightTrimFraction = 0.1, LeftTrimFraction = 0.1, Hmin=0.1)
outliers_OF <- OF$results$LocusName[OF$results$OutlierFlag == TRUE]
print(outliers_OF)
length(outliers_OF)

### grafico con Het et Fst + outliers
punti_out <- pOutlierFinderChiSqNoCorr(FstDataFrame, Fstbar = OF$FSTNoCorrbar, dfInferred = OF$dfInferred, qthreshold = 0.05, Hmin=0.1)
my_out <- punti_out[punti_out$OutlierFlag==TRUE,]  ### qui metto gli outliers, nel caso ce ne siano
#plot(punti_out$He, punti_out$FST, pch=19, col=rgb(0,0,0,0.1))
#points(my_out$He, my_out$FST, col="blue")
##

### grafico p-values
hist(punti_out$pvaluesRightTail) ### grafico p-values
###

### Manhattan plot + outliers
#plot(punti_out$LocusName[punti_out$He>0.1], punti_out$FST[punti_out$He>0.1],xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
#points(my_out$LocusName, my_out$FST, col="magenta", pch=20) 
###

### output classico di OutFlanck
###OutFLANKResultsPlotter(OF, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
########

lista_finale<-list(OF, punti_out, my_out)

return(lista_finale)

}
########################################################################################################################################################



################# FUNZIONE 11  #############################
#################################################################################
scrivi_arp<-function(percorso, lista_pop) {  ###


nome<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(percorso))
path<-dirname(percorso)
nome_output<-paste(path,"/", nome, ".arp", sep="")



data<-read.vcfR(percorso, verbose = TRUE)

#### ordino i dati secondo la lista pop ####
format<-c("FORMAT")
dati<-data[,c(format,unlist(lista_pop))]
###########################################

####### calcolo la matrice dati per alrequin##################
mat_temp<-extract.haps(dati,unphased_as_NA = FALSE)
mat_temp<-t(mat_temp)
### trasformo NA in ? per arlequin
mat_temp[is.na(mat_temp)] <- 0
for (i in 1:nrow(mat_temp)) {
for (j in 1:ncol(mat_temp)) {
if (mat_temp[i,j] == "0") {mat_temp[i,j]="?"}
}
}

nomi_ind<-dimnames(mat_temp)[1] ## é ancora una lista
nomi_ind_new<-nomi_ind[[1]][] ### ora é vettore
### siccome con extract.haps ogni individuo é su due righe, devo toglier il nome x ciascuno nella seconda riga, con il ciclo seguente
for (i in 1:length(nomi_ind_new)) {
if (i %% 2 == 0) {nomi_ind_new[i]=c("")}
}
nomi_ind_new<-as.matrix(nomi_ind_new, ncol=1, nrow=length(nomi_ind_new))


### ciclo per mettere gli 1 nel file arp
ddd<-matrix(rep("Na",length(nomi_ind_new)), ncol=1, nrow=length(nomi_ind_new))#### aggiungo frequenza per file arp
for (i in 1:length(nomi_ind_new))   
{
if (i %% 2 == 0) {ddd[i]=c("")}
else {ddd[i]=1}
}
ddd<-as.matrix(ddd, ncol=1, nrow=length(ddd))

dat_final<-cbind(nomi_ind_new,ddd,mat_temp)
#####################################################################

colonne<-dim(dat_final)[2]
n_pop<-length(lista_pop)

con <- file(nome_output, "w")
cat("[Profile]", file=con,sep="\n")
cat("  Title = \"\"", file=con,sep="\n")
cat("  NbSamples = ", file=con)
cat(n_pop, file=con,sep="\n")
cat("  DataType = DNA", file=con,sep="\n")
cat("  GenotypicData = 1", file=con,sep="\n")
cat("  LocusSeparator = WHITESPACE", file=con,sep="\n")
cat("  MissingData = \".\"", file=con,sep="\n")
cat("  GameticPhase = 0", file=con,sep="\n")
cat("  RecessiveData = 0", file=con,sep="\n")
cat("", file=con,sep="\n")
cat("", file=con,sep="\n")
cat("", file=con,sep="\n")
cat("  [Data]", file=con,sep="\n")
cat("    [[Samples]]", file=con,sep="\n")

for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])*2
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)*2
born_sup=sum(vett_lungh_sup)*2
}
##
dat_temp<-dat_final[born_inf:born_sup,]
nome_pop<-paste("pop_", z, sep="")
cat("SampleName =", file=con)
cat("\"", file=con)
cat(nome_pop, file=con)
cat("\"", file=con,sep="\n")

cat("SampleSize = ", file=con)
cat(length(lista_pop[[z]]),  file=con,sep="\n")
cat("SampleData = {", file=con,sep="\n")


write(t(dat_temp), file=con,ncolumns=colonne, append=T)
cat("}", file=con,sep="\n")
cat("", file=con,sep="\n")
cat("", file=con,sep="\n")
}
### fine ciclo

cat("[[Structure]]", file=con,sep="\n")

cat("	StructureName=\"Simulated data\"", file=con,sep="\n")
cat("		NbGroups=1", file=con,sep="\n")

cat("			Group={", file=con,sep="\n")


for (z in 1:n_pop){
nome_pop<-paste("pop_", z, sep="")
cat("\"", file=con)
cat(nome_pop, file=con)
cat("\"", file=con,sep="\n")
}
cat("			}", file=con,sep="\n")
	

close(con)


}
########################################################################################################################################################



################# FUNZIONE 12  #############################
#################################################################################

elimina_HW<-function(data, vett_ind=c()){ ##### dati é un vcf letto da vcfR

if (length(vett_ind)==0) {
vettore_nomi<-c(names(data@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
individui<-vettore_nomi
}
else {
individui<-vett_ind
}

format<-c("FORMAT")
dati<-data[,c(format,individui)]

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = TRUE, IDtoRowNames = TRUE, extract = TRUE, convertNA = F) ### prendo i genotipi in 1/0
geno_temp<-t(genotype) ##### devo trasporre perché l'oggetto genind ha per colonne gli snps e righe gli individui
for (i in 1:ncol(geno_temp)) { ### trasformo i dati mancanti in NA perché se no non gli piace
geno_temp[geno_temp[,i]=="./.",i]<-NA 
}
mat_data=df2genind(geno_temp, ploidy = 2, sep="/") #### creo oggetto genind dove ho specificato anche le pop di appartenenza

hw_test<-hw.test(mat_data, B = 0)
hw_res<-subset(hw_test,hw_test[,3]<0.05)
nome_snp_noneq<-rownames(hw_res)

snp_togli<-c()
for (i in 1:length(nome_snp_noneq)){
temp<-which(rownames(genotype)==nome_snp_noneq[i])
snp_togli<-c(snp_togli,temp)
}

dati_clean<-dati[-snp_togli,]



output<-list(dati_clean, hw_res)

return(output)
}
###################################################################################################################




################# FUNZIONE 12a  #############################
#################################################################################

pulisci_HW_piu_pop<-function(percorso, lista_pop) {  ### dati é un vcf letto con vcfR!

data<-read.vcfR(percorso, verbose = TRUE)

nome<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(percorso))
path<-dirname(percorso)
nome_output<-paste(path,"/", nome, "_sitiHW.txt", sep="")
nome_output_vcf<-paste(path,"/", nome, "_clean.vcf.gz", sep="")

n_pop<-length(lista_pop)
format<-c("FORMAT")


snp_nonHW_totali<-c() ### qui metto i nomi SNPs eliminati in ogni popolazione
snps_eliminati<-c()  ### numero di SNPs eliminati in ogni pop
snps_tot_per_pop<-c()  ### numero di SNPs in ogni pop (esclusi i fissi)

for (r in 1:n_pop){

dati<-data[,c(format,lista_pop[[r]])]
step1<-elimina_snp_fixed_1(dati)
step2<-elimina_snp_fixed_0(step1)
hw_elim<-elimina_HW(step2)

vettore_snps<-length(getPOS(dati))

nomi_snps<-rownames(hw_elim[[2]])
snp_nonHW_totali<-c(snp_nonHW_totali,nomi_snps)

snps_eliminati<-cbind(snps_eliminati,dim(hw_elim[[2]])[1])

snps_tot_per_pop<-cbind(snps_tot_per_pop,length(getPOS(dati)))

}



vett_nomi_righe_eli<-c()
vett_nomi_righe_tot<-c()

for (k in 1:n_pop){
aaa<-paste("pop_",k,"_eli",sep="")
bbb<-paste("pop_",k,"_tot",sep="")
vett_nomi_righe_eli<-rbind(vett_nomi_righe_eli,aaa)
vett_nomi_righe_tot<-rbind(vett_nomi_righe_tot,bbb)
}

colnames(snps_eliminati)<-vett_nomi_righe_eli
colnames(snps_tot_per_pop)<-vett_nomi_righe_tot
finale_snps<-cbind(snps_eliminati,snps_tot_per_pop)


nomi_snps_totali2<-snp_nonHW_totali[!duplicated(snp_nonHW_totali)] ### metto solo i nomi SNPs singomi (per alcune pop potrei aver tolto lo stesso SNP

genotype<-extract.gt(data, element = "GT", mask = FALSE, as.numeric=F,return.alleles = TRUE, IDtoRowNames = TRUE, extract = TRUE, convertNA = F) ### prendo i genotipi in 1/0
snp_togli<-c()
for (i in 1:length(nomi_snps_totali2)){
temp<-which(rownames(genotype)==nomi_snps_totali2[i])
snp_togli<-c(snp_togli,temp)
}

dati_clean<-data[-snp_togli,]


write.table(finale_snps,nome_output)
write.vcf(dati_clean,nome_output_vcf)

#lista_finale<-list(dati_clean,finale_snps)


#return(lista_finale)

}
############################################################################################################



################# FUNZIONE 13  #############################
#################################################################################
calcola_Fis_per_pop<-function(data,lista_pop, bootstrap=1000){


format<-c("FORMAT")
dati<-data[,c(format,unlist(lista_pop))]

genotype<-extract.gt(dati, element = "GT", mask = FALSE, as.numeric=F,return.alleles = TRUE, IDtoRowNames = TRUE, extract = TRUE, convertNA = F) ### prendo i genotipi in 1/0
geno_temp<-t(genotype) ##### devo trasporre perché l'oggetto genind ha per colonne gli snps e righe gli individui
for (i in 1:ncol(geno_temp)) { ### trasformo i dati mancanti in NA perché se no non gli piace
geno_temp[geno_temp[,i]=="./.",i]<-NA 
}


n_pop<-length(lista_pop)
vett_pop<-c()
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop e creare la nuova lista di ind random per fare bootstrap
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
##
ss<-rep(z,(born_sup-born_inf+1))
vett_pop<-c(vett_pop,ss)
vett_pop<-as.factor(vett_pop)
}
### fine ciclo che mi da il vettore di fattori per identificare la pop d'appartenenza di ciascun individuo

mat_data<-df2genind(geno_temp, ploidy = 2, sep="/", pop=vett_pop) #### creo oggetto genind dove ho specificato anche le pop di appartenenza
###risultati<-basic.stats(mat_data,diploid=TRUE,digits=4)

valori<-boot.ppfis(mat_data,nboot=bootstrap,quant=c(0.025,0.975),diploid=TRUE,dig=4)
valori_CI<-valori$fis.ci


return(valori_CI)

}
#####################################################################################


################# FUNZIONE 14  #############################
##############################################################################################################################

calcola_freq_maf<-function(data,MAF=0.05){

genotype<-extract.gt(data, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0


snps_tot=dim(genotype)[1]
individui<-dim(genotype)[2]


vector_missing<-matrix(rep(0,snps_tot), ncol=snps_tot, nrow=1)

### questo primo ciclo serve per creare un vettore dove metto i missing data per snp per poter calcolare frequenze corrette
for (i in 1:snps_tot) {
temp<-0
for (j in 1:ncol(genotype)) {
if (genotype[i,j]=="./.") {temp<-temp+1}
}
vector_missing[i]<-temp
}

for (i in 1:ncol(genotype)) #coding all 0/0 as 0, 0/1 as 1 and 1/1 as 2 (0=ref, 1=alt)
{
genotype[genotype[,i]=="0/0",i]<-0 #homo anc
genotype[genotype[,i]=="0/1",i]<-1 #het
genotype[genotype[,i]=="1/1",i]<-2#homo alt
genotype[genotype[,i]=="./.",i]<-0# aggiunto da me perché ho i missing data
genotype[genotype[,i]=="./0",i]<-0# aggiunto da me perché ho i missing data
genotype[genotype[,i]=="./1",i]<-0# aggiunto da me perché ho i missing data
#print (i)
}
genotype<-matrix(as.numeric(genotype),nrow=snps_tot, ncol=individui)

#### creo vettore delle frequenze alleliche MAF 
frequenze<-matrix(rep(0,snps_tot), ncol=snps_tot, nrow=1)
for (i in 1:snps_tot) {
frequenze[i]<-sum(genotype[i,])/((individui-vector_missing[i])*2)
if (frequenze[i]>0.5) {frequenze[i]<-1-frequenze[i]} ### x come l'ho scritto, non ci sarebbe bisogno di questo controllo. ma meglio essere sicuri
}

elimina_snps<-c()  #### 
for (i in 1:snps_tot) {
if (frequenze[i]<=MAF) {elimina_snps<-rbind(elimina_snps,i)} 
}


if (length(elimina_snps)>0) {
dati_clean<-dati[-elimina_snps,]
}
else {
dati_clean<-dati
}

lista_finale<-list(frequenze,dati_clean)

return(lista_finale)

}

##########################################################################################################################

################# FUNZIONE 15  #############################
##############################################################################################################################

calcola_normalized_foldedSFS<-function(vettore_sfs){
somma<-sum(vettore_sfs)
eta_2<-vettore_sfs/somma # normalizzo per numero di SNPs
ind<-length(eta_2)

### trasformo per ottenere una curva flat per una pop costante (formula Lapierre et al. 2017, Genetics)
for (i in 1:(ind-1)){
eta_2[i]<-eta_2[i]*i*((2*ind)-i)/(2*ind)
}
eta_2[ind]<-eta_2[ind]*ind
eta_2<-as.vector(eta_2)

asse_x<-as.vector(seq(1,ind, by=1)/(2*ind)) #### calcolo gli i/2N per plottare il folded normalizzato
ssss=as.matrix(t(rbind(asse_x,eta_2)))

return(ssss)
}
##########################################################################################################################





################# FUNZIONE 16  #############################
##############################################################################################################################
library(rlist)



crea_reads_simcoal<-function(percorso, lunghezza_read=110, cov_mean=10, cov_sd=5, err_illumina=0.001, versione="fsc26") {  ###

nome<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(percorso))
path<-dirname(percorso)

sequenze<-grep("_",readLines(percorso), value=T)

#### estraggo solo la matrice dei dati da un file .arp: é la parte variabile delle future reads

n_ind<-length(sequenze)
if (versione=="fsc26") {temp_seq=strsplit(sequenze, split="\t")}
else {temp_seq=strsplit(sequenze, split=" ")}

matrice_nuc<-c()
for (i in 1:n_ind){
campo<-length(temp_seq[[i]])
matrice_nuc<-rbind(matrice_nuc,temp_seq[[i]][campo])
}
###############


####### le reads devono essere diciamo di n_base=110. creo la parte costante delle future reads 

n_base<-lunghezza_read ### da spostare come argomento se faccio funzione
da_estrarre<-length(strsplit(matrice_nuc[1,], split="")[[1]])
#vett_freq<-rmultinom(1, n_base-da_estrarre, prob=c(0.25,0.25,0.25,0.25)) ### ogni base ha 1/4 di prob ovviamente
#parte_conservata<-c(rep("A",vett_freq[1]), rep("G",vett_freq[2]), rep("C",vett_freq[3]), rep("T",vett_freq[4]) )
parte_conservata<-c()
basi<-c("A","G","C","T")
for (i in 1:(n_base-da_estrarre)){
a=sample(basi,1)
parte_conservata<-c(parte_conservata,a)
}
#######################################


###### ora produco le reads per simulare coverage. per ogni cromosoma (e non individuo, quindi per ogni linea della matrice_nuc)
###### replico delle reads, il cui numero sarà estratto da una normale. Se voglio simulare un coverage medio di 40X, ogni allele
###### avra un coverage di 20X. Metto una sd di 10, tutti e due parametri da mettere nella funzione eventuale
###### creo una matrice per individuo, prendendo a caso due linee di matrice_nuc per replicarle
###### conviene quindi scrivere per ogni locus n_ind nuove matrici, perché le devo salvare in file diversi ed é + pratico
###### simulo almeno 1 reads per locus, per una questione di semplicità con gli indici

mean_cov<-cov_mean
sd_cov<-cov_sd



### creo locus con tutte le reads
locus<-c()
indici_alleli<-c() ### contiene quante reads ci sono per individuo
for (i in 1:nrow(matrice_nuc)){
aggiungi<-floor(abs(rnorm(1,mean_cov,sd_cov)))
if (aggiungi==0){aggiungi=1}
#print(aggiungi)
indici_alleli<-c(indici_alleli, aggiungi)
for (j in 1:aggiungi){
locus<-rbind(locus,matrice_nuc[i,])
}
}

######


### mi scrivo gli indici per sapere a che cromosoma corrisponde ciascuna reads
born_inf<-c()
born_sup<-c()
for (z in 1:n_ind){
####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'é in almeno un individuo
if (z==1) {
born_inf=cbind(born_inf,1)
born_sup=cbind(born_sup,indici_alleli[z])
}
else {

born_inf=cbind(born_inf,(sum(indici_alleli[1:(z-1)])+1))
born_sup=cbind(born_sup,(sum(indici_alleli[1:z])))
}
}

####


#### creo la lista dove ogni elemento é un individo (i due alleli son stati creati indipendentemente)
finale<-list()
for (z in 1:(n_ind/2)){
mat_temp<-locus[born_inf[(z*2-1)]:born_sup[(z*2)],]
finale<-list.append(finale,mat_temp)
}
##### finale é una lista con n_ind/2 individui (n_ind é il numero di cromosomi).

###########
###########
###########

########################################################## aggiungo parte conservata a ciascun locus +  mutazioni
reads_totali_per_ind<-list()

for (m in 1:length(finale)){



reads_intera<-c()
for (i in 1:length(finale[[m]])){
reads_intera<-rbind(reads_intera,strsplit(finale[[m]][i], split="")[[1]])
}


matrice_conservata<-c()
for (i in 1:length(finale[[m]])){
matrice_conservata<-rbind(matrice_conservata,parte_conservata)
}

locus_finale<-cbind(reads_intera,matrice_conservata) ##### creo data frame reads

########## aggiungo mutazioni randon, simulando errore Illumina

epsilon=err_illumina #### il valore é espresso come probabilità per base. posso metterlo come parametro dell'eventuale funzione
n_tot_mut<-floor(dim(locus_finale)[1]*dim(locus_finale)[2]*epsilon)
if (n_tot_mut==0) {n_tot_mut=1}

for (z in 1:n_tot_mut){

i<-floor(runif(1,1,dim(locus_finale)[1]))
j<-floor(runif(1,1,dim(locus_finale)[2]))

basi<-c("A","G","C","T")
a=sample(basi,1)
while (a == locus_finale[i,j]) {
a=sample(basi,1)
}

locus_finale[i,j]=a
}


reads_totali_per_ind<-list.append(reads_totali_per_ind,locus_finale)
}

########  reads_totali_per_ind  é una lista dove ogni elemento sono le reads con mut per individuo da scrivere in file fastq


####### scrivi ogni elemento della lista in un file diverso come fosse una vera read in fastq

#### scrivo nome dei file in una nuova lista
nome_output<-list()
for (i in 1:(n_ind/2)) {
nome_output<-list.append(nome_output,paste(path,"/", nome, "_", i , ".fq", sep=""))
}

for (i in 1:(n_ind/2)) {
con <- file(nome_output[[i]], "w")

for (j in 1:nrow(reads_totali_per_ind[[i]])){
cat("@", file=con)
cat(rnorm(1,12,13), file=con)
cat(rnorm(1,12,13), file=con)
cat("/1", file=con, sep="\n")
cat(reads_totali_per_ind[[i]][j,], file=con, sep="")
cat("", file=con, sep="\n")
cat("+", file=con, sep="\n")
cat(rep("F",n_base), file=con, sep="")
cat("", file=con, sep="\n")
}
close(con)
}
}

#########################################################################################################################

################# FUNZIONE 17  #############################
##############################################################################################################################

mpd_from_sfs<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 1:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

################# FUNZIONE 17a  #############################
##############################################################################################################################

mpd_from_sfs_singleton<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 2:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

################# FUNZIONE 17b  #############################
##############################################################################################################################

mpd_from_sfs_doubleton<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 3:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

################# FUNZIONE 20  ############################# calcola la D di tajima a partire d'un SFS folded
##############################################################################################################################
calcola_TD_folded<-function(folded_sfs){  #### ha bisogno della funzione calcola_MPD
sample_size<-2*length(folded_sfs)
theta_P<-mpd_from_sfs(folded_sfs)
S<-sum(folded_sfs)

###### 
calcola_a1<-function(sample_size){
a<-0
for (i in 1:(sample_size-1)){
a<-a+(1/i)
}
return(a)
}
##################

###### 
calcola_a2<-function(sample_size){
a<-0
for (i in 1:(sample_size-1)){
a<-a+(1/(i*i))
}
return(a)
}
##################


###### 
calcola_b1<-function(sample_size){

a<-(sample_size+1)/(3*(sample_size-1))

return(a)
}
##################

###### 
calcola_b2<-function(sample_size){

a<-(2*((sample_size*sample_size)+sample_size+3))/(9*sample_size*(sample_size-1))

return(a)
}
##################

###### 
calcola_c1<-function(sample_size){
a<-calcola_b1(sample_size)-(1/calcola_a1(sample_size))
return(a)
}
##################


###### 
calcola_c2<-function(sample_size){
a<-calcola_b2(sample_size)-((sample_size+2)/(calcola_a1(sample_size)*sample_size))+(calcola_a2(sample_size)/(calcola_a1(sample_size)^2))

return(a)
}
##################

calcola_e1<-function(sample_size){
a<-calcola_c1(sample_size)/calcola_a1(sample_size)
return(a)
}

calcola_e2<-function(sample_size){
a<-calcola_c2(sample_size)/((calcola_a1(sample_size)^2)+calcola_a2(sample_size))
return(a)
}

TD<-(theta_P-(S/calcola_a1(sample_size)))/(((calcola_e1(sample_size)*S)+(calcola_e2(sample_size)*S*(S-1)))^0.5)
risultati<-list(TD,theta_P,S/calcola_a1(sample_size),S)
return(risultati)
}
#########################################################################################################################


################# FUNZIONE 21  ############################# calcola la Fst di Reynolds et al. 1983
##############################################################################################################################
fst_reynolds_sfs2D<-function(sfs_2D,maff){ ### input uno spettro 2D e un valore di maf. Restituisce un vettore con Fst di Reynolds et al. 1983(unweigheted, weighted)
dati<-sfs_2D
maf<-maff ### definisco maf
n_pop1<-(dim(dati)[2]-1)/2 ## diploidi
n_pop2<-(dim(dati)[1]-1)/2 ## diploidi
somma_as<-0
somma_bs<-0
fst_unw<-c()
for (i in 1:nrow(dati)){
for (j in 1:ncol(dati)){

if (((i+j-2)/(2*(n_pop1+n_pop2)))>=maf){
freq_min_pop1<-(j-1)/(2*n_pop1)
freq_min_pop2<-(i-1)/(2*n_pop2)
freq_min_all<-(i+j-2)/(2*(n_pop1+n_pop2))

bs<-((n_pop1*2*freq_min_pop1*(1-freq_min_pop1))+(n_pop2*2*freq_min_pop2*(1-freq_min_pop2)))/(n_pop1+n_pop2-1)
as<-((4*n_pop1*((freq_min_pop1-freq_min_all)^2))+(4*n_pop2*((freq_min_pop2-freq_min_all)^2))-bs)/(2*((2*n_pop1*n_pop2)/(n_pop1+n_pop2)))

somma_as<-somma_as+(as*dati[i,j])
somma_bs<-somma_bs+(bs*dati[i,j])
if (as+bs==0) {
temp_fst=NA
} else {
temp_fst<-as/(as+bs)
}
fst_unw<-c(fst_unw, rep(temp_fst,dati[i,j]))
}

}
}
fst<-somma_as/(somma_as+somma_bs) # weighted
fst_unw_finale<-mean(fst_unw, na.rm=T) # unweighted
risultato<-c(fst_unw_finale,fst)
return(risultato)
}
####################################################################

################# FUNZIONE 22  ############################# scrivi file GENO
############################################################################################################################## 
fun_geno<-function(data){
for (i in 1:length(data)) {
if (data[i]=="./.") {data[i]<-9}
if (data[i]=="0/0") {data[i]<-0}
if (data[i]=="0/1") {data[i]<-1}
if (data[i]=="1/1") {data[i]<-2}
}
return(data)
}
##############################################################################################################################


################# FUNZIONE 22a ############################# scrivi file GENO mod per NaN
############################################################################################################################## 
fun_geno_mod<-function(data){
for (i in 1:length(data)) {
if (data[i]=="./.") {data[i]<-NaN}
if (data[i]=="0/0") {data[i]<-0}
if (data[i]=="0/1") {data[i]<-1}
if (data[i]=="1/1") {data[i]<-2}
}
return(data)
}
##############################################################################################################################

################# FUNZIONE 22b ############################# scrivi file GENO mod per NaN
############################################################################################################################## 
fun_geno_mod_na<-function(data){
for (i in 1:length(data)) {
if (is.na(data[i])) {data[i]<-NaN}
if (data[i]=="0/0") {data[i]<-0}
if (data[i]=="0/1") {data[i]<-1}
if (data[i]=="1/1") {data[i]<-2}
}
return(data)
}
##############################################################################################################################


################# FUNZIONE 22 c  ############################# scrivi file GENO
############################################################################################################################## 
fun_conta_genotipi<-function(data){
a=0
b=0
c=0
d=0
for (i in 1:length(data)) {
if (data[i]=="0/0") {a<-a+1}
else if (data[i]=="0/1") {b<-b+1}
else if (data[i]=="1/1") {c<-c+1}
else if (data[i]=="./.") {d<-d+1}
}
risultato<-c(a,b,c,d)
return(risultato)
}
##############################################################################################################################



################# FUNZIONE 23  ############################# conta allele ALT in ogni popolazione
############################################################################################################################## 
################## conta allele ALT in ogni pop. Argomento: lista pop. Oggetto: matrice geno (0,1,2)
fun_conta_ALT_2<-function(data, lista_pop){
##################  prendo gli estremi di ogni popolazione per poter poi calcolare i dati mancanti per pop per filtrare dopo
born_inf_fin<-c()
born_sup_fin<-c()
n_pop<-length(lista_pop)
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'é in almeno un individuo
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
born_inf_fin<-c(born_inf_fin,born_inf) ## i due vettori con gli indici degli estremi delle popolazioni
born_sup_fin<-c(born_sup_fin,born_sup)
}
####################################
risultato<-c()
#for (i in 1:length(data)) {
for (j in 1:length(born_inf_fin)) {
risultato<-c(risultato, sum(data[born_inf_fin[j]:born_sup_fin[j]], na.rm=T))
}
risultato<-c(risultato,sum(data, na.rm=T))
#}
return(risultato)
}
#################################### FINE FUNZ CONTA ALT ####
#########################################################################################################################################

################# FUNZIONE 24  ############################# calcola 2S-SFS pairwise per fastsimcoal. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_sfs2D_pairwise<-function(dati,lista_pop){
n_ind<-2*length(unlist(lista_pop)) ### num of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF

########## creo matrici di output
new_list<-list()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
righe=2*length(lista_pop[[i]])+1
colonne=2*length(lista_pop[[j]])+1
matrice<-matrix(rep(0,righe*colonne),nrow=righe, ncol=colonne)
print(i)
print(j)
new_list<-list.append(new_list,matrice)
}
}
}
########

for (z in 1:nrow(mat_con_alt)){
if (mat_con_alt[z,ncol(mat_con_alt)]<(n_ind/2)) { ### il MAF é l'allele ALT
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else if (mat_con_alt[z,ncol(mat_con_alt)]>(n_ind/2)) { ### il MAF é l"allele REF
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else {  ### i due alleli hanno freq 0.5 nella pop totale, quindi seguo Laurent che aumenta di 0.5 la entries per entrambi
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+0.5
new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+0.5
indice_list<-indice_list+1
}
}
}
}
}


return(new_list)
}
####################################################################################################################################################################################


################# FUNZIONE 25  ############################# calcola fst di reynolds pairwise. Ha busogno di matrice prodotta da FUNZIONE 23 come input.
##############################################################################################
calcola_fst_da_geno<-function(ogg_geno, lista_pop, maf){
########## creo matrici di output
new_list<-list()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
righe=2*length(lista_pop[[i]])+1
colonne=2*length(lista_pop[[j]])+1
matrice<-matrix(rep(0,righe*colonne),nrow=righe, ncol=colonne)
print(i)
print(j)
new_list<-list.append(new_list,matrice)
}
}
}
########
n_ind<-ncol(ogg_geno)-1
for (z in 1:nrow(ogg_geno)){
if (ogg_geno[z,ncol(ogg_geno)]<(n_ind/2)) { ### il MAF é l'allele ALT
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]<-new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else if (ogg_geno[z,ncol(ogg_geno)]>(n_ind/2)) { ### il MAF é l"allele REF
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else {  ### i due alleli hanno freq 0.5 nella pop totale, quindi seguo Laurent che aumenta di 0.5 la entries per entrambi
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]+0.5
new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]<-new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]+0.5
indice_list<-indice_list+1
}
}
}
}
}
fst_parwise_obs<-c()
for (i in 1: length(new_list)){
temp=fst_reynolds_sfs2D(new_list[[i]],maf)[2]
fst_parwise_obs<-c(fst_parwise_obs,temp)
}
return(fst_parwise_obs)
}
##############################################################################################""


################# FUNZIONE 26  ############################# calcola fst pairwise e bootstrap. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_fst_pairwise_bootstrap<-function(dati,lista_pop, maf, boot){
n_ind<-2*length(unlist(lista_pop)) ### num of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
colnames(genotype_num)<-unlist(lista_pop)
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
pairwise_fst_obs<-calcola_fst_da_geno(mat_con_alt, lista_pop, maf)

fst_bootstrapped<-c()
for (z in 1: boot) {
pairwise_fst_temp<-c()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_lista_pop_temp<-sample(unlist(list(lista_pop[[i]],lista_pop[[j]])))
new_lista_pop<-list(new_lista_pop_temp[1:length(lista_pop[[i]])],new_lista_pop_temp[(length(lista_pop[[i]])+1):length(new_lista_pop_temp)])
print(new_lista_pop)
mat_con_alt_temp<-t(apply(genotype_num[,new_lista_pop_temp], 1, fun_conta_ALT_2, new_lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
temp<-calcola_fst_da_geno(mat_con_alt_temp, new_lista_pop, maf)###metti maf
pairwise_fst_temp<-c(pairwise_fst_temp, temp)
}
}
}
fst_bootstrapped<-rbind(fst_bootstrapped,pairwise_fst_temp)
}

lista_finale<-list(pairwise_fst_obs, fst_bootstrapped)
return(lista_finale)
}
####################################################################################################################################################################################


################# FUNZIONE 26 con is.na per STACKS  ############################# calcola fst pairwise e bootstrap. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_fst_pairwise_bootstrap_na<-function(dati,lista_pop, maf, boot){
n_ind<-2*length(unlist(lista_pop)) ### num of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod_na)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
colnames(genotype_num)<-unlist(lista_pop)
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
pairwise_fst_obs<-calcola_fst_da_geno(mat_con_alt, lista_pop, maf)

fst_bootstrapped<-c()
for (z in 1: boot) {
pairwise_fst_temp<-c()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_lista_pop_temp<-sample(unlist(list(lista_pop[[i]],lista_pop[[j]])))
new_lista_pop<-list(new_lista_pop_temp[1:length(lista_pop[[i]])],new_lista_pop_temp[(length(lista_pop[[i]])+1):length(new_lista_pop_temp)])
print(new_lista_pop)
mat_con_alt_temp<-t(apply(genotype_num[,new_lista_pop_temp], 1, fun_conta_ALT_2, new_lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
temp<-calcola_fst_da_geno(mat_con_alt_temp, new_lista_pop, maf)###metti maf
pairwise_fst_temp<-c(pairwise_fst_temp, temp)
}
}
}
fst_bootstrapped<-rbind(fst_bootstrapped,pairwise_fst_temp)
}

lista_finale<-list(pairwise_fst_obs, fst_bootstrapped)
return(lista_finale)
}
####################################################################################################################################################################################






################# FUNZIONE 26 a ############################# calcola fst pairwise senza bootstrap. PER SLIDING WINDOWS. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_fst_pairwise_nobootstrap<-function(dati,lista_pop, maf){
n_ind<-2*length(unlist(lista_pop)) ### num of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
colnames(genotype_num)<-unlist(lista_pop)
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
pairwise_fst_obs<-calcola_fst_da_geno(mat_con_alt, lista_pop, maf)

return(pairwise_fst_obs)
}
####################################################################################################################################################################################

################# FUNZIONE 27  ############################# SLIDING WINDOWS Fst pairwise e theta x tutte pop de lista_pop
############################################################################################################################## 
sliding_fst_spectre<-function(dati, lista_pop, window, slide, maf){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
n_pop<-length(lista_pop)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

ris_fst_sliding<-c()
ris_spectre_sliding<-c()
for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-dati[sss,]
fst_temp<-calcola_fst_pairwise_nobootstrap(dati_temp,lista_pop, maf) ### METTERE MAF IN FUNZIONE!
ris_fst_sliding<-rbind(ris_fst_sliding, fst_temp)

temp_spectre<-c()
for (z in 1:length(lista_pop)){
dati_bin<-vcfR2DNAbin(dati_temp[, c("FORMAT",lista_pop[[z]])], extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
temp_spectre<-c(temp_spectre,sfs_tot)
}
ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
} else {
ris_fst_sliding<-rbind(ris_fst_sliding, rep(NaN,n_pop*(n_pop-1)/2))
ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,length(unlist(lista_pop))))
}
}

finale_fst<-cbind(per_plot,ris_fst_sliding)
finale_fst<-na.omit(finale_fst)

finale_spectre<-cbind(per_plot,ris_spectre_sliding)
finale_spectre<-na.omit(finale_spectre)

lista_risultati<-list(finale_fst, finale_spectre)

return(lista_risultati)
}
##############################################################################################################################



################# FUNZIONE 27 a  ############################# SLIDING WINDOWS Fst pairwise e theta x tutte pop de lista_pop
############################################################################################################################## 
sliding_spectre_onepop<-function(dati, window, slide){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
n_ind<-length(vettore_nomi)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

ris_spectre_sliding<-c()
for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-dati[sss,]
temp_spectre<-c()
dati_bin<-vcfR2DNAbin(dati_temp, extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
temp_spectre<-c(temp_spectre,sfs_tot)
ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
} else {
ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,n_ind))
}
}

finale_spectre<-cbind(per_plot,ris_spectre_sliding)
finale_spectre<-na.omit(finale_spectre)

return(finale_spectre)
}
##############################################################################################################################


################# FUNZIONE 27 B ############################# SLIDING WINDOWS Fst pairwise e theta x tutte pop de lista_pop
###### in finale_spectre la 1 colonna e il punto mediano della finestra, la 2 é il numero di basi reali presenti. poi cominciano gli spettri
############################################################################################################################## 
sliding_fst_spectre_stand<-function(dati, lista_pop, window, slide, maf, posizioni){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
n_pop<-length(lista_pop)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

## sliding fst
ris_fst_sliding<-c()
ris_spectre_sliding<-c()
temp_siti_tot<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-dati[sss,]
fst_temp<-calcola_fst_pairwise_nobootstrap(dati_temp,lista_pop, maf) ### METTERE MAF IN FUNZIONE!
ris_fst_sliding<-rbind(ris_fst_sliding, fst_temp)
##

## sliding siti totali considerati
sss_scale=which(posizioni > low_bound[i] & posizioni < upper_bound[i] )
temp_siti_tot<-rbind(temp_siti_tot,length(sss_scale))
##

## sliding spectre
temp_spectre<-c()
for (z in 1:length(lista_pop)){
dati_bin<-vcfR2DNAbin(dati_temp[, c("FORMAT",lista_pop[[z]])], extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
temp_spectre<-c(temp_spectre,sfs_tot)
}
ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
##


} else {
ris_fst_sliding<-rbind(ris_fst_sliding, rep(NaN,n_pop*(n_pop-1)/2))
ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,length(unlist(lista_pop))))
temp_siti_tot<-rbind(temp_siti_tot, NaN)
}
}

finale_fst<-cbind(per_plot,ris_fst_sliding)
finale_fst<-na.omit(finale_fst)

finale_spectre<-cbind(per_plot,temp_siti_tot,ris_spectre_sliding)
finale_spectre<-na.omit(finale_spectre)

lista_risultati<-list(finale_fst, finale_spectre)

return(lista_risultati)
}
##############################################################################################################################

################# FUNZIONE 27 C  ############################# SLIDING WINDOWS Fst pairwise e theta x tutte pop de lista_pop
###### in finale_spectre la 1 colonna e il punto mediano della finestra, la 2 é il numero di basi reali presenti. poi cominciano gli spettri
############################################################################################################################## 
sliding_spectre_onepop_stand<-function(dati, window, slide, posizioni){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
n_ind<-length(vettore_nomi)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

ris_spectre_sliding<-c()
temp_siti_tot<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {

## sliding siti totali considerati
sss_scale=which(posizioni > low_bound[i] & posizioni < upper_bound[i] )
temp_siti_tot<-rbind(temp_siti_tot,length(sss_scale))
##

dati_temp<-dati[sss,]
temp_spectre<-c()
dati_bin<-vcfR2DNAbin(dati_temp, extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
temp_spectre<-c(temp_spectre,sfs_tot)
ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
} else {
ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,n_ind))
temp_siti_tot<-rbind(temp_siti_tot, NaN)
}
}

finale_spectre<-cbind(per_plot,temp_siti_tot,ris_spectre_sliding)
finale_spectre<-na.omit(finale_spectre)

return(finale_spectre)
}
##############################################################################################################################


################# FUNZIONE 28  ############################# Estrai un locus ogni xxx basi. restituisce vcf clean e numero di basi chiamate totali
############################################################################################################################## 
################################################################################################################################# dati = object vcfR
### window = size of the locus (for example 1000 bp)
### spazio = size of the distanze from the next locus (for exampke 100000 bp)
### posizioni = like always, the vector of all called positions to standardize the theta
### results = list with the first object: a new vcfR with the extracted positions in the loci previously defined; second: the lenght of the total bp called 
### after the fonction, three lines to show how to use it

estrai_finestre_noLD<-function(dati, window, spazio, posizioni){

dim_locus<-window
interval<-spazio
vettore_snps_pos<-getPOS(dati)
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
n_ind<-length(vettore_nomi)

####### prendo limiti inf e sup di ciascuna finestra, assicurandomi di essere all'interno della lunghezza del cromosoma #####
n_finestre<-floor((vettore_snps_pos[length(vettore_snps_pos)])/(dim_locus+interval))
inf_value<-c()
sup_value<-c()
for (i in 1:n_finestre){
temp_inf<-1+(i-1)*interval
temp_sup<-dim_locus+(i-1)*interval
inf_value<-c(inf_value,temp_inf)
sup_value<-c(sup_value,temp_sup)
}
#######################################################################################################################

siti_da_estrarre<-c()
temp_siti_tot<-c()

for (i in 1:length(inf_value)) {
sss=which(vettore_snps_pos > inf_value[i] & vettore_snps_pos < sup_value[i] )
if (length(sss)>0) {
## sliding siti totali considerati
sss_scale=which(posizioni > inf_value[i] & posizioni < sup_value[i] )
temp_siti_tot<-rbind(temp_siti_tot,length(sss_scale))
siti_da_estrarre<-c(siti_da_estrarre,sss)
} else {
temp_siti_tot<-rbind(temp_siti_tot,NaN)
}
}

dati_noLD<-dati[siti_da_estrarre,] #### nuovo vcf
lunghezza_tot<-sum(temp_siti_tot, na.rm=T) #### lunghezza tot prendendo in conto la referenza

da_salvare<-list(dati_noLD, lunghezza_tot)

return(da_salvare)

}
#######################################################################################################################
#######################################################################################################################


################# FUNZIONE 29 ############################# SLIDING WINDOWS to extract SNP all het in pop 1 and all fixed in pop 2
#### BE CAREFUL ABOUT THE ORDER OF POP!!! IT WILL COMPUTE THE HET in POP1 and the HOMO in POP2. So be careful to the order in lista_pop
############################################################################################################################## 
sliding_het_fix<-function(dati, lista_pop, window, slide){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
n_pop<-length(lista_pop)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

dati_temp2<-c()
for (i in 1:length(lista_pop)){
dati_moment<-dati[,c("FORMAT", lista_pop[[i]])]
genotype<-extract.gt(dati_moment, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
aaa=t(apply(genotype, 1, fun_conta_genotipi))
dati_temp2<-cbind(dati_temp2,aaa)
}


n_siti_h<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-matrix(dati_temp2[sss,], ncol=ncol(dati_temp2), nrow=length(sss))
snp=0
for (h in 1:nrow(dati_temp)){
if (dati_temp[h,2] == (length(lista_pop[[1]])-dati_temp[h,4]) && (dati_temp[h,5] + dati_temp[h,7])  == (length(lista_pop[[2]])-dati_temp[h,8]) ) {snp=snp+1} 
}
n_siti_h<-c(n_siti_h,snp)
} else {n_siti_h<-c(n_siti_h,NA)}

}

finale<-cbind(per_plot,n_siti_h)
finale<-na.omit(finale)


return(finale)
}
##############################################################################################################################

################# FUNZIONE 29 B ############################# SLIDING WINDOWS to extract SNP all het in pop 1 and all fixed in pop 2
#### BE CAREFUL ABOUT THE ORDER OF POP!!! IT WILL COMPUTE THE HET in POP1 and the HOMO in POP2. So be careful to the order in lista_pop
############################################################################################################################## 
sliding_het_fix_soloref<-function(dati, lista_pop, window, slide){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
n_pop<-length(lista_pop)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

dati_temp2<-c()
for (i in 1:length(lista_pop)){
dati_moment<-dati[,c("FORMAT", lista_pop[[i]])]
genotype<-extract.gt(dati_moment, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
aaa=t(apply(genotype, 1, fun_conta_genotipi))
dati_temp2<-cbind(dati_temp2,aaa)
}


n_siti_h<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-matrix(dati_temp2[sss,], ncol=ncol(dati_temp2), nrow=length(sss))
snp=0
for (h in 1:nrow(dati_temp)){
if (dati_temp[h,2] == (length(lista_pop[[1]])-dati_temp[h,4]) && dati_temp[h,5]  == (length(lista_pop[[2]])-dati_temp[h,8]) ) {snp=snp+1} 
}
n_siti_h<-c(n_siti_h,snp)
} else {n_siti_h<-c(n_siti_h,NA)}

}

finale<-cbind(per_plot,n_siti_h)
finale<-na.omit(finale)


return(finale)
}
##############################################################################################################################


################# FUNZIONE 29 C ############################# SLIDING WINDOWS to extract SNP all het in pop 1 and all fixed in pop 2
####  IT WILL COMPUTE THE SNPs fixed between POP1 & POP2
############################################################################################################################## 
###################################################################################################################################################
sliding_SNP_homo_fix<-function(dati, lista_pop, window, slide){
dim_locus<-window
interval<-slide
vettore_snps_pos<-getPOS(dati)
n_pop<-length(lista_pop)
low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
per_plot<-(upper_bound+low_bound)/2

dati_temp2<-c()
for (i in 1:length(lista_pop)){
dati_moment<-dati[,c("FORMAT", lista_pop[[i]])]
genotype<-extract.gt(dati_moment, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
aaa=t(apply(genotype, 1, fun_conta_genotipi))
dati_temp2<-cbind(dati_temp2,aaa)
}


n_siti_h<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>0) {
dati_temp<-matrix(dati_temp2[sss,], ncol=ncol(dati_temp2), nrow=length(sss))
snp_a=0
snp_b=0
for (h in 1:nrow(dati_temp)){
if (dati_temp[h,1] == (length(lista_pop[[1]])-dati_temp[h,4]) && dati_temp[h,7]  == (length(lista_pop[[2]])-dati_temp[h,8]) ) {snp_a=snp_a+1} 
if (dati_temp[h,3] == (length(lista_pop[[1]])-dati_temp[h,4]) && dati_temp[h,5]  == (length(lista_pop[[2]])-dati_temp[h,8]) ) {snp_b=snp_b+1} 
}
snp<-snp_a+snp_b
n_siti_h<-c(n_siti_h,snp)
} else {n_siti_h<-c(n_siti_h,NA)}

}

finale<-cbind(per_plot,n_siti_h)
finale<-na.omit(finale)


return(finale)
}
###################################################################################################################################################
###WINGS<- c("13_Q_A_S12","14_Q_A_S13","15_Q_A_S14","16_Q_A_S15","17_Q_A_S16","18_Q_A_S17","19_Q_A_S18", "20_Q_A_S19","21_Q_A_S20","22_Q_A_S21","23_Q_A_S22","28_Q_A_S24", "2_Q_S_A_S2","6_Q_S_A_S6") 
###NOWINGS<-c("1_Q_S_A_S1","3_Q_S_A_S3","4_Q_S_A_S4","5_Q_S_A_S5","7_Q_S_A_S7","8_Q_S_A_S8" ,"10_Q_S_A_S9","11_Q_S_A_S10","12_Q_S_A_S11", "26_Q_A_S23" ) 
###lista_pop<-list(WINGS, NOWINGS)
###dati<-read.vcfR("C:/Users/SMona/Documents/antoine/VCF_stefano/fst/Structure_Scaffold17_GATK_TAG_Flowqual_Noindels_NoRef_Na20.vcf", verbose = TRUE)
###dati<-dati[,c("FORMAT", unlist(lista_pop))]
###risultato<-sliding_SNP_homo_fix(dati, lista_pop, 100000, 15000)
###plot(risultato, type="l")
###################################################################################################################################################


################# FUNZIONE 30 a ############################# SLIDING WINDOWS PCA
#### slides the number of SNPs, not the base pair as the previous functions!
##############################################################################################################################
sliding_pca_adegenet<-function(dati, n_snp, slide){ ### 1st arg: your good vcf (with filters for structure analyses); 2nd: the Num of SNPs in the window (not the total lenght as in the other function); 3rd: Num of SNPs to slide.

vettore_snps_pos<-getPOS(dati)
snp_totali<-length(vettore_snps_pos)
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf

low_bound<-seq(1, snp_totali, slide)
upper_bound<-low_bound+n_snp


pca_asse1<-c()
pca_asse2<-c()
finestra_low<-c()
finestra_up<-c()

for (i in 1:length(low_bound)) {

if (upper_bound[i] <= snp_totali) {

dati_moment<-dati[low_bound[i]:upper_bound[i],]
genotype<-extract.gt(dati_moment, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

geno_temp<-t(genotype) ##### devo trasporre perché l'oggetto genind ha per colonne gli snps e righe gli individui
for (j in 1:ncol(geno_temp)) { ### trasformo i dati mancanti in NA perché se no non gli piace
geno_temp[geno_temp[,j]=="./.",j]<-NA 
}


mat_data<-df2genind(geno_temp, ploidy = 2, sep="/", pop=vettore_nomi) #### creo oggetto genind dove ho specificato anche le pop di appartenenza

x<-indpca(mat_data)
new_eigen<-x$ipca$eig/sum(x$ipca$eig)
pca_asse1<-c(pca_asse1,new_eigen[1])
pca_asse2<-c(pca_asse2,new_eigen[2])

finestra_low<-c(finestra_low,vettore_snps_pos[low_bound[i]])
finestra_up<-c(finestra_up,vettore_snps_pos[upper_bound[i]])
#print(i)
}
}

finale<-cbind(finestra_low,finestra_up,pca_asse1, pca_asse2)

return(finale)
}
###dati<-read.vcfR("C:/Users/SMona/Documents/antoine/VCF_stefano/fst/Structure_Scaffold17_GATK_TAG_Flowqual_Noindels_NoRef_Na20.vcf", verbose = TRUE) 
###risultati<-sliding_pca_adegenet(dati,n_snp=1000, slide=50) 
### in risultati you have a matrix of four col: 
## col1 = lower value of the window
## col2 = upper value of the window
## col3 = % of variance of the first axe
## col4 = % of variance of the second axe
#### you may plot for example the average values between col1 & col2 (i.e.: more or less the midpoint of the window) agains either col3 or col3+col4
#### I tested on scaffold 17 and even not hyper powerful it identified what we know is the supergene
#### the slide is, contrary to the previous sliding window analayses, based on the Num of SNPs and not the base pair. So is the slide. 
#### I thought it was a good idea to have enough SNP in each window. of course I can modify it but i think is it fair
#### in the above example I put a window of 300 SNP sliding 100 SNP. You can play with the two parameters! dati is as always the vcf (the one you with the filters for structure analyses)
############################################################################################################################################################################################################################################################


################# FUNZIONE 30 b ############################# SLIDING WINDOWS PCA
#### slides the base pair as the previous functions!
##############################################################################################################################
##############################################################################################################################
sliding_pca_adegenet_basepair<-function(dati, window, slide, min_n_snp){ 
### 1st arg: your good vcf (with filters for structure analyses); 
### 2nd arg: the size of the window
### 3rd arg: the size of the slide
### 4th arg: minimum Num of SNPs to compute PCA (10 should be enough, it is jsut to avoid error at the bounds of the scaffold)

vettore_snps_pos<-getPOS(dati)
vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf

dim_locus<-window
interval<-slide
#low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
#upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
low_bound<-seq(1, max(vettore_snps_pos)-interval, interval)
upper_bound<-low_bound+window


pca_asse1<-c()
pca_asse2<-c()

for (i in 1:length(low_bound)) {
sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
if (length(sss)>min_n_snp) {
print(length(sss))
dati_moment<-dati[sss,]
genotype<-extract.gt(dati_moment, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

geno_temp<-t(genotype) ##### devo trasporre perché l'oggetto genind ha per colonne gli snps e righe gli individui
for (j in 1:ncol(geno_temp)) { ### trasformo i dati mancanti in NA perché se no non gli piace
geno_temp[geno_temp[,j]=="./.",j]<-NA 
}
mat_data<-df2genind(geno_temp, ploidy = 2, sep="/", pop=vettore_nomi) #### creo oggetto genind dove ho specificato anche le pop di appartenenza
x<-indpca(mat_data)
new_eigen<-x$ipca$eig/sum(x$ipca$eig)
pca_asse1<-c(pca_asse1,new_eigen[1])
pca_asse2<-c(pca_asse2,new_eigen[2])
} else {
pca_asse1<-c(pca_asse1,NA)
pca_asse2<-c(pca_asse2,NA)
}
#print(i)
}

finale<-cbind(low_bound,upper_bound,pca_asse1, pca_asse2)
return(finale)
}
##############################################################################################################################
###dati<-read.vcfR("C:/Users/SMona/Documents/antoine/VCF_stefano/fst/Structure_Scaffold17_GATK_TAG_Flowqual_Noindels_NoRef_Na20.vcf", verbose = TRUE) 
###risultati<-sliding_pca_adegenet_basepair(dati,window=50000, slide=10000, min_n_snp=10)
###plot(apply(risultati[,1:2],1,mean),risultati[,3])
###
### in risultati you have a matrix of four col: 
## col1 = lower value of the window
## col2 = upper value of the window
## col3 = % of variance of the first axe
## col4 = % of variance of the second axe
#### you may plot for example the average values between col1 & col2 (i.e.: more or less the midpoint of the window) agains either col3 or col3+col4
##############################################################################################################################
##############################################################################################################################

################# FUNZIONE 31 ############################# FST by SNP
##############################################################################################################################
##############################################################################################################################

fst.hudson <- function(vcf, pop_list){
  require(rlist)
  require(vcfR)
  
  gt2alleles = function(gt){
    alleles = unlist(strsplit(gt,"[/|]"))
    alleles[which(alleles==".")] = NA
    while(length(alleles) != length(gt)*2){
      alleles <- c(alleles,NA)
    }
    return(as.numeric(alleles))
  }
  
  
  mpd = function(an,ac){
    n_pairs = an * (an - 1) / 2
    n_same = rowSums(ac * (ac - 1) / 2)
    
    n_diff = n_pairs - n_same
    
    return(n_diff / n_pairs)
  }
  
  
  mpd.between = function(an1,an2,ac1,ac2){
    n_pairs = an1 * an2 
    n_same = rowSums(ac1 * ac2)
    
    n_diff = n_pairs - n_same
    
    return(n_diff / n_pairs)
  }
  
  GT = extract.gt(vcf,as.numeric = F,convertNA = F)
  
  mx = max(as.numeric(unlist(strsplit(max(GT),"[/|]"))))
  
  a.counts = function(site){
    c <- c()
    for(i in 1:(mx+1)){
      c <- c(c,length(which(site == (i-1))))
    }
    return(c)
  }
  
  alleles_counts = allele_num = list()
  for(i in 1:length(pop_list)){
    alleles = t(apply(GT[,pop_list[[i]]],1,gt2alleles))
    
    alleles_counts[[i]] = t(apply(alleles,1,a.counts))
    
    allele_num[[i]] = rowSums(alleles_counts[[i]])
  }
  
  
  fst_nuc = list()
  fst = nms = c()
  for(i in 1:(length(pop_list)-1)){
    for(j in (i+1):length(pop_list)){
      within = (mpd(allele_num[[i]],alleles_counts[[i]]) +
                  mpd(allele_num[[j]],alleles_counts[[j]]))/2
      between = mpd.between(allele_num[[i]],allele_num[[j]],alleles_counts[[i]],alleles_counts[[j]])
      fst_nuc = list.append(fst_nuc,(between - within)/between)
      fst = c(fst,(sum(between,na.rm = T) - sum(within,na.rm = T))/sum(between,na.rm = T))
      nms = c(nms,paste0(names(pop_list)[i],"/",names(pop_list)[j]))
    }
  }
  names(fst_nuc) = nms
  fst = data.frame(FST=fst)
  rownames(fst) = nms
  
  result = list(fst,fst_nuc)
  names(result) = c("pairwise_fst","pairwise_fst_nuc")
  return(result)
}