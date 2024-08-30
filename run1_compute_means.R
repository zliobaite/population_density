# 2022 10 21 I.Zliobaite
# finally hexagons
# updated 2023 09 11

data_sites <- read.csv('data/ISEA3H09_small/Centroids_ISEA3H09_Geodetic_V_WGS84.txt', header = TRUE, sep = "\t")
data_traits <- read.csv('data/traits_modern_historic_4_copy.csv', header = TRUE, sep = "\t")
data_clim <- read.csv('data/ISEA3H09_small/WorldClim30AS_V02/ISEA3H09_WorldClim30AS_V02_BIO_Mean.txt', header = TRUE, sep = "\t")
data_realms <- read.csv('data/ISEA3H09_small/WWFTE_V02/ISEA3H09_WWFTE_V02_Realm_Mode.txt', header = TRUE, sep = "\t")
data_koppen <- read.csv('data/ISEA3H09_small/GLOH2O_V01/ISEA3H09_GLOH2O_V01_KoppenGeiger_Mode.txt', header = TRUE, sep = "\t")
data_den <- read.csv('data/data_densities_mam_copy.csv', header = TRUE, sep = "\t")

ind <- which(data_traits[,'Order']!="Cetacea")
data_traits <- data_traits[ind,]

densityTetra <- c()
for (sk in 1:dim(data_traits)[1]){
  sp_now <- data_traits[sk,'Taxon']
  ind <- which(data_den[,'taxon']==sp_now)
  if (length(ind)==0){
    densityTetra <- c(densityTetra,0)
  }else{
    densityTetra <- c(densityTetra,data_den[ind,'md_density'])
    #if (data_traits[sk,'MASS']!=data_den[ind,'MASSALL']){
    #  print('troblem')
    #  print(sp_now)
    #  print(data_traits[sk,'MASS'])
    #  print(data_den[ind,'MASSALL'])
    #}
  }
}

data_traits <- cbind(data_traits,densityTetra)

density_pred <- 10^(1.335925 -0.7309009*log10(data_traits[,'MASS']))

ingestion_pred <- 0.63*data_traits[,'MASS']^0.71
# from Shipley94 ? g/min

data_traits <- cbind(data_traits,density_pred,ingestion_pred)

data_Art <- read.csv('data/ISEA3H09_small/PHYLACINE_V010201/ISEA3H09_PHYLACINE_V010201_Present_Cetartiodactyla_Centroid.txt', header = TRUE, sep = "\t")
ind_art <- 1
print(dim(data_Art))
for (sk in 2:dim(data_Art)[2]){
  if (colnames(data_Art)[sk] %in% data_traits[,'Phyname']){
    ind_art <- c(ind_art,sk)
  }
}
data_Art <- data_Art[,ind_art]
print(dim(data_Art))
data_Per <- read.csv('data/ISEA3H09_small/PHYLACINE_V010201/ISEA3H09_PHYLACINE_V010201_Present_Perissodactyla_Centroid.txt', header = TRUE, sep = "\t")
print(dim(data_Per))
data_Pro <- read.csv('data/ISEA3H09_small/PHYLACINE_V010201/ISEA3H09_PHYLACINE_V010201_Present_Proboscidea_Centroid.txt', header = TRUE, sep = "\t")
print(dim(data_Pro))
data_Pri <- read.csv('data/ISEA3H09_small/PHYLACINE_V010201/ISEA3H09_PHYLACINE_V010201_Present_Primates_Centroid.txt', header = TRUE, sep = "\t")
print(dim(data_Pri))

#checksums
if (sum(data_Art[,1]!=data_sites[,1])){
  print('problem Artiodactyla')
}
if (sum(data_Per[,1]!=data_sites[,1])){
  print('problem Perissodactyla')
}
if (sum(data_Pro[,1]!=data_sites[,1])){
  print('problem Proboscidea')
}
if (sum(data_Pri[,1]!=data_sites[,1])){
  print('problem Primates')
}

sites_Art <- data_Art[,1]

data_Art <- data_Art[,2:dim(data_Art)[2]]
data_Per <- data_Per[,2:dim(data_Per)[2]]
data_Pro <- data_Pro[,2:dim(data_Pro)[2]]
data_Pri <- data_Pri[,2:dim(data_Pri)[2]]

names_Art <- colnames(data_Art)
names_Per <- colnames(data_Per)
names_Pro <- colnames(data_Pro)
names_Pri <- colnames(data_Pri)

sites_count_Art <- apply(data_Art,1,sum)
sites_count_Per <- apply(data_Per,1,sum)
sites_count_Pro <- apply(data_Pro,1,sum)
sites_count_Pri <- apply(data_Pri,1,sum)

sites_sum <- sites_count_Art + sites_count_Per + sites_count_Pro + sites_count_Pri

ind_present <- which(sites_sum>2)

pdf('plots/fig_map_sites_all.pdf',width = 10, height = 7)
plot(data_sites[ind_present,'X'],data_sites[ind_present,'Y'],pch=16,cex = 0.2)
dev.off()


data_occ <- cbind(data_Art,data_Per,data_Pro,data_Pri)
data_occ <- data_occ[ind_present,]

print("removing Homo sapiens")
ind <- which(colnames(data_occ)!='Homo_sapiens_Centroid')
data_occ <- data_occ[,ind]

data_HYP <- data_occ
data_DEN <- data_occ
data_DENo <- data_occ
data_MASS <- data_occ
data_log_MASS <- data_occ
data_ING <- data_occ
for (sk in 1:dim(data_HYP)[2]){
  ind <- which(data_traits[,'Phyname']==colnames(data_HYP)[sk])
  if (length(ind)==0){
    print('trouble')
    print(colnames(data_HYP)[sk])
  }else{
    data_HYP[,sk] <- data_HYP[,sk]*data_traits[ind,'HYP']
    data_DEN[,sk] <- data_DEN[,sk]*data_traits[ind,'density_pred']
    data_DENo[,sk] <- data_DENo[,sk]*data_traits[ind,'densityTetra']
    data_MASS[,sk] <- data_MASS[,sk]*data_traits[ind,'MASS']
    data_log_MASS[,sk] <- data_log_MASS[,sk]*log10(data_traits[ind,'MASS'])
    data_ING[,sk] <- data_ING[,sk]*data_traits[ind,'ingestion_pred']
  }
}

data_HYP_DEN <- data_DEN*data_HYP
data_ING_DEN <- data_DEN*data_ING
data_MASS_DEN <- data_DEN*data_MASS
data_HYP_MASS_DEN <- data_MASS_DEN*data_HYP
data_HYP_ING_DEN <- data_ING_DEN*data_HYP

data_HYP_DENo <- data_DENo*data_HYP
data_ING_DENo <- data_DENo*data_ING
data_MASS_DENo <- data_DENo*data_MASS
data_HYP_MASS_DENo <- data_MASS_DENo*data_HYP
data_HYP_ING_DENo <- data_ING_DENo*data_HYP

data_sites <- data_sites[ind_present,]
prec <- round(data_clim[ind_present,'BIO12_Mean'])
temp <- round(data_clim[ind_present,'BIO01_Mean'],digits = 1)
data_koppen <- data_koppen[ind_present,]
data_realms <- data_realms[ind_present,]
  
mean_hyp <- c()
mean_mass <- c()
mean_log_mass <- c()
mean_hyp_den <- c()
mean_hyp_mass_den <- c()
mean_hyp_ing_den <- c()
occ_count <- c()

mean_hyp_o <- c()
mean_hyp_deno <- c()
mean_hyp_mass_deno <- c()
mean_hyp_ing_deno <- c()
den_count <- c()

for (sk in 1:dim(data_HYP)[1]){
  ind <- which(data_HYP[sk,]>0)
  if (length(ind)<2){
    print('trouble')
  }
  mn <- mean(as.numeric(data_HYP[sk,ind]))
  mn_mass <- mean(as.numeric(data_MASS[sk,ind]))
  mn_log_mass <- mean(as.numeric(data_MASS[sk,ind]))
  mn_den <- sum(data_HYP_DEN[sk,ind])/sum(data_DEN[sk,ind])
  mn_mass_den <- sum(data_HYP_MASS_DEN[sk,ind])/sum(data_MASS_DEN[sk,ind])
  mn_ing_den <- sum(data_HYP_ING_DEN[sk,ind])/sum(data_ING_DEN[sk,ind])
  
  mean_hyp <- c(mean_hyp,mn)
  mean_mass <- c(mean_mass,mn_mass)
  mean_log_mass <- c(mean_log_mass,mn_log_mass)
  mean_hyp_den <- c(mean_hyp_den,mn_den)
  mean_hyp_mass_den <- c(mean_hyp_mass_den,mn_mass_den)
  mean_hyp_ing_den <- c(mean_hyp_ing_den,mn_ing_den)
  
  occ_count <- c(occ_count,length(ind))
  
  ind <- which(data_HYP_DENo[sk,]>0)
  if (length(ind)>1){
    mn_o <- mean(as.numeric(data_HYP[sk,ind]))
    mn_deno <- sum(data_HYP_DENo[sk,ind])/sum(data_DENo[sk,ind])
    mn_mass_deno <- sum(data_HYP_MASS_DENo[sk,ind])/sum(data_MASS_DENo[sk,ind])
    mn_ing_deno <- sum(data_HYP_ING_DENo[sk,ind])/sum(data_ING_DENo[sk,ind])
    den_count <- c(den_count,length(ind))
  }else{
    if (length(ind)==1){
      mn_o <- as.numeric(data_HYP[sk,ind])
      mn_deno <- mn_o
      mn_mass_deno <- mn_o
      mn_ing_deno <- mn_o
      den_count <- c(den_count,1)
    }else{
      mn_o <- NA
      mn_deno <- mn_o
      mn_mass_deno <- mn_o
      mn_ing_deno <- mn_o
      den_count <- c(den_count,0)
    }
  }
  mean_hyp_o <- c(mean_hyp_o,mn_o)
  mean_hyp_deno <- c(mean_hyp_deno,mn_deno)
  mean_hyp_mass_deno <- c(mean_hyp_mass_deno,mn_mass_deno)
  mean_hyp_ing_deno <- c(mean_hyp_ing_deno,mn_ing_deno)
}

data_sites <- cbind(data_sites,data_koppen,data_realms,prec,temp,mean_hyp,mean_hyp_den,mean_hyp_mass_den,mean_hyp_ing_den,occ_count,mean_hyp_o,mean_hyp_deno,mean_hyp_mass_deno,mean_hyp_ing_deno,den_count,mean_mass,mean_log_mass)

ind <- which(data_sites[,'prec'] > 0)
data_sites <- data_sites[ind,]

ind <- which(data_sites[,'temp'] > (-99) )
data_sites <- data_sites[ind,]

lim_temp <- 3000 / (1 + exp(1.315 - 0.119 * (data_sites[,'temp'])))
lim_prec <- 3000 * (1 - exp(-0.000664*data_sites[,'prec']))
npp <- round(apply(cbind(lim_temp,lim_prec),1,min))
prec_lim <- npp*0
ind <- which(lim_prec<lim_temp)
prec_lim[ind] <- 1

data_sites <- cbind(data_sites,npp,prec_lim)

write.table(data_sites, file = "data_working/data_meanhyp_sites.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   

pdf('fig_box_mass.pdf',width = 10, height = 7)
boxplot(mean_mass ~ prec_lim, data = data_sites)
dev.off()

pdf('fig_box_log_mass.pdf',width = 10, height = 7)
boxplot(mean_log_mass ~ prec_lim, data = data_sites)
dev.off()

ind <- which(data_sites[,'prec_lim']==0)
data_sites_temp <- data_sites[ind,]

minx <- min(data_sites_temp[,'X'])
maxx <- max(data_sites_temp[,'X'])
miny <- min(data_sites_temp[,'Y'])
maxy <- max(data_sites_temp[,'Y'])
