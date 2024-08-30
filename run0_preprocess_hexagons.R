# 2022 10 21 I.Zliobaite
# finally hexagons

data_sites <- read.csv('data/ISEA3H09_small/Centroids_ISEA3H09_Geodetic_V_WGS84.txt', header = TRUE, sep = "\t")
data_traits <- read.csv('data/traits_modern_historic_4_copy.csv', header = TRUE, sep = "\t")
data_clim <- read.csv('data/ISEA3H09_small/WorldClim30AS_V02/ISEA3H09_WorldClim30AS_V02_BIO_Mean.txt', header = TRUE, sep = "\t")
data_koppen <- read.csv('data/ISEA3H09_small/GLOH2O_V01/ISEA3H09_GLOH2O_V01_KoppenGeiger_Mode.txt', header = TRUE, sep = "\t")
data_den <- read.csv('data/data_densities_mam_copy.csv', header = TRUE, sep = "\t")

# this is for body mass box plots
#TRUE works only after run2_analysis.R has been executed at least once
do_subset <- FALSE 
which_subset <- 'ALL'
# 'TEMP' 'PREC' 


ind <- which(data_traits[,'Order']!="Cetacea")
data_traits <- data_traits[ind,]

density_pred <- 10^(1.317435 -0.7287368*log10(data_traits[,'MASS']))

ingestion_pred <- 0.71*data_traits[,'MASS']^0.7
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

ind_present <- which(sites_sum>=10)

pdf('plots/fig_map_sites.pdf',width = 10, height = 6)
plot(data_sites[ind_present,'X'],data_sites[ind_present,'Y'],pch=16,cex = 0.2)
dev.off()

#data_koppen

data_occ <- cbind(data_Art,data_Per,data_Pro,data_Pri)
data_occ <- data_occ[ind_present,]

data_sites <- cbind(data_sites,data_koppen)
data_sites <- data_sites[ind_present,]

prec <- round(data_clim[ind_present,'BIO12_Mean'])
temp <- round(data_clim[ind_present,'BIO01_Mean'],digits = 1)

if (do_subset){
  # plot mass box
  
  hid_all <- read.csv('data_working/hid_used_all.csv', header = TRUE, sep = "\t")
  hid_all <- hid_all[,'HID']
  ind <- which(data_sites[,'HID'] %in% hid_all)
  print('mass by used sites')
  
  #ind_temp <- which(hid_all %in% hid_temp)
  
  file_name_mass <- 'plots/fig_mass_hyp_ALL.pdf'
  file_plot_title <- 'All study sites'
  
  if (which_subset == 'TEMP'){
    hid_temp <- read.csv('data_working/hid_used_tlim.csv', header = TRUE, sep = "\t")
    ind <- which(data_sites[,'HID'] %in% hid_temp[,1])
    print('mass by temperature')
    file_name_mass <- 'plots/fig_mass_hyp_TEMP.pdf'
    file_plot_title <- 'Temperature limited sites'
    
  }
  
  if (which_subset == 'PREC'){
    hid_prec <- read.csv('data_working/hid_used_plim.csv', header = TRUE, sep = "\t")
    ind <- which(data_sites[,'HID'] %in% hid_prec[,1])
    print('mass by precipitation')
    file_name_mass <- 'plots/fig_mass_hyp_PREC.pdf'
    file_plot_title <- 'Precipitation limited sites'
  }
  
  data_occ_for_mass <- data_occ[ind,]
  print(dim(data_occ_for_mass))
  sm_occ <- apply(data_occ_for_mass,2,sum)
  ind <- which(sm_occ>0)
  print(length(ind))
  data_occ_for_mass <- data_occ[,ind]
  
  
  data_mass_hyp <- c()
  for (sk in 1:dim(data_occ_for_mass)[2]){
    ind <- which(data_traits[,'Phyname']==colnames(data_occ_for_mass)[sk])  
    if (length(ind)==0){
      print('trouble')
    }else{
      data_mass_hyp <- rbind(data_mass_hyp,data_traits[ind,c('HYP','MASS','density_pred')])
    }
  }
  
  log_mass <- log10(data_mass_hyp[,'MASS'])
  data_mass_hyp <- cbind(data_mass_hyp,log_mass)
  log_density <- log10(data_mass_hyp[,'density_pred'])
  data_mass_hyp <- cbind(data_mass_hyp,log_density)
  
  
  pdf(file_name_mass,width = 4.5, height = 5)
  boxplot(log_mass ~ HYP, data = data_mass_hyp,xlab = 'Hypsodonty category',ylab = 'Log10 individual mass, kg',ylim = c(-2,5) , xaxt="n", main = file_plot_title)
  axis(1, labels=c('Brachydont','Mesodont','Hypsodont'), at = c(1,2,3))
  dev.off()
  
}