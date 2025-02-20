# 2023 09 13 I.Zliobaite

data_sites <- read.csv('data_working/data_meanhyp_sites.csv', header = TRUE, sep = "\t")

dd_variant <- c('all','plim','tlim','Af','Am','Aw','BSh','BSk','BWh','BWk','Cfa','Cfb','Csb','Cwa','Cwb','Dwb','Dwc','ET')

wd1 <- 4.5
ht1 <- 5

ind_trop <- which(data_sites[,'Y']< 23.4362)
ind_trop <- intersect(ind_trop,which(data_sites[,'Y'] > (-23.4362)))

pdf('plots/fig_map_sites_start.pdf',width = 10, height = 7)
plot(data_sites[,'X'],data_sites[,'Y'],pch=16,cex = 0.2)
points(data_sites[ind_trop,'X'],data_sites[ind_trop,'Y'],pch=16,cex = 0.2,col = 'red')
dev.off()

data_sites[,'Madagascar'] <- 0
ind <- which(data_sites[,'X']>42)
ind  <- intersect(ind,which(data_sites[,'X']<60))
ind  <- intersect(ind,which(data_sites[,'Y']< (-11)))
data_sites[ind,'Madagascar'] <- 1


data_sites[,'Indonesia'] <- 0
ind <- which(data_sites[,'Y']<7)
ind  <- intersect(ind,which(data_sites[,'X']>92))
ind  <- intersect(ind,which(data_sites[,'X']< 132))
data_sites[ind,'Indonesia'] <- 1


ind <- which(data_sites[,'occ_count']>=10)
data_sites <- data_sites[ind,]

#ind <- which(data_sites[,'den_count']>1)
#data_sites <- data_sites[ind,]

ind <- which(data_sites[,'Madagascar'] == 0)
data_sites <- data_sites[ind,]

ind <- which(data_sites[,'Indonesia'] == 0)
data_sites <- data_sites[ind,]

ind <- which(data_sites[,'X']> (-30))
data_sites <- data_sites[ind,]

pdf('plots/fig_map_sites_filtered.pdf',width = 10, height = 7)
plot(data_sites[,'X'],data_sites[,'Y'],pch=16,cex = 0.2)
dev.off()

ind_trop <- which(data_sites[,'Y']< 23.4362)
ind_trop <- intersect(ind_trop,which(data_sites[,'Y'] > (-23.4362)))

pdf('plots/fig_box_mass.pdf',width = 10, height = 7)
boxplot(mean_mass ~ prec_lim, data = data_sites)
dev.off()

pdf('plots/fig_box_log_mass.pdf',width = 10, height = 7)
boxplot(mean_log_mass ~ prec_lim, data = data_sites)
dev.off()

#write.table(data_sites[,1], file = "data_working/HID_npp.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
#write.table(data_sites[,1], file = "data_working/HID_prec.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   

data_sites_all <- data_sites

res_sites <-c()

for (sk in 1:length(dd_variant)){
  dd_now <- dd_variant[sk]
  
  if (dd_now=='all'){
    data_sites <- data_sites_all
    write.table(data_sites[,1:3],file = "data_working/hid_used_all.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
  }else{
    if (dd_now=='plim'){
      ind <- which(data_sites_all[,'prec_lim']==1)
      data_sites <- data_sites_all[ind,] 
      write.table(data_sites[,1:3],file = "data_working/hid_used_plim.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
    }else{
      if (dd_now=='tlim'){
        ind <- which(data_sites_all[,'prec_lim']==0)
        data_sites <- data_sites_all[ind,] 
        write.table(data_sites[,1:3],file = "data_working/hid_used_tlim.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
      }else{
        ind <- which(data_sites_all[,'KoppenGeiger_Mode']==dd_now)
        data_sites <- data_sites_all[ind,] 
      }
    }
  }
  
  print(dd_now)
  print(dim(data_sites))
  
  #NPP
  
  cc <- 1.3
  coltrans <- '#80808080'
  ll <- 2
  llci <- 0.5
  
  fit_npp <- lm(npp ~ mean_hyp,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_npp_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp'],data_sites[,'npp'],xlab = 'mean HYP sp.', ylab = 'NPP',cex=cc, main = 'Mean HYP by species',ylim = c(0,3000),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  abline(fit_npp,lwd=ll)
  tx <- paste('R2 =',round(summary(fit_npp)$r.squared,digits = 2))
  text(2.7,2600,tx,cex=cc)
  dev.off()
  
  fit_npp_den <- lm(npp ~ mean_hyp_den,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_npp_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_den'],data_sites[,'npp'],xlab = 'mean HYP rel. ab.', ylab = 'NPP',cex=cc, main = 'Mean HYP by rel. abundances',ylim = c(0,3000),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  abline(fit_npp_den,lwd=ll)
  tx <- paste('R2 =',round(summary(fit_npp_den)$r.squared,digits = 2))
  text(2.7,2600,tx,cex=cc)
  dev.off()
  
  fit_npp_mass_den <- lm(npp ~ mean_hyp_mass_den,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_npp_mass_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_mass_den'],data_sites[,'npp'],xlab = 'mean HYP biomass', ylab = 'NPP',cex=cc, main = 'Mean HYP by biomass',ylim = c(0,3000),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  abline(fit_npp_mass_den,lwd=ll)
  tx <- paste('R2 =',round(summary(fit_npp_mass_den)$r.squared,digits = 2))
  text(2.7,2600,tx,cex=cc)
  dev.off()
  
  fit_npp_ing_den <- lm(npp ~ mean_hyp_ing_den,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_npp_ing_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_ing_den'],data_sites[,'npp'],xlab = 'mean HYP en. int.', ylab = 'NPP',cex=cc, main = 'Mean HYP by energy intake',ylim = c(0,3000),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_npp_ing_den)$r.squared,digits = 2))
  abline(fit_npp_ing_den,lwd=ll)
  text(2.7,2600,tx,cex=cc)
  dev.off()
  
  res_npp <-c(round(summary(fit_npp)$r.squared,digits = 2),round(summary(fit_npp_den)$r.squared,digits = 2),round(summary(fit_npp_mass_den)$r.squared,digits = 2),round(summary(fit_npp_ing_den)$r.squared,digits = 2))
  
  #PREC
  
  fit_prec <- lm(prec ~ mean_hyp ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_prec_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp'],data_sites[,'prec'],xlab = 'mean HYP sp.', ylab = 'PREC',cex=cc, main = 'Mean HYP by species',ylim = c(0,3500),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_prec)$r.squared,digits = 2))
  abline(fit_prec,lwd=ll)
  text(2.7,3000,tx,cex=cc)
  dev.off()
  
  fit_prec_log_mass <- lm(prec ~ mean_log_mass, data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_log_mass_prec_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_log_mass'],data_sites[,'prec'],xlab = 'mean log MASS', ylab = 'PREC',cex=cc, main = 'Mean log MASS by species',ylim = c(0,3500),xlim = c(0,800),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_prec)$r.squared,digits = 2))
  tx <- paste('R2 =',round(summary(fit_prec_log_mass)$r.squared,digits = 2))
  abline(fit_prec_log_mass,lwd=ll)
  text(2.7,3000,tx,cex=cc)
  dev.off()
  
  fit_prec_den <- lm(prec ~ mean_hyp_den ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_prec_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_den'],data_sites[,'prec'],xlab = 'mean HYP rel. ab.', ylab = 'PREC',cex=cc, main = 'Mean HYP by rel. abundances',ylim = c(0,3500),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_prec_den)$r.squared,digits = 2))
  text(2.7,3000,tx,cex=cc)
  abline(fit_prec_den,lwd=ll)
  dev.off()
  
  fit_prec_mass_den <- lm(prec ~ mean_hyp_mass_den,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_prec_mass_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_mass_den'],data_sites[,'prec'],xlab = 'mean HYP biomass', ylab = 'PREC',cex=cc, main = 'Mean HYP by biomass',ylim = c(0,3500),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_prec_mass_den)$r.squared,digits = 2))
  text(2.7,3000,tx,cex=cc)
  abline(fit_prec_mass_den,lwd=ll)
  dev.off()
  
  fit_prec_ing_den <- lm(prec ~ mean_hyp_ing_den,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_prec_ing_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_ing_den'],data_sites[,'prec'],xlab = 'mean HYP en. int.', ylab = 'PREC',cex=cc, main = 'Mean HYP by energy intake',ylim = c(0,3500),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_prec_ing_den)$r.squared,digits = 2))
  text(2.7,3000,tx,cex=cc)
  abline(fit_prec_ing_den,lwd=ll)
  dev.off()
  
  res_prec <-c(round(summary(fit_prec)$r.squared,digits = 2),round(summary(fit_prec_den)$r.squared,digits = 2),round(summary(fit_prec_mass_den)$r.squared,digits = 2),round(summary(fit_prec_ing_den)$r.squared,digits = 2))
  
  #TEMP
  
  fit_temp <- lm(temp ~ mean_hyp ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_temp_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp'],data_sites[,'temp'],xlab = 'mean HYP sp.', ylab = 'TEMP',cex=cc, main = 'Mean HYP by species',ylim = c(-15,35),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_temp)$r.squared,digits = 2))
  abline(fit_temp,lwd=ll)
  text(1.5,-12,tx,cex=cc)
  dev.off()
  
  fit_temp_log_mass <- lm(temp ~ mean_log_mass, data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_log_mass_temp_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_log_mass'],data_sites[,'temp'],xlab = 'mean log MASS', ylab = 'TEMP',cex=cc, main = 'Mean log MASS by species',ylim = c(-15,35),xlim = c(0,800),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_temp_log_mass)$r.squared,digits = 2))
  abline(fit_temp_log_mass,lwd=ll)
  text(600,-12,tx,cex=cc)
  dev.off()
  
  fit_temp_den <- lm(temp ~ mean_hyp_den ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_temp_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_den'],data_sites[,'temp'],xlab = 'mean HYP rel. ab.', ylab = 'TEMP',cex=cc, main = 'Mean HYP by rel. abundances',ylim = c(-15,35),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_temp_den)$r.squared,digits = 2))
  abline(fit_temp_den,lwd=ll)
  text(1.5,-12,tx,cex=cc)
  dev.off()
  
  fit_temp_mass_den <- lm(temp ~ mean_hyp_mass_den ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_temp_mass_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_mass_den'],data_sites[,'temp'],xlab = 'mean HYP biomass', ylab = 'TEMP',cex=cc, main = 'Mean HYP by biomass',ylim = c(-15,35),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_temp_mass_den)$r.squared,digits = 2))
  abline(fit_temp_mass_den,lwd=ll)
  text(1.5,-12,tx,cex=cc)
  dev.off()
  
  fit_temp_ing_den <- lm(temp ~ mean_hyp_ing_den ,data = data_sites)
  
  file_name <- paste('plots/',dd_now,'/fig_hyp_temp_ing_den_',dd_now,'.pdf',sep = '')
  pdf(file_name,width = wd1, height = ht1)
  plot(data_sites[,'mean_hyp_ing_den'],data_sites[,'temp'],xlab = 'mean HYP en. int.', ylab = 'TEMP',cex=cc, main = 'Mean HYP by energy intake',ylim = c(-15,35),xlim = c(1,3),col = coltrans,cex.axis = cc,cex.lab = cc,lwd = llci)
  tx <- paste('R2 =',round(summary(fit_temp_ing_den)$r.squared,digits = 2))
  abline(fit_temp_ing_den,lwd=ll)
  text(1.5,-12,tx,cex=cc)
  dev.off()
  
  res_temp <-c(round(summary(fit_temp)$r.squared,digits = 2),round(summary(fit_temp_den)$r.squared,digits = 2),round(summary(fit_temp_mass_den)$r.squared,digits = 2),round(summary(fit_temp_ing_den)$r.squared,digits = 2))
  
  res_sites <- rbind(res_sites,c(dd_now,dim(data_sites)[1],'NPP',res_npp,'PREC',res_prec,'TEMP',res_temp))
}

write.table(res_sites, file = "outputs/res_sites.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
