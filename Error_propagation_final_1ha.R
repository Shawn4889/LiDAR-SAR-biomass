#Author: Xiaoxuan Li
#02082023
library(corrplot)
library(tidyverse)
library(psych)
library(xlsx)
library(car)
library(MASS)
library(rsq)
library(ggpmisc)
library(ggrepel)
library(dplyr)
library(rlist)
library(data.table)
library(lemon)
library(lmtest)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(mltools)
library(tidyr)
library(stats)
library(lidR)
library(tools) 
library(raster)
library(ggpointdensity)
library(viridis)
library(grid)
library(readxl)
library(ehaGoF)
library(Metrics)
library(rgdal)
library(caret)
library(randomForest)
library(DAAG)
library(matrixStats)
library(FactoMineR)
library(ggcorrplot)
library(devtools)
library(factoextra)
library(purrr)
library(ggridges)
library(terra)


filedir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\CSIR_AGBD_1ha.csv"
Data = read.csv(filedir)

dir_hcc_gedi <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_gedi_1ha.csv"
data_hcc_gedi = read.csv(dir_hcc_gedi)
data_hcc_gedi$rh50 <- data_hcc_gedi$rh50 + 100
df_HCC <- data.frame(data_hcc_gedi[,"HCC"])
colnames(df_HCC) <- c("HCC")

field_sim<-list()
Data_final<-list()
lidar_model<-list()
lidar_est<-list()
CL95<-list()

error_field_func <- function(){
  field_AGB <- Data[,3]
  error_measurement <- 0.1
  error_allometry <- 0.22*0.25
  error_sampling <- 0.1
  error_field <- sqrt((error_measurement)^2+(error_allometry)^2+(error_sampling)^2)
  MC <- rnorm(length(field_AGB), mean = 0, sd = 1)
  field_AGB_sim <- field_AGB*(1+MC*error_field)
  return(field_AGB_sim)
}

a<- list()
b<-list()
for(i in 1:100) {
  field_sim[[i]] <- error_field_func()
  Data_final <- data.frame(cbind(field_sim[[i]],Data[,2]))
  colnames(Data_final) <- c("AGBD","HCC")
  lidar_model[[i]] <- lm(formula = AGBD ~ HCC, Data_final)
  lidar_est[[i]] <- predict(lidar_model[[i]], newdata = df_HCC)
  a[[i]] <- lidar_model[[i]]$coefficients[[1]]
  b[[i]] <- lidar_model[[i]]$coefficients[[2]]
}

lidar_est_df <- data.frame(lidar_est)
lidar_est_df_t <- t(lidar_est_df)
#lidar agbd samples in total
for(i in 1:nrow(lidar_est_df)){
  CL95[[i]] <- (quantile(lidar_est_df_t[,i], 0.975)[[1]]- 
                  quantile(lidar_est_df_t[,i], 0.025)[[1]])/2
}

# A list of prediction uncertainty for each LiDAR subplot, 175 in total.
pred_error <- (data.frame(unlist(CL95))/data.frame(colMeans((lidar_est_df_t))))*100
colnames(pred_error) <- c("Uncertainty")
Data_sel_2 <- data.frame(pred_error[pred_error$Uncertainty > 0 & pred_error$Uncertainty < 50,])
list <- c(row.names(Data_sel_2))
lidar_est_df_good <- data.frame(cbind(data_hcc_gedi[list, ],
                                      data.frame(lidar_est)[list, ]))

#free some memory
rm(field_sim)
rm(lidar_est)
rm(lidar_est_df)
rm(lidar_est_df_t)
rm(lidar_model)
rm(df_HCC)
rm(Data_sel_2)
rm(pred_error)
rm(CL95)


for(i in 1:100){
  colnames(lidar_est_df_good)[i+12] <- paste0("ALS_AGBD_",i)
}

cor(lidar_est_df_good$agbd,lidar_est_df_good$ALS_AGBD_5)

#lidar agbd to GEDI agbd
#------------------------------------------------------------------------------------------------#
#rf
gedi_model<-list()
gedi_est <-list()
CL95_gedi<-list()
dir_gedi_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\gedi_sar_100m.csv"
data_gedi_sar = read.csv(dir_gedi_sar)
data_gedi_sar$rh50 <- data_gedi_sar$rh50 + 100
data_gedi_sar <- data_gedi_sar[data_gedi_sar$rh75 >0,]
data_gedi_sar_var_rf = data_gedi_sar[,c("s1","GEDI_cover","fhd","rh50","rh75","rh90","rh98")]

for(i in 1:100){
  ALS_AGBD_name <- paste0("ALS_AGBD_",i)
  formula <- as.formula(paste(ALS_AGBD_name, 
                              "~ s1+GEDI_cover+fhd+rh50+rh75+rh90+rh98"))
  repeat_cv <- trainControl(method='cv', number=5)
  gedi_model[[i]] <- train(
    formula,
    data=lidar_est_df_good, 
    method='rf',
    ntree=100,
    nodesize = 4,
    trControl=repeat_cv)
  gedi_est[[i]] <- predict(gedi_model[[i]], newdata = data_gedi_sar_var_rf)}

#hist(lidar_est_df_good$ALS_AGBD_1)
#file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\table_MC\\ALS.csv")
#write.csv(lidar_est_df_good$ALS_AGBD_1, file = file_path)


gedi_est_df <- data.frame(gedi_est)
gedi_est_df_t <- t(gedi_est_df)
gedi_rf_hist <- data.frame(colMeans((gedi_est_df_t)))
#lidar agbd samples in total
for(i in 1:nrow(gedi_est_df)){
  CL95_gedi[[i]] <- (quantile(gedi_est_df_t[,i], 0.975)[[1]]- 
                       quantile(gedi_est_df_t[,i], 0.025)[[1]])/2}

pred_error_gedi <- (data.frame(unlist(CL95_gedi))/data.frame(colMeans((gedi_est_df_t))))*100
colnames(pred_error_gedi) <- c("Uncertainty")
Data_sel_3 <- data.frame(pred_error_gedi[pred_error_gedi$Uncertainty >= 0 & pred_error_gedi$Uncertainty < 50,])
list_gedi <- c(row.names(Data_sel_3))
gedi_est_df_good <- data.frame(cbind(data_gedi_sar[list_gedi, ],
                                     data.frame(gedi_est)[list_gedi, ]))

#hist(gedi_est_df_good$c..1....13.9660530193603...2....9.35513310157466...3....29.1795572180158..)
#file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\table_MC\\GEDI_RF.csv")
#write.csv(gedi_est_df_good$c..1....13.9660530193603...2....9.35513310157466...3....29.1795572180158.., file = file_path)

rm(gedi_est)
rm(gedi_est_df)
rm(gedi_est_df_t)
rm(gedi_model)
rm(data_gedi_sar)
rm(Data_sel_3)
rm(pred_error_gedi)
rm(CL95_gedi)

#GEDI agbd to sar agbd snp 
chunck_size <- c(0,2000001,4000001,6000001)

chunck_size_1 <- c(2000000,4000000,6000000,7364178)

for(m in 1:length(chunck_size)){
  sar_model<-list()
  sar_est <-list()
  CL95_sar<-list()
  
  raster_df_1m <- 
    data.frame(
      data.frame(
        rasterToPoints(
          raster("E:\\ScanSAR\\ScanSAR\\Result_1004\\sar_knp\\SAR_100_all.tif"))[,3])[chunck_size[[m]]:chunck_size_1[[m]],])
  raster_xy_df_1m <- 
    data.frame(
      data.frame(
        rasterToPoints(
          raster("E:\\ScanSAR\\ScanSAR\\Result_1004\\sar_knp\\SAR_100_all.tif"))[,1:2])[chunck_size[[m]]:chunck_size_1[[m]],])
  
  colnames(raster_df_1m) <- c("SAR")
  
  for(i in 1:100){
    gedi_est_df_good <- data.frame(
      gedi_est_df_good[gedi_est_df_good[,i+11] > 0,])
    sar_model[[i]] <- lm(formula = log(gedi_est_df_good[,i+11]) ~ SAR, 
                         gedi_est_df_good)
    sar_est[[i]] <- predict(sar_model[[i]], newdata = raster_df_1m)}
  
  sar_est_df <- exp(data.frame(sar_est))
  sar_est_df_t <- t(sar_est_df)
  
  #lidar agbd samples in total
  for(i in 1:nrow(sar_est_df)){
    CL95_sar[[i]] <- (quantile(sar_est_df_t[,i], 0.975)[[1]]- 
                        quantile(sar_est_df_t[,i], 0.025)[[1]])/2}
  
  result_ras <- cbind(raster_xy_df_1m,data.frame(colMeans(sar_est_df_t)))
  file_path_ras <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_rf",paste0("SAR_AGBD_",m,".csv"))
  write.csv(result_ras, file = file_path_ras)
  
  #result_uncertainty <- cbind(raster_xy_df_1m,(data.frame(unlist(CL95_sar))/data.frame(colMeans(sar_est_df_t)))*100)
  result_uncertainty_abs <- cbind(raster_xy_df_1m,data.frame(unlist(CL95_sar)))
  
  file_path_uncertainty <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_rf",paste0("SAR_Uncertainty_",m,".csv"))
  #write.csv(result_uncertainty, file = file_path_uncertainty)
  write.csv(result_uncertainty_abs, file = file_path_uncertainty)
}


#------------------------------------------------------------------------------------------------#
#glm
gedi_model<-list()
gedi_est <-list()
CL95_gedi<-list()
dir_gedi_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\gedi_sar_100m.csv"
data_gedi_sar = read.csv(dir_gedi_sar)
data_gedi_sar$rh50 <- data_gedi_sar$rh50 + 100
data_gedi_sar_var_glm = data_gedi_sar[,c("rh50","rh98","fhd")]

for(i in 1:100){
  ALS_AGBD_name <- paste0("ALS_AGBD_",i)
  formula <- as.formula(paste(ALS_AGBD_name, "~ fhd + rh50 + rh98"))
  repeat_cv <- trainControl(method='cv', number=5)
  lidar_est_df_good <- data.frame(
    lidar_est_df_good[lidar_est_df_good[,i+12] > 0,])
  gedi_model[[i]] <- train(
    formula, 
    data=lidar_est_df_good, 
    method='glm', 
    family = Gamma(link = "identity"),
    start=c(-600,-2,6,2),
    trControl=repeat_cv)
  gedi_est[[i]] <- predict(gedi_model[[i]], newdata = data_gedi_sar_var_glm)}


gedi_est_df <- data.frame(gedi_est)
gedi_est_df_t <- t(gedi_est_df)

#lidar agbd samples in total
for(i in 1:nrow(gedi_est_df)){
  CL95_gedi[[i]] <- (quantile(gedi_est_df_t[,i], 0.975)[[1]]- 
                       quantile(gedi_est_df_t[,i], 0.025)[[1]])/2}

pred_error_gedi <- (data.frame(unlist(CL95_gedi))/data.frame(colMeans((gedi_est_df_t))))*100
colnames(pred_error_gedi) <- c("Uncertainty")
Data_sel_3 <- data.frame(pred_error_gedi[pred_error_gedi$Uncertainty >= 0 & pred_error_gedi$Uncertainty < 50,])
list_gedi <- c(row.names(Data_sel_3))
gedi_est_df_good <- data.frame(cbind(data_gedi_sar[list_gedi, ],
                                     data.frame(gedi_est)[list_gedi, ]))

rm(gedi_est)
rm(gedi_est_df)
rm(gedi_est_df_t)
rm(gedi_model)
rm(data_gedi_sar)
rm(Data_sel_3)
rm(pred_error_gedi)
rm(CL95_gedi)

#GEDI agbd to sar agbd snp 
chunck_size <- c(0,2000001,4000001,6000001)

chunck_size_1 <- c(2000000,4000000,6000000,7364178)

for(m in 1:length(chunck_size)){
  gc()
  sar_model<-list()
  sar_est <-list()
  CL95_sar<-list()
  raster_df_1m <- 0
  raster_xy_df_1m <- 0
  sar_est_df <- 0
  sar_est_df_t <- 0 
  result_ras<- 0
  result_uncertainty_abs <- 0 
  
  raster_df_1m <- 
    data.frame(
      data.frame(
        rasterToPoints(
          raster("E:\\ScanSAR\\ScanSAR\\Result_1004\\sar_knp\\SAR_100.tif"))[,3])[chunck_size[[m]]:chunck_size_1[[m]],])
  
  raster_xy_df_1m <- 
    data.frame(
      data.frame(
        rasterToPoints(
          raster("E:\\ScanSAR\\ScanSAR\\Result_1004\\sar_knp\\SAR_100.tif"))[,1:2])[chunck_size[[m]]:chunck_size_1[[m]],])
  
  colnames(raster_df_1m) <- c("SAR")
  
  for(i in 1:100){
    gedi_est_df_good <- data.frame(
      gedi_est_df_good[gedi_est_df_good[,i+11] > 0,])
    sar_model[[i]] <- lm(formula = log(gedi_est_df_good[,i+11]) ~ SAR, 
                         gedi_est_df_good)
    sar_est[[i]] <- predict(sar_model[[i]], newdata = raster_df_1m)}
  
  sar_est_df <- exp(data.frame(sar_est))
  sar_est_df_t <- t(sar_est_df)
  
  #lidar agbd samples in total
  for(i in 1:nrow(sar_est_df)){
    CL95_sar[[i]] <- (quantile(sar_est_df_t[,i], 0.975)[[1]]- 
                        quantile(sar_est_df_t[,i], 0.025)[[1]])/2}
  
  result_ras <- cbind(raster_xy_df_1m,data.frame(colMeans(sar_est_df_t)))
  file_path_ras <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_glm",paste0("SAR_AGBD_",m,".csv"))
  write.csv(result_ras, file = file_path_ras)
  
  #result_uncertainty <- cbind(raster_xy_df_1m,(data.frame(unlist(CL95_sar))/data.frame(colMeans(sar_est_df_t)))*100)
  result_uncertainty_abs <- cbind(raster_xy_df_1m,data.frame(unlist(CL95_sar)))
  
  file_path_uncertainty <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_glm",paste0("SAR_Uncertainty_",m,".csv"))
  #write.csv(result_uncertainty, file = file_path_uncertainty)
  write.csv(result_uncertainty_abs, file = file_path_uncertainty)
}


for(m in 1:16){
  dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_glm"
  #dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_GEDI_SAR_rf"

  filedir1 <- file.path(dir,paste0("SAR_AGBD_",m,".csv"))
  Data = read.csv(filedir1)
  colnames(Data) <- c("Index","X","Y","AGBD")
  a <- rasterFromXYZ(Data[,c("X","Y","AGBD")])
  output <- file.path(dir,paste0("SAR_AGBD_",m,".tif"))
  writeRaster(a, output, format = "GTiff",overwrite=TRUE )
  
  filedir2 <- file.path(dir,paste0("SAR_Uncertainty_",m,".csv"))
  Data = read.csv(filedir2)
  colnames(Data) <- c("Index","X","Y","Uncertainty")
  a <- rasterFromXYZ(Data[,c("X","Y","Uncertainty")])
  output <- file.path(dir,paste0("SAR_Uncertainty_",m,".tif"))
  writeRaster(a, output, format = "GTiff", overwrite=TRUE)
}
