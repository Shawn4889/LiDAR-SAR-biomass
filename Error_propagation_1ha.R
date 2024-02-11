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

dir_hcc_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_1ha.csv"
data_hcc_sar = read.csv(dir_hcc_sar)
df_HCC <-data.frame(data_hcc_sar[,5])
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
b<- list()
for(i in 1:100) {
  field_sim[[i]] <- error_field_func()
  Data_final <- data.frame(cbind(field_sim[[i]],Data[,2]))
  colnames(Data_final) <- c("AGBD","HCC")
  lidar_model[[i]] <- lm(formula = AGBD ~ HCC, Data_final)
  lidar_est[[i]] <- predict(lidar_model[[i]], newdata = df_HCC)
  a[[i]] <- as.numeric(lidar_model[[i]]$coefficients[[1]])
  b[[i]] <- as.numeric(lidar_model[[i]]$coefficients[[2]])
}
mean(as.numeric(a))
mean(as.numeric(b))

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
lidar_est_df_good <- data.frame(cbind(data_hcc_sar[list, ],
                                      data.frame(lidar_est)[list, ]))

#free some memory
rm(field_sim)
rm(lidar_est)
rm(lidar_est_df)
rm(lidar_est_df_t)
rm(lidar_model)
rm(data_hcc_sar)
rm(df_HCC)
rm(Data_sel_2)
rm(pred_error)
rm(CL95)

#lidar agbd to snp sar agbd

chunck_size <- c(0,2000001,4000001,6000001)

chunck_size_1 <- c(2000000,4000000,6000000,7364178)

for(m in 1:length(chunck_size)){
  gc()
  sar_model<-list()
  sar_est <-list()
  CL95_sar<-list()
  sar_est_df<-0
  sar_est_df_t<-0
  
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
    lidar_est_df_good <- data.frame(
      lidar_est_df_good[
        lidar_est_df_good[,i+5] > 0,])
    sar_model[[i]] <- lm(formula = log(lidar_est_df_good[,i+5]) ~ SAR, 
                         lidar_est_df_good)
    sar_est[[i]] <- predict(sar_model[[i]], newdata = raster_df_1m)
  }
  
  sar_est_df <- exp(data.frame(sar_est))
  sar_est_df_t <- t(sar_est_df)
  
  #lidar agbd samples in total
  for(i in 1:nrow(sar_est_df)){
    CL95_sar[[i]] <- (quantile(sar_est_df_t[,i], 0.975)[[1]]- 
                        quantile(sar_est_df_t[,i], 0.025)[[1]])/2
  }
  
  result_ras <- cbind(raster_xy_df_1m,data.frame(colMeans(sar_est_df_t)))
  file_path_ras <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_SAR",paste0("SAR_AGBD_",m,".csv"))
  write.csv(result_ras, file = file_path_ras)
  
  #result_uncertainty <- cbind(raster_xy_df_1m,(data.frame(unlist(CL95_sar))/data.frame(colMeans(sar_est_df_t)))*100)
  result_uncertainty_abs <- cbind(raster_xy_df_1m,data.frame(unlist(CL95_sar)))
  
  file_path_uncertainty <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_SAR",paste0("SAR_Uncertainty_",m,".csv"))
  write.csv(result_uncertainty_abs, file = file_path_uncertainty)
  
}

#Output results
for(m in 1:4){
  dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\Field_ALS_SAR"
  filedir1 <- file.path(dir,paste0("SAR_AGBD_",m,".csv"))
  Data = read.csv(filedir1)
  colnames(Data) <- c("Index","X","Y","AGBD")
  a <- rasterFromXYZ(Data[,c("X","Y","AGBD")])
  output <- file.path(dir,paste0("SAR_AGBD_",m,".tif"))
  writeRaster(a, output, format = "GTiff", )
  
  filedir2 <- file.path(dir,paste0("SAR_Uncertainty_",m,".csv"))
  Data = read.csv(filedir2)
  colnames(Data) <- c("Index","X","Y","Uncertainty")
  a <- rasterFromXYZ(Data[,c("X","Y","Uncertainty")])
  output <- file.path(dir,paste0("SAR_Uncertainty_",m,".tif"))
  writeRaster(a, output, format = "GTiff", )
}


