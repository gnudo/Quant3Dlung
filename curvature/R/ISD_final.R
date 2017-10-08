library(ks)            # for weighted kernel density
library(pracma)        # for trapz formula
library(matrixStats)   # for colSds function (column std)
library(fields)        # for image.plot

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())  # remove whole workspace

source('R_functions.R')

start.time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)

# densitiesISD       = list()     ## initialize list
# ISD_data_mean  = list()     ##
# ISD_data_std   = list()     ##
# ISD_data_up    = list()     ##
# ISD_data_dn    = list()     ##
#---------------------------------------------------------------------
# (0) Constants
#---------------------------------------------------------------------
datasrc = '/home/lovric_g/DATA_HD/quant-paper/curvature/'
curvdata_info = '/home/lovric_g/application-data/eclipse-workspace/lungs-3D-quant/curvature/R/curvdatainfo.txt'

px_size = 1.1E-6             ## Pixel size (10x optics)

colors = c('blue','green','red')
nbin = 4096
wished_range = 320
# n_reduced = 1e5

density_bw_Mean  = 1e-5               ## Density bandwidth
density_bw_Gauss = 1e-8               ## Density bandwidth
# density_bw_Mean  = 1e-5               ## Density bandwidth for undersampled datasets
# density_bw_Gauss = 1e-5               ## Density bandwidth for undersampled datasets

#---------------------------------------------------------------------
# (1) Load datasets from "data_info" file
#---------------------------------------------------------------------
returnlist = loadDatasetsFromTextfile(curvdata_info)
datasets = returnlist[[1]]
datasets_numbers = returnlist[[2]]

#---------------------------------------------------------------------
# (2) Input data
#---------------------------------------------------------------------

## Loop for Datasets under different pressures
for ( ll in 1:(length(datasets_numbers)-1) ) {
  
  ## Compose here a different datasets_per_volume than in thickmap
  datasets_per_vol = c( datasets[[ll]][[1]] )
  for ( ind in 2:(length(datasets[[ll]])) ) {
    datasets_per_vol = c(datasets_per_vol,datasets[[ll]][[ind]])
  }

  
  #---------------------------------------------------------------------
  # (3) Input data with SAME pressure, but different segmentations
  #---------------------------------------------------------------------
  datasets_tmp = composeDatasets(args,datasrc,datasets_per_vol)
  
  for (ii in 1:length(datasets_tmp$csv_file) ) {
    #foreach (ii =1:length(csv_file) ) %dopar% {
    print(sprintf("--------- starting loop ------------"))
    print(datasets_per_vol[ii])
    X_sample = loadFile(datasets_tmp$csv_file[ii],datasets_tmp$csv_file_red[ii],datasets_tmp$csv_file_bin[ii])
    
    # kappa_1 <- X_sample[,1] * px_size / 1E-6
    # kappa_2 <- X_sample[,2] * px_size / 1E-6
    kappa_1 <- -X_sample[,2] * px_size / 1E-6  ## IF THE RAW DATA WAS INVERTED
    kappa_2 <- -X_sample[,1] * px_size / 1E-6  ## IF THE RAW DATA WAS INVERTED
    
    areas_tmp <- X_sample[,3]
    classifier <- X_sample[,4]
    
    rm(X_sample)
    

    # n_vertices = length(kappa_1)
    # A_total = sum(areas)

    ## Filter high-frequency noise (ultra-high curvature values are not realistic)
    regions <- kappa_2 <= 0.6 & kappa_1 >= -0.6
    areas = areas_tmp[regions]
    x <- cbind(kappa_1[regions],kappa_2[regions])
    
    # x <- cbind(kappa_1,kappa_2)

    rm(areas_tmp)
    rm(regions)
    rm(kappa_1)
    rm(kappa_2)
    #---------------------------------------------------------------------
    # (4) ISD Plot
    #---------------------------------------------------------------------
    map2 <- kde(x,
                gridsize=nbin,  ## changing has no effect
                #H=kernel,
                xmin=c(-0.6,-0.6),
                xmax=c(0.6,0.6),   ## IF this interval too big, integral doesn't yield 100%
                # xmin=c(-0.3,-0.3),
                # xmax=c(0.3,0.3),
                binned=TRUE,
                bgridsize=c(nbin,nbin),
                w=areas*length(areas)/sum(areas)
                # w=areas/sum(areas)
    )

    int2d <- trapz2D( map2$estimate,map2$eval.points[[1]],map2$eval.points[[2]],-0.6,0.6,-0.6,0.6 )
    print(sprintf("Integral (KDE): %4.3f", int2d*100))
    
    ## Integrate testing if we cut all values higher than wished-range
    # higher_cut = map2$estimate
    # higher_cut[map2$estimate > wished_range] = 0
    # int2d_alt <- trapz2D( higher_cut,x1,x2,-0.6,0.6,-0.6,0.6 )
    # print(sprintf("Integral when higher values omitted(KDE): %4.3f", int2d_alt*100))  # This tests whether bandwidth/kernel is chosen well
    # rm(higher_cut)
    
    # map2 <- kde(gauss_curv,
    #             gridsize = nbin,
    #             #h=density_bw_Gauss,
    #             xmin=-0.05,
    #             xmax=0.05,
    #             binned=TRUE,
    #             bgridsize = nbin,
    #             w=areas*length(areas)/sum(areas)
    # )
    print('------ finished ISD-KDE // saving results // and restarting Loop ------')
    # densitiesISD[[ii]] <- list("x" = map2$eval.points, "y" = map2$estimate)
    
    if (ii == 1) {
      mean_ISD = map2$estimate / length(datasets_tmp$csv_file)
    } else {
      mean_ISD = mean_ISD + ( map2$estimate / length(datasets_tmp$csv_file) )
    }
    
    x1 = map2$eval.points[[1]]
    x2 = map2$eval.points[[2]]
    
    rm(map2)
      
  } ### >> FOR-LOOP datasets
  
  # if (ll == 1) {
  #   mean1 = mean_ISD
  # }  else if (ll == 2) {
  #   mean2 = mean_ISD
  #   plotDensityMat(x1,x2,abs(mean2-mean1),100)
  # } else if (ll == 3) {
  #   mean3 = mean_ISD
  #   plotDensityMat(x1,x2,abs(mean3-mean2),100)
  # }
  
  # plotDensityMat(map2$eval.points[[1]],map2$eval.points[[2]],map2$estimate,)
  plotDensityMat(x1,x2,mean_ISD,wished_range)
  
  
  
  
  
  #---------------------------------------------------------------------
  # (6) Calculate confidence intervals
  #---------------------------------------------------------------------
  # print(sprintf("--------- calculating confidence ------------"))
  # 
  # ## Mean
  # curvMean_data_mat <- matrix(nrow = length(datasets_per_vol), ncol = size(densitiesMean[[1]]$x)[2], byrow = FALSE,
  #                             dimnames = NULL)
  # for (ii in 1:length(datasets_per_vol) ) {
  #   curvMean_data_mat[ii,] = densitiesMean[[ii]]$y
  # }
  # 
  # curvMean_data_mean[[ll]] <- colMeans(curvMean_data_mat)
  # curvMean_data_std[[ll]]  <- colSds(curvMean_data_mat)
  # 
  # curvMean_data_up[[ll]] = curvMean_data_mean[[ll]]+curvMean_data_std[[ll]]
  # curvMean_data_dn[[ll]] = curvMean_data_mean[[ll]]-curvMean_data_std[[ll]]
  # 
  # ## Gauss
  # curvGauss_data_mat <- matrix(nrow = length(datasets_per_vol), ncol = size(densitiesGauss[[1]]$x)[2], byrow = FALSE,
  #                              dimnames = NULL)
  # for (ii in 1:length(datasets_per_vol) ) {
  #   curvGauss_data_mat[ii,] = densitiesGauss[[ii]]$y
  # }
  # 
  # curvGauss_data_mean[[ll]] <- colMeans(curvGauss_data_mat)
  # curvGauss_data_std[[ll]]  <- colSds(curvGauss_data_mat)
  # 
  # curvGauss_data_up[[ll]] = curvGauss_data_mean[[ll]]+curvGauss_data_std[[ll]]
  # curvGauss_data_dn[[ll]] = curvGauss_data_mean[[ll]]-curvGauss_data_std[[ll]]
  # 
  # 
  # rm(datasets_per_vol)
  
  #---------------------------------------------------------------------
  # (7) Calculate integrals
  #---------------------------------------------------------------------
  # ## Mean
  # Mean_mean_data = curvMean_data_mean[[ll]]
  # Mean_up_data   = curvMean_data_up[[ll]]
  # Mean_dn_data   = curvMean_data_dn[[ll]]
  # Mean_integral_mean    = trapz(densitiesMean[[1]]$x,Mean_mean_data)
  # print(sprintf("Integral (Mean curvature): %4.3f", Mean_integral_mean*100))
  # Mean_integral_up = trapz(densitiesMean[[1]]$x,Mean_up_data)
  # print(sprintf("Integral (Mean curvature) (+/-): %4.3f", abs(Mean_integral_up - Mean_integral_mean)*100 ))
  # 
  # ## Gauss
  # Gauss_mean_data = curvGauss_data_mean[[ll]]
  # Gauss_up_data   = curvGauss_data_up[[ll]]
  # Gauss_dn_data   = curvGauss_data_dn[[ll]]
  # Gauss_integral_mean    = trapz(densitiesGauss[[1]]$x,Gauss_mean_data)
  # print(sprintf("Integral (Gauss curvature): %4.3f", Gauss_integral_mean*100))
  # Gauss_integral_up = trapz(densitiesGauss[[1]]$x,Gauss_up_data)
  # print(sprintf("Integral (Gauss curvature) (+/-): %4.3f", abs(Gauss_integral_up - Gauss_integral_mean)*100 ))
  
} ### >> FOR-LOOP datasets




#---------------------------------------------------------------------
# (7) Plots
#---------------------------------------------------------------------
# print(sprintf("--------- plotting ------------"))
# 
# # 1 - Mean curvature density
# output = "~/Desktop/mean_curv.eps"
# #cairo_ps(output,height=6, width=6)
# 
# par(mfrow=c(1,1))
# plot(densitiesMean[[1]]$x,
#      curvMean_data_mean[[1]],
#      col=4,
#      xlab=expression(paste("Mean curvature [", mu, "m"^-1,"]")),
#      ylab="Propability density distribution (PDF)",
#      type='l',
#      xlim=c(-0.2,0.1),
#      ylim=c(0,15))
# polygon(c(densitiesMean[[1]]$x,rev(densitiesMean[[1]]$x)),
#         c(curvMean_data_dn[[1]],rev(curvMean_data_up[[1]])),
#         col=rgb(0, 0, 1,0.25),
#         border = FALSE)
# 
# if (length(datasets) > 1) {
#   lines(densitiesMean[[1]]$x,curvMean_data_mean[[2]],col=rgb(0, 1,0))
#   polygon(c(densitiesMean[[1]]$x,rev(densitiesMean[[1]]$x)),
#           c(curvMean_data_dn[[2]],rev(curvMean_data_up[[2]])),
#           col=rgb(0, 1, 0,0.25),
#           border = FALSE) 
# }
# 
# if (length(datasets) > 2) {
#   lines(densitiesMean[[1]]$x,curvMean_data_mean[[3]],col=rgb(1, 0,0))
#   polygon(c(densitiesMean[[1]]$x,rev(densitiesMean[[1]]$x)),
#           c(curvMean_data_dn[[3]],rev(curvMean_data_up[[3]])),
#           col=rgb(1, 0, 0,0.25),
#           border = FALSE)
# }
# 
# legend(-0.15, 10, c(expression("10 cmH"[2]*"O"),expression("20 cmH"[2]*"O"), expression("30 cmH"[2]*"O")), col = c(4,3,2),
#        text.col = "black", lty = c(1, 1, 1),# pch = c(NA, NA, 1,1,1),
#        merge = TRUE, bg = "white")
# 
# # dev.off()
# 
# 
# # 2 - Gauss curvature density
# output = "~/Desktop/gauss_curv.eps"
# #cairo_ps(output,height=6, width=6)
# 
# par(mfrow=c(1,1))
# plot(densitiesGauss[[1]]$x,
#      curvGauss_data_mean[[1]],
#      col=4,
#      xlab=expression(paste("Gauss curvature [", mu, "m"^-2,"]")),
#      ylab="Propability density distribution (PDF)",
#      type='l',
#      xlim=c(-0.01,0.01),
#      ylim=c(0,1200))
# polygon(c(densitiesGauss[[1]]$x,rev(densitiesGauss[[1]]$x)),
#         c(curvGauss_data_dn[[1]],rev(curvGauss_data_up[[1]])),
#         col=rgb(0, 0, 1,0.25),
#         border = FALSE)
# 
# if (length(datasets) > 1) {
#   lines(densitiesGauss[[1]]$x,curvGauss_data_mean[[2]],col=rgb(0, 1,0))
#   polygon(c(densitiesGauss[[1]]$x,rev(densitiesGauss[[1]]$x)),
#           c(curvGauss_data_dn[[2]],rev(curvGauss_data_up[[2]])),
#           col=rgb(0, 1, 0,0.25),
#           border = FALSE)
# }
# 
# if (length(datasets) > 2) {
#   lines(densitiesGauss[[1]]$x,curvGauss_data_mean[[3]],col=rgb(1, 0,0))
#   polygon(c(densitiesGauss[[1]]$x,rev(densitiesGauss[[1]]$x)),
#           c(curvGauss_data_dn[[3]],rev(curvGauss_data_up[[3]])),
#           col=rgb(1, 0, 0,0.25),
#           border = FALSE)
# }
# 
# legend(0.005, 1200, c(expression("10 cmH"[2]*"O"),expression("20 cmH"[2]*"O"), expression("30 cmH"[2]*"O")), col = c(4,3,2),
#        text.col = "black", lty = c(1, 1, 1),# pch = c(NA, NA, 1,1,1),
#        merge = TRUE, bg = "white")

#dev.off()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)