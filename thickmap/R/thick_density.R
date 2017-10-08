library(ks)            # for weighted kernel density
library(pracma)        # for trapz formula
library(matrixStats)   # for colSds function (column std)
library(pryr)          # for evaluating memory usage in R

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())  # remove whole workspace

start.time <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)

datasets = list()
datasets_numbers = list()
densities = list()           ## initialize list
thick_data_mean = list()     ## 
thick_data_std = list()      ##
thick_data_up = list()       ##
thick_data_dn = list()       ##
#---------------------------------------------------------------------
# (0) Constants
#---------------------------------------------------------------------
#thickdata_info = '/home/lovric_g/application-data/eclipse-workspace/lungs-3D-quant/thickmap/R/data_info_10x_optics.txt' ## >> when using 10x optics
thickdata_info = '/home/lovric_g/application-data/eclipse-workspace/lungs-3D-quant/thickmap/R/data_info_2-4x_optics.txt' ## >> when using 2-4x optics
#data_xyz = c(1079,1201,449)    ## Data size
data_xyz = c(1479,1565,729)    ## Data size
#density_bw = 2                 ## Density bandwidth >> when using 10x optics
density_bw = 4                 ## >> when using 2-4x optics
#px_size = 1.1E-6             ## Pixel size (10x optics)
px_size = 2.9E-6             ## Pixel size (2-4x optics)

#---------------------------------------------------------------------
# (1) Load datasets from "data_info" file
#---------------------------------------------------------------------
conn <- file(thickdata_info,open="r")
linn <-readLines(conn)

k = 1

## Find how many datasets and note Linenumbers
for (i in 1:length(linn)){
  #print(linn[i])
  if (substr(linn[i], 1, 10) == '## dataset') {
    datasets_numbers[[k]] = i
    k = k+1
  }
  if (i == length(linn)) {
    datasets_numbers[[k]] = length(linn)+2
  }
}

## Remove non-existing datasets
tmp_length = length(datasets_numbers)-1
tmp_datasets_numbers = datasets_numbers
for ( ii in 1:tmp_length ) {
  upper = tmp_datasets_numbers[[tmp_length-ii+2]]
  lower = tmp_datasets_numbers[[tmp_length-ii+1]]
  if ( (upper - lower) < 3 ) {
    ind = - (tmp_length-ii+2)
    datasets_numbers = datasets_numbers[ind]
    }
  }

#stop("stopped by GORRRAAAANNN")
## Loop through each Dataset
for ( ii in 1:(length(datasets_numbers)-1) ) {
  n_datasets = datasets_numbers[[ii+1]] - datasets_numbers[[ii]] -2
  start_i = datasets_numbers[[ii]]+1
  end_i   = datasets_numbers[[ii+1]]-2
  
  ## Loop through multiple Datasets_per_volume for each Dataset
  if ( n_datasets == 1 ) {
    datasets_per_vol = linn[start_i]
  } else {
    datasets_per_vol = c(linn[start_i],linn[start_i+1])
  }
  
  if ( n_datasets >= 2 ) {
    for (kk in (start_i+2):end_i) {
      datasets_per_vol = c(datasets_per_vol,linn[kk])
    }
  }
  datasets[[ii]] = datasets_per_vol
}

close(conn)

#---------------------------------------------------------------------
# (1) Input data and scale to px size
#---------------------------------------------------------------------

## Loop for Datasets under different pressures
for ( ll in 1:(length(datasets_numbers)-1) ) {

  datasets_per_vol = datasets[[ll]]
  
  ## Loop for Datasets of same pressure but with different segmentation algorithms  
  for (ii in 1:length(datasets_per_vol) ) {
    print(sprintf("--------- starting loop ------------"))
    
    print(datasets_per_vol[ii])
    myfile <- file(datasets_per_vol[ii],'rb')
    thickness_dat <-readBin(myfile,double(),data_xyz[1]*data_xyz[2]*data_xyz[3],size=4,endian='big')
    thickness_dat <- thickness_dat[thickness_dat != 0]  # remove zeros (masking)
    print('------ file read ------')
    
    thickness_dat_mu <- thickness_dat * px_size / 1E-6

    close(myfile)
    rm(thickness_dat)

    #---------------------------------------------------------------------
    # (2) Density calculation
    #---------------------------------------------------------------------

    print('------ starting KDE ------')
    
    # Check whether temporary file exists - if YES, use it - if NO, recalculate
    
    tmp_string = strsplit(datasets_per_vol[ii], '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)
    temp_file = tmp_string[[1]][length(tmp_string[[1]])]
    temp_file = substr(temp_file,1,nchar(temp_file)-4)
    temp_file = paste(temp_file,'temp.bin',sep = "_")
    orig_path = paste(tmp_string[[1]][1:(length(tmp_string[[1]])-1)], collapse = "/")
    temp_file = paste(orig_path,temp_file,sep = "/")
    
    if ( file.exists(temp_file) ) {
      load(temp_file)
      } else {
        map1 <- kde(thickness_dat_mu,
                    h=density_bw,
                    gridsize = 512,
                    xmin=0,
                    xmax=170,
                    binned=TRUE,
                    bgridsize = 512 )
        #save(map1, file=temp_file, compress=FALSE)
    }
    
    print('------ finished KDE // saving results // and restarting Loop ------')
    # density <- list("x" = map1$eval.points, "y" = map1$estimate)
    # densities[[ii]] <- density
    densities[[ii]] <- list("x" = map1$eval.points, "y" = map1$estimate)
    rm(thickness_dat_mu)
    rm(map1)

  } ### >> FOR-LOOP datasets_per_vol


  #---------------------------------------------------------------------
  # (3) Calculate confidence intervals
  #---------------------------------------------------------------------
  print(sprintf("--------- calculating confidence ------------"))

  thick_data_mat <- matrix(nrow = length(datasets_per_vol), ncol = size(densities[[1]]$x)[2], byrow = FALSE,
                           dimnames = NULL)
  for (ii in 1:length(datasets_per_vol) ) {
    thick_data_mat[ii,] = densities[[ii]]$y
  }

  thick_data_mean[[ll]] <- colMeans(thick_data_mat)
  thick_data_std[[ll]]  <- colSds(thick_data_mat)

  thick_data_up[[ll]] = thick_data_mean[[ll]]+thick_data_std[[ll]]
  thick_data_dn[[ll]] = thick_data_mean[[ll]]-thick_data_std[[ll]]
  
  
  
  #---------------------------------------------------------------------
  # (4) Calculate integrals
  #---------------------------------------------------------------------
  
  region_1 <- densities[[1]]$x >= 20 & densities[[1]]$x < 50
  region_2 <- densities[[1]]$x >= 50 & densities[[1]]$x < 80
  region_3 <- densities[[1]]$x >= 80 & densities[[1]]$x < 110
  region_4 <- densities[[1]]$x >= 110
  
  mean_data = thick_data_mean[[ll]]
  up_data   = thick_data_up[[ll]]
  dn_data   = thick_data_dn[[ll]]
  
  integral_mean = trapz(densities[[1]]$x[region_1],mean_data[region_1])
  print(sprintf("Integral 1 (mean): %4.3f", integral_mean*100))
  integral_up = trapz(densities[[1]]$x[region_1],up_data[region_1])
  print(sprintf("Integral 1 (+/-): %4.3f", abs(integral_up - integral_mean)*100 ))
  
  integral_mean = trapz(densities[[1]]$x[region_2],mean_data[region_2])
  print(sprintf("Integral 2 (mean): %4.3f", integral_mean*100))
  integral_up = trapz(densities[[1]]$x[region_2],up_data[region_2])
  print(sprintf("Integral 2 (+/-): %4.3f", abs(integral_up - integral_mean)*100))
  
  integral_mean = trapz(densities[[1]]$x[region_3],mean_data[region_3])
  print(sprintf("Integral 3 (mean): %4.3f", integral_mean*100))
  integral_up = trapz(densities[[1]]$x[region_3],up_data[region_3])
  print(sprintf("Integral 3 (+/-): %4.3f", abs(integral_up - integral_mean)*100))
  
  integral_mean = trapz(densities[[1]]$x[region_4],mean_data[region_4])
  print(sprintf("Integral 4 (mean): %4.3f", integral_mean*100))
  integral_up = trapz(densities[[1]]$x[region_4],up_data[region_4])
  print(sprintf("Integral 4 (up): %4.3f", abs(integral_up - integral_mean)*100))

} ### >> FOR-LOOP datasets

#---------------------------------------------------------------------
# (5) Plots
#---------------------------------------------------------------------
output = "/home/lovric_g/application-data/eclipse-workspace/lungs-quantitative/img/thickness_density_2-4x.eps"
#output = "~/Desktop/thickness_density_10x.eps"
output = "~/Desktop/thickness_density_2-4x.eps"
#output = "~/Desktop/image.eps"
cairo_ps(output,height=5, width=6)

# ALTERNATIVE METHOD
#setEPS()
#postscript("~/Desktop/image.eps")

print(sprintf("--------- plotting ------------"))

# Density
par(mfrow=c(1,1))
plot(densities[[1]]$x,
     thick_data_mean[[1]],
     col=4,
     xlab=expression(paste("Structure diameter [", mu, "m]")),
     ylab="Propability density distribution (PDF)",
     type='l',
     xlim=c(0,170),
     ylim=c(0,0.020))
polygon(c(densities[[1]]$x,rev(densities[[1]]$x)),c(thick_data_dn[[1]],rev(thick_data_up[[1]])),col=rgb(0, 0, 1,0.25), border = FALSE)

if (length(datasets) > 1) {
  lines(densities[[1]]$x,thick_data_mean[[2]],col=rgb(0, 1,0))
  polygon(c(densities[[1]]$x,rev(densities[[1]]$x)),c(thick_data_dn[[2]],rev(thick_data_up[[2]])),col=rgb(0, 1, 0,0.25), border = FALSE) 
}

if (length(datasets) > 2) {
  lines(densities[[1]]$x,thick_data_mean[[3]],col=rgb(1, 0,0))
  polygon(c(densities[[1]]$x,rev(densities[[1]]$x)),c(thick_data_dn[[3]],rev(thick_data_up[[3]])),col=rgb(1, 0, 0,0.25), border = FALSE)
}

# for (ii in 2:length(datasets) ) {
#   lines(densities[[ii]]$x,densities[[ii]]$y,col=ii)
# }
#lines(densities[[ii]]$x,thick_data_mean,col=10)
legend(115, 0.02, c("10 cmH2O","20 cmH2O", "30 cmH2O"), col = c(4,3,2),
       text.col = "black", lty = c(1, 1, 1),# pch = c(NA, NA, 1,1,1),
       merge = TRUE, bg = "white")

dev.off()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)