library(ks)            # for weighted kernel density
library(pracma)        # for trapz formula
library(matrixStats)   # for colSds function (column std)
library(pryr)          # for evaluating memory usage in R

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())  # remove whole workspace

start.time <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)

#---------------------------------------------------------------------
# (0) Constants
#---------------------------------------------------------------------
datasets = c(
  '/home/lovric_g/DATA_HD/quant-paper/2x-4x/mouseB01crop_seg1_LocThk.raw',
  '/home/lovric_g/DATA_HD/quant-paper/2x-4x/mouseB01crop_seg2_LocThk.raw',
  '/home/lovric_g/DATA_HD/quant-paper/2x-4x/mouseB01crop_seg3_LocThk.raw'
  #'PSI-M2E-exv25sc_faces'
)
data_xyz = c(1479,1565,729)    ## Data size
density_bw = 2               ## Density bandwidth
#px_size = 1.1E-6             ## Pixel size (10x optics)
px_size = 2.9E-6             ## Pixel size (2-4x optics)

#---------------------------------------------------------------------
# (1) Load datasets
#---------------------------------------------------------------------
print(sprintf("--------- starting loop ------------"))

  for ( ii in 1:length(datasets) ) {
    print(datasets[ii])
    myfile <- file(datasets[ii],'rb')
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
    
    tmp_string = strsplit(datasets[ii], '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)
    temp_file = tmp_string[[1]][length(tmp_string[[1]])]
    temp_file = substr(temp_file,1,nchar(temp_file)-4)
    temp_file = paste(temp_file,'temp.bin',sep = "_")
        
    map1 <- kde(thickness_dat_mu,
                h=density_bw,
                gridsize = 512,
                xmin=0,
                xmax=170,
                binned=TRUE,
                bgridsize = 512 )
    save(map1, file=temp_file, compress=FALSE)
    
    print('------ finished KDE // saving results // and restarting Loop ------')
    # density <- list("x" = map1$eval.points, "y" = map1$estimate)
    # densities[[ii]] <- density
    #densities[[ii]] <- list("x" = map1$eval.points, "y" = map1$estimate)
    rm(thickness_dat_mu)
    rm(map1)

  } ### >> FOR-LOOP datasets_per_vol

#---------------------------------------------------------------------
#---------------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)