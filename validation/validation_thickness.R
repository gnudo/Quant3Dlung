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
thickdata_info = 'fullValidation_LocThk.raw'
data_xyz = c(512,512,512)    ## Data size
density_bw = 1
px_size = 1.1E-6             ## Pixel size (10x optics)

#---------------------------------------------------------------------
# (1) Load ground truth data
#---------------------------------------------------------------------
gt_data = read.csv("quant_data.txt", fill = TRUE, sep = "")
diameters = gt_data$R *2 * px_size / 1E-6
Vol = 4/3 * pi * (gt_data$R)^3
percent_Vol = Vol * gt_data$X.px.

FWHM = 2.3548
percent_Vol = percent_Vol/sum(percent_Vol) *(1 /(FWHM*density_bw*px_size/1E-6)  )

#---------------------------------------------------------------------
# (1) Input data and scale to px size
#---------------------------------------------------------------------
sprintf("--------- starting calc ------------")
myfile <- file(thickdata_info,'rb')
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

tmp_string = strsplit(thickdata_info, '/', fixed = FALSE, perl = FALSE, useBytes = FALSE)
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
}
print('------ finished KDE ------')

densities <- list("x" = map1$eval.points, "y" = map1$estimate)
rm(thickness_dat_mu)
rm(map1)

#---------------------------------------------------------------------
# (5) Plots
#---------------------------------------------------------------------
#output = "~/Desktop/image.eps"
#cairo_ps(output,height=5, width=6)

print(sprintf("--------- plotting ------------"))

# Density
par(mfrow=c(1,1))
plot(densities$x,
     densities$y,
     col=4,
     xlab=expression(paste("Structure diameter [", mu, "m]")),
     ylab="Propability density distribution (PDF)",
     type='l',
     xlim=c(0,170))

points(diameters,percent_Vol,col=rgb(0, 1,0))

legend(100, 0.10, c("calculated from thickness map","ground truth"), col = c(4,3),
       text.col = "black", lty = c(1, 1, 1),# pch = c(NA, NA, 1,1,1),
       merge = TRUE, bg = "white")

#dev.off()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)