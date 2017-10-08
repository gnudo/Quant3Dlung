library(ks)            # for weighted kernel density
library(smoother)     # for smth function (column std)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())  # remove whole workspace

source(paste(dirname(getwd()),'/curvature/R/R_functions.R',sep = ''))

start.time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)

#---------------------------------------------------------------------
# (0) Constants
#---------------------------------------------------------------------
px_size = 1.1E-6             ## Pixel size (10x optics)

nbin = 4096

density_radius  = 0.75       ## Density bandwidth

#---------------------------------------------------------------------
# (1) Load ground truth data
#---------------------------------------------------------------------
curvdata = 'testdata_vertices1_NEW.csv'
curvdata_bin = 'testdata_vertices1_NEW.bin'  ### just use any random name
curvdata_red = 'testdata_vertices1_NEW_red.bin' ### just use any random name
gt_data = read.csv("quant_data.txt", fill = TRUE, sep = "")

R_values = gt_data$R * px_size / 1E-6
n_balls = gt_data$X.px.

Surface = 4 * pi * (R_values )^2
percent_Surface = Surface * n_balls
percent_Surface = percent_Surface/sum(percent_Surface)

#---------------------------------------------------------------------
# (2) Input data
#---------------------------------------------------------------------
print(sprintf("--------- loading dataset ------------"))

X_sample = loadFile(curvdata,curvdata_red,curvdata_bin)
kappa_1 <- -X_sample[,2] * px_size / 1E-6  ## IF THE RAW DATA WAS INVERTED
kappa_2 <- -X_sample[,1] * px_size / 1E-6  ## IF THE RAW DATA WAS INVERTED

areas <- X_sample[,3]
classifier <- X_sample[,4]

rm(X_sample)

#---------------------------------------------------------------------
# (3) Mean curvatures
#---------------------------------------------------------------------
mean_curv_par = ( kappa_1 + kappa_2 )/2

mean_R = (1/mean_curv_par) * px_size / 1E-6

map1 <- kde(mean_R,
            h=density_radius,
            gridsize = nbin,
            xmin=0,
            xmax=100,
            binned=TRUE,
            bgridsize = nbin,
            w=areas/sum(areas)
            )
print('------ finished MeanCurv-KDE // saving results // and restarting Loop ------')
rm(mean_curv_par)
    
#---------------------------------------------------------------------
# (4) Mimic uncertainty due to imprecise curvature calculation
#---------------------------------------------------------------------
x_R = linspace(1,100,nbin)
Surface_R = approx(R_values, percent_Surface, x_R,yleft=0, yright=0)
ys = smth(Surface_R$y,window = 0.06,method = "gaussian") #SMOOTHING

#---------------------------------------------------------------------
# (7) Plots
#---------------------------------------------------------------------
print(sprintf("--------- plotting ------------"))

# output = "~/Desktop/mean_curv.eps"
#cairo_ps(output,height=6, width=6)

par(mfrow=c(1,1))
plot(map1$eval.points,
     map1$estimat,
     col=4,
     xlab=expression(paste("Sphere radius [", mu, "m]")),
     ylab="Propability density distribution (PDF)",
     type='l',
     xlim=c(0,60))
lines(x_R,ys,col=rgb(0, 1,0))

legend(30, 0.12, c("sphere radii from curvature values","ground truth sphere radii (Gauss convolved)"), col = c(4,3,2),
       text.col = "black", lty = c(1, 1, 1),# pch = c(NA, NA, 1,1,1),
       merge = TRUE, bg = "white")

#dev.off()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)