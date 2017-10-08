#---------------------------------------------------------------------
# Load Datasets from a text file
#---------------------------------------------------------------------
# The file must be an ASCII-txt file with with each dataset heading
# preceeded by "##". At the end an empty newline must be added.
# Sample file:
#    ## Datasets 1
#    Dataset1/Path
#    Dataset2/Path
#
#    ## Datasets 2
#    ...
#
# TODO: use a more generic fileformat (XML or similar)
#---------------------------------------------------------------------
loadDatasetsFromTextfile <- function(filepath) {
  datasets = list()  # from other file
  datasets_numbers = list()  # from other file
  
  conn <- file(filepath,open="r")
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
  } ## LOOP through each Dataset
  close(conn)
  return( list(datasets,datasets_numbers) )
}


#---------------------------------------------------------------------
# Compose datasets
#---------------------------------------------------------------------
# Either from command line arguments or hardcoded
#---------------------------------------------------------------------
composeDatasets <- function(args,datasrc,datasets) {
  if ( length(args) > 0 ) {
    csv_file = args[1]
    csv_file_red_tmp = strtrim(csv_file, nchar(csv_file)-4)
    csv_file_red = paste(csv_file_red_tmp,'_reduced.bin', sep = '', collapse = NULL)
    csv_file_bin = paste(csv_file_red_tmp,'.bin', sep = '', collapse = NULL)
    if ( length(args) > 1 ) {
      for ( ii in 2:length(args) ) {
        csv_file = c(csv_file,args[ii])
        csv_file_red_tmp = strtrim(args[ii], nchar(args[ii])-4)
        csv_file_red = c(csv_file_red,paste(csv_file_red_tmp,'_reduced.bin', sep = '', collapse = NULL))
        csv_file_bin = c(csv_file_bin,paste(csv_file_red_tmp,'.bin', sep = '', collapse = NULL))
      } 
    }
  } else {
    csv_file = paste(datasrc,datasets,'.csv', sep = '', collapse = NULL)
    csv_file_red = paste(datasrc,datasets,'_reduced.bin', sep = '', collapse = NULL)
    csv_file_bin = paste(datasrc,datasets,'.bin', sep = '', collapse = NULL)
  }
  list("csv_file" = csv_file, "csv_file_red" = csv_file_red, "csv_file_bin" = csv_file_bin)
}


#---------------------------------------------------------------------
# Reduce data if necessary
#---------------------------------------------------------------------
# We treat 3 cases and 1 condition:
# (i)   we have a reduced binary --> we load this
# (ii)  we have a binary-CSV --> then we load that one
# (iii) we have no binary files --> then we load original CSV 
# (1)   we want to save a reduced file
#---------------------------------------------------------------------
loadFile <- function(csv_file,csv_file_red,csv_file_bin) {
  if ( file.exists(csv_file_red) ) {
    # Case 1
    load(csv_file_red)
    print("loaded reduced CSV file")
  } else {
    if ( file.exists(csv_file_bin) ) {
      # Case 2
      load(csv_file_bin)
      print("loaded binary CSV file")
    } else {
      # Case 3
      dat = read.csv(csv_file, header = FALSE)
      X_sample <- cbind( dat[,1], dat[,2], dat[,3], dat[,4] )
      print("loaded original CSV file")
    
      # Always save binary CSV
      save(X_sample, file=csv_file_bin, compress=FALSE)
      print("saved binary CSV file")
    }
  
    if ( exists("n_reduced") ) {
      # Condition 1
      s <- sample( nrow(X_sample), size=n_reduced )  # take "size" rows out of nrow rows from matrix
      X_sample <- X_sample[s,]  # create a new matrix
    
      # Write CSV in R
      save(X_sample, file=csv_file_red, compress=FALSE)
      print("saved reduced CSV file")
    }
  }

  if ( exists("n_plotting") ) {
    X_sample <- X_sample[1:n_plotting,]
  }
  X_sample
}


#---------------------------------------------------------------------
# 2D trapezoidal integration
#---------------------------------------------------------------------
trapz2D <- function(f_num, x, y, a_x, b_x, a_y, b_y)  {
  x_reg <- x >= a_x & x <= b_x
  y_reg <- y >= a_y & y <= b_y
  
  fun <- function(x,x_tmp) trapz(x_tmp,x)
  
  int_x <- apply(f_num[y_reg,x_reg], 1, fun,x_tmp=x[x_reg] )
  
  #int_x <- apply(f_num[y_reg,x_reg], 1, function(var) trapz(x[x_reg],var) )
  int2d <- trapz(y[y_reg],rev(int_x))
  int2d
}


#---------------------------------------------------------------------
# Plot kernel density
#---------------------------------------------------------------------
plotDensityMat <- function(x,y,fhat,wished_range) {
  
  if(missing(wished_range)) {
    wished_range = max(fhat)
  }
  
  XYmesh <- meshgrid( x, y )
  kappa_1 = c(XYmesh$Y)
  kappa_2 = c(XYmesh$X)
  
  df <- data.frame(kappa_1,kappa_2)
  max_counts = max(fhat)
  colramp=colorRampPalette(c("black", "white"))
  xy1 <- xy.coords(kappa_1, kappa_2)
  select <- is.finite(xy1$x) & is.finite(xy1$y)
  x1_tmp <- cbind(xy1$x, xy1$y)[select, ]
  xbin <- cut(x1_tmp[,1], mkBreaks(x), labels = FALSE)
  ybin <- cut(x1_tmp[,2], mkBreaks(y), labels = FALSE)
  dens <- fhat[cbind(xbin, ybin)] #* X_sample[,3]
  dens[is.na(dens)] <- 0
  len1 = length(dens)
  
  ## transform densities to colors
  print("transforming densities to colors")
  
  ## SCALE COLORBAR
  if ( exists("wished_range") & (max_counts > wished_range) ) {
    dens[dens>wished_range] = 0
  }
  
  colpal <- cut(dens, length(dens), labels = FALSE)
  
  cols   <- rep(NA_character_, length(select))
  cols[select] <- colramp(len1)[colpal]
  
  if (length(nbin) == 1)
    nbin <- c(nbin, nbin)
  if (!is.numeric(nbin) || length(nbin) != 2)
    stop("'nbin' must be numeric of length 1 or 2")
  
  
  df$dens <- round((c(fhat)/wished_range)*255)+1
  
  if ( exists("wished_range") & (max_counts < wished_range) ) {
    fac_scal = wished_range/max_counts
  } else {
    fac_scal = 1
  }
  
  print(sprintf("Max counts: %4.1f", max_counts ))  ## is not really color-coded because we do mbreaks before
  
  cols <-  colorRampPalette(c("#FFFFFF","#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))( as.integer(256)*fac_scal) #)  ## must be 256 because of col2rgb
  
  df$col <- cols[df$dens]
  
  par(mfrow=c(1,1))
  plot(df$kappa_2[order(df$dens)]~df$kappa_1[order(df$dens)],
       data=df[order(df$dens),],
       pch=20, col=col,
       cex=0.5,  # size of markers
       # cex=5,  # >>> to visualize very small bright peaks
       #xlim=c(-0.4,0.1), ylim=c(-0.2,0.3),
       #xlim=c(-0.2,0.1), ylim=c(-0.1,0.2),
       # xlim=c(-0.25,0.07), ylim=c(-0.05,0.27), ## probably best one
       xlim=c(-0.20,0.07), ylim=c(-0.05,0.22), ## probably best one
       # xlim=c(-0.3,0.3), ylim=c(-0.3,0.3),  ## just test
       asp=1,    # aspect ratio
       main="Interface Shape Distributions",
       # xlab=expression( kappa[1]~ textstyle(over(1,paste("[", mu, "m]"))) ),
       xlab=expression(paste(kappa[1]," [", mu, "m"^-1,"]")),
       # ylab=expression( kappa[2]~ textstyle(over(1,paste("[", mu, "m]"))) ),
       ylab=expression(paste(kappa[2]," [", mu, "m"^-1,"]")),
       cex.lab=1.5,  # axis labels >> bigger
       xaxs = "i", yaxs = "i"  # remove margins
  )
  
  par(pty='s')
  par(mar = c(5,4,4,5) + .1)  # this is the margin where to place the bar
  
  # Add lines for the mean of X and Y
  abline(h=0, v=0, col="black", lwd=1.5)
  abline(0, 1, col="black", lwd=1.5, lty=2)  # lty = dashed
  abline(0, -1, col="black", lwd=1.5, lty=2)  # lty = dashed
  
  # Add colorbar to the image
  image.plot(fhat,
             legend.only=T,
             col = cols,
             nlevel = 256, ## if changing here doesn't effect scale bar
             #legend.args=list(text='Counts', side=4, font=2, line=2.5, cex=1),
             add = FALSE,  ## important! otherwise colorbar is in the plot
             #legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8)
             zlim=c(0,wished_range)  ## gives an explicit limit for the colorbar
  )
}


#---------------------------------------------------------------------
# Small helper function for calculating middle points after CUT
#---------------------------------------------------------------------
mkBreaks <- function(u) {
  u - diff(range(u))/(length(u)-1)/2
}