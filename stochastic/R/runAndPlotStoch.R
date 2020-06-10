setwd("/Users/flo/Documents/Work/Projects/1_EnCours/2019_LeoDriveInSpace/scripts/stochastic/R")

# Compile the C code
system("gcc ../stochwave.c -lm -o stochwave")

runSim <- function(K, s, r, mig, tmax, nreps = 1){
  # Prepare command to run the simulation
  cmd <- paste0("./stochwave ", K, " ", s, " ", r, " ", mig, " ", tmax, " ", nreps, " > ../data/sim_K-", K, "_s-", s, "_r-", r, "_mig-", mig, ".csv" )
  # Run the simulation
  system(cmd)
  
  # Load the data 
  data <- as.matrix(read.csv(paste0("../data/sim_K-", K, "_s-", s, "_r-", r, "_mig-", mig, ".csv"), skip = 1, header = FALSE))
  
  return(data)
}



colWT <- "black"
colD <- "red"

plotEndSim <- function(data){
  # If running this while the simulation is running, there will be NAs at the end
  # -> remove them
  if(any(is.na(data))) data <- data[-nrow(data),]
  
  # Extract number of sites
  nsites <- (ncol(data) - 1)/2
  # Indices of the two types (O: WT, D: drive)
  iO <- 2*(1:nsites)
  iD <- iO + 1
  
  # Max pop size ever, to scale the plots
  ymax <- max(data[,-1])
  
  # Plot the data
  itime <- nrow(data)
  plot(1:nsites, data[itime, iO], ylim = c(0, ymax), type = "l", col = colWT, xlab = "Space", ylab = "Pop size")
  lines(1:nsites, data[itime, iD], col = colD)
}


for(i in dev.list()) dev.off()
mydata <- runTSN(1000, 0.7, 0.5, 10000)
plot(mydata[nrow(mydata), -1]/1000, ylim = c(0, 1))


plotMovie <- function(data, fname){
  # If running this while the simulation is running, there will be NAs at the end
  # -> remove them
  if(any(is.na(data))) data <- data[-nrow(data),]
  
  # Extract number of sites
  nsites <- (ncol(data) - 1)/2
  # Indices of the two types (O: WT, D: drive)
  iO <- 2*(1:nsites)
  iD <- iO + 1
  
  # Max pop size ever, to scale the plots
  ymax <- max(data[,-1])
  
  # Function to plot the data
  plotData <- function(i){
    par(las = 1)
    xpar <- 3.25
    par(mar = c(xpar, xpar, xpar, 1))
    par(mgp = c(2, 0.5, 0))
    lwdOD <- 4
    plot(1:nsites, data[i, iO], ylim = c(0, ymax), type = "l", 
         col = colWT, lwd = lwdOD,
         xlab = "Space", ylab = "Pop size")
    lines(1:nsites, data[i, iD], col = "red", lwd = lwdOD)
    title(main = paste0("t = ", data[i, 1]))
  }
  
  wmax <- nchar(as.character(nrow(data)))
  for(itime in 1:nrow(data)){
    png(paste0("Pics/", fname, "_", formatC(itime, width = wmax, format = "d", flag = "0")
, ".png"), width = 600, height = 450, pointsize = 16)
    plotData(itime)
    dev.off()
  }
}

data <- runSim(1000, 0.6, 1.0, 0.1, 10000)

png("Pics/plot.png", width = 300, height = 225)
par(las = 1)
xpar <- 3.25
par(mar = c(xpar, xpar, xpar, 1))
par(mgp = c(2, 0.5, 0))
plotEndSim(data)
title("x")
dev.off()


plotMovie(data, "test")



runAndPlotMovie <- function(K, s, r, mig, tmax, nreps = 1){
  cat('Running the simulation...\n')
  # Run the simulation
  data <- runSim(K, s, r, mig, tmax, nreps)
  
  # Output name for the figures
  fname <- paste0("simplot_K-", K, "_s-", s, "_r-", r, "_mig-", mig)

  # Remove previous files to avoid problems when creating movie
  system(paste0("rm Pics/", fname, "*.png"))
  
  cat('Plotting the frames...\n')
  # Plot the frames
  plotMovie(data, fname)
  
  cat('Converting into movie...\n')
  # Convert the frames into a movie
  cmd <- paste0("convert -delay 2 -loop 0 Pics/", fname, "_*.png Movies/anim_", fname, ".mpeg")
  system(cmd)
  cat('done!\n')
  
  # Open the movie
  system(paste0("open Movies/anim_", fname, ".mpeg") )
}

runAndPlotMovie(1000, 0.6, 1.08, 0.1, 10000, 1)
runAndPlotMovie(1000, 0.6583, 2.04, 0.1, 10000, 1)

runAndPlotMovie(1000, 0.522, 1.5, 0.1, 10000, 1)

runAndPlotMovie(1000, 0.522, 0.24, 0.1, 10000, 1)

runAndPlotMovie(1000, 0.522, 0.84, 10000, 1)



runAndPlotMovie(100, 0.6, 1.08, 0.1, 10000, 1)
runAndPlotMovie(100, 0.6583, 2.04, 0.1, 10000, 1)
runAndPlotMovie(100, 0.522, 1.5, 0.1, 10000, 1)
runAndPlotMovie(100, 0.522, 0.24, 0.1, 10000, 1)
runAndPlotMovie(100, 0.522, 0.84, 0.1, 10000, 1)

runAndPlotMovie(10, 0.6, 1.08, 0.1, 10000, 1)
runAndPlotMovie(10, 0.6583, 2.04, 0.1, 10000, 1)
runAndPlotMovie(10, 0.522, 1.5, 0.1, 10000, 1)
runAndPlotMovie(10, 0.522, 0.24, 0.1, 10000, 1)
runAndPlotMovie(10, 0.522, 0.84, 0.1, 10000, 1)

# Plot snapshots
plotData(1)
plotData(floor(nrow(data)/3))
plotData(floor(2*nrow(data)/3))
plotData(nrow(data))

# Plot all tMoviesimes for movie
for(i in 1:nrow(data)){
  plotData(i)
  cat(i )
}


#########################################################################
# TSN
system("gcc ../stochTSN.c -lm -o stochTSN")

runTSN <- function(K, s, mig, tmax, nreps = 1){
  # Prepare command to run the simulation
  cmd <- paste0("./stochTSN ", K, " ", s, " ", mig, " ", tmax, " ", nreps, " > ../data/TSN_K-", K, "_s-", s, "_mig-", mig, ".csv" )
  # Run the simulation
  system(cmd)
  
  # Load the data 
  data <- as.matrix(read.csv(paste0("../data/TSN_K-", K, "_s-", s, "_mig-", mig, ".csv"), skip = 1, header = FALSE))
  
  return(data)
}

#--------------------------------------------

plotMovieTSN <- function(data, fname, K = 1000){
  # If running this while the simulation is running, there will be NAs at the end
  # -> remove them
  if(any(is.na(data))) data <- data[-nrow(data),]
  
  # Extract number of sites
  nsites <- ncol(data) - 1
  # Indices of the drive indiv
  iD <- 1 + (1:nsites)

  # Max pop size ever, to scale the plots
  #ymax <- max(data[,-1])
  ymax <- K
  
  # Function to plot the data
  plotData <- function(i){
    par(las = 1)
    xpar <- 3.25
    par(mar = c(xpar, xpar, xpar, 1))
    par(mgp = c(2, 0.5, 0))
    lwdOD <- 4
    plot(1:nsites, K - data[i, iD], ylim = c(0, ymax), type = "l", 
         col = colWT, lwd = lwdOD,
         xlab = "Space", ylab = "Pop size")
    lines(1:nsites, data[i, iD], col = "red", lwd = lwdOD)
    title(main = paste0("t = ", data[i, 1]))
  }
  
  wmax <- nchar(as.character(nrow(data)))
  for(itime in 1:nrow(data)){
    png(paste0("Pics/", fname, "_", formatC(itime, width = wmax, format = "d", flag = "0")
               , ".png"), width = 600, height = 450, pointsize = 16)
    plotData(itime)
    dev.off()
  }
}

#-----------------------------------------------

runAndPlotMovieTSN <- function(K, s, mig, tmax, nreps = 1){
  # Run the simulation
  data <- runTSN(K, s, mig, tmax, nreps)
  
  # Output name for the figures
  fname <- paste0("TSNplot_K-", K, "_s-", s, "_mig-", mig)
  
  # Remove previous files to avoid problems when creating movie
  system(paste0("rm Pics/", fname, "*.png"))
  
  # Plot the frames
  plotMovieTSN(data, fname, K)
  
  # Convert the frames into a movie
  cmd <- paste0("convert -delay 2 -loop 0 Pics/", fname, "_*.png Movies/anim_", fname, ".mpeg")
  system(cmd)
  
  # Open the movie
  system(paste0("open Movies/anim_", fname, ".mpeg") )
}

#------------------------------------------------------
runAndPlotMovieTSN(1000, 0.522, 0.5, 10000, 1)
runAndPlotMovieTSN(1000, 0.65, 0.5, 10000, 1)
runAndPlotMovieTSN(1000, 0.8, 0.5, 10000, 1)

