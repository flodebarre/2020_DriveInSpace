# Load the data
data <- as.matrix(read.csv("../simtext.csv", skip = 1, header = FALSE))

plotPic <- function(filename){
  # Load data
  data <- as.matrix(read.csv(filename, skip = 1, header = FALSE))
  # If running this while the simulation is running, there will be NAs at the end
  # -> remove them
  if(any(is.na(data))) data <- data[-nrow(data),]
  
  # Extract number of sites, and indices of the two types
  nsites <- (ncol(data) - 1)/2
  iO <- 2*(1:nsites)
  iD <- iO + 1
  
  # Max pop size ever, to scale the plots
  ymax <- max(data[,-1])
  
  # Function to plot the data
  plotData <- function(i){
    par(las = 1)
    plot(1:nsites, data[i, iO], ylim = c(0, ymax), type = "l", col = "black", xlab = "Space", ylab = "Pop size")
    lines(1:nsites, data[i, iD], col = "red")
    lines(1:nsites, data[i, iD] + data[i, iO], col = "grey", lty = 3)
    title(main = paste0("t = ", data[i, 1]))
  }
  
  # Plot snapshots
  layout(matrix(1:4, ncol = 1, byrow = TRUE))
  plotData(1)
  plotData(floor(nrow(data)/3))
  plotData(floor(2*nrow(data)/3))
  plotData(nrow(data))
}

# Parameter set 2) 
# fw = 0.7, omegaH = 0.1, r = 5 (no extinction)
plotPic("../simCI_2.csv")

# 1) 
# fw = 0.7, omegaH = 0.1, r = 0.5 (no extinction)
plotPic("../simCI_1.csv")

# 3)
# fw = 0.6, omegaH = 0.1, r = 0.5 (extinction)
plotPic("../simCI_3.csv")

# Plot all times for movie
for(i in 1:nrow(data)){
  plotData(i)
  cat(i )
}

