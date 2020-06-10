options(stringsAsFactors = FALSE)

# Range of values that will be simulated
npts <- 25
slist <- seq(0.45, 0.7, length.out = 25)
rlist <- seq(0, 3, length.out = npts + 1)
rlist <- rlist[-1]
miglist <- c(0.1, 1)

# Note: The output of these lines is pasted into the script that runs the stochastic model
slist
rlist
miglist

# Need to save the vectors in the same format as they are printed
# (trailing zeros), to be able to write file names correctly
slist.text <- sprintf("%0.7f", slist)
rlist.text <- sprintf("%0.2f", rlist)
miglist.text <- c("0.1", "1")

# Write all compbinations of the parameters
prms <- expand.grid(slist, rlist, miglist)
prms.text <- expand.grid(stext = slist.text, rtext = rlist.text, migtext = miglist.text, stringsAsFactors = FALSE)

# Function to code the output of a simulation
#  0 coexistence at the end of the simulation
# +1 WT won (no Drive in the end)
# -1 Drive won (no WT in the end)
# The function also returns the time at which the simulation stopped, 
# i.e. extinction time, if no coexistence
characterizeSim <- function(stext, rtext, migtext){
  a <- as.matrix(read.csv(paste0("../data/stoch_s-", stext, "_r-", rtext, "_mig-", migtext, ".csv"), skip = 1, header = FALSE))
  nsites <- (ncol(a) - 1)/2
  iO <- 2*(1:nsites)
  iD <- iO + 1
  
  # Simulation stops when one type disappears, or when maxtime is reached
  # Compute abundances of the two types
  totO <- sum(a[iO])
  totD <- sum(a[iD])
  
  if(totO > 0 && totD > 0){
    # Coexistence
    typeOutcome <- 0
  }else{
    if(totO > 0){
      # WT won
      typeOutcome <- +1
    }else{
      # Drive won
      typeOutcome <- -1
    }
  }
  # Return type of outcome and final time
  out <- c(typeOutcome, c(a[1]))
  return(out)
}

# Characterize all the simulations
# Initiate output vector
v <- matrix(0, ncol = 2, nrow = nrow(prms))
# Loop on all simulations
for(i in seq_len(nrow(prms))){
  cat(i, "")
  pprm <- prms.text[i,]
  
  v[i,] <- do.call(characterizeSim, as.list(pprm[1:3]))
}

# Add the result to the table of parameters
vp <- cbind(prms, v, 1/v[,2])
vp.df <- as.data.frame(vp)
names(vp.df) <- c("s", "r", "mig", "type", "time", "speed")

# PLOTTING

# Function to make color transparent
# Modified from 
# Source: https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent <- function(col, alpha=0.5){
  # We accept alpha greater than 1, and turn it into 1
  if(alpha > 1){
    alph <- 1
  }else{
    alph <- alpha
  }  

  # Turn alpha into value in (0, 255)
  alph <- floor(255*alph)  
  newColor <- col2rgb(col = col, alpha = FALSE)
  
  # Apply to our colors
  newColor <- rgb(red = newColor[1], green = newColor[2], blue = newColor[3], 
                  alpha = alph, maxColorValue = 255)
  return(newColor)
}

# Colors of the winning type
colWTwins <- "#F0B400" # Blue
colDwins <- "#000082" # Yellow

themig <- 0.1
subvp <- vp.df[vp.df$mig == themig,]

# Defining max for colors using quantiles (or not, 1 = 100%)
thrspeed <- quantile(abs(subvp$type*subvp$speed), probs = c(1))
thrspeed

# Compute vector of colors
# We color-code the speeds using color intensities
# and direction (left or right) using color value
cols <- rep(0, nrow(subvp))
cols[subvp$type == 0] <- "gray"
cols[subvp$type > 0] <- vapply(subvp[subvp$type > 0, "speed"]/thrspeed, FUN.VALUE = "a", FUN = makeTransparent, col = colWTwins)
cols[subvp$type < 0] <- vapply(subvp[subvp$type < 0, "speed"]/thrspeed, FUN.VALUE = "a", FUN = makeTransparent, col = colDwins)

# Plot!
par(las = 1)
plot(subvp$s, subvp$r, col = cols, pch = 16, cex = 1.2, axes = FALSE, 
     xlab = "s", ylab = "r")
axis(1)
axis(2, at = seq(0, max(rlist), by = 0.5))

# Add curve for eradication drive
ss <- seq(min(slist), max(slist), by = 0.01)
lines(ss, ss/(1-ss), lwd = 2)

# Add limits for threshold drives
abline(v = 0.5)
abline(v = 2/3)


#########################################################
# TSN

# Change parameters
slist <- seq(0.45, 0.9, length.out = 45)
slist.text <- sprintf("%0.7f", slist)
slist
miglist <- c(0.1, 0.5, 1)
miglist.text <- c("0.1", "0.5", "1")

prms <- expand.grid(slist, miglist)

prms.text <- expand.grid(stext = slist.text, migtext = miglist.text, stringsAsFactors = FALSE)

# Function to code the output of a simulation
# (TSN)
characterizeSimTSN <- function(stext, migtext, K = 1000){
  a <- as.matrix(read.csv(paste0("../data/TSN_s-", stext, "_mig-", migtext, ".csv"), skip = 1, header = FALSE))
  nsites <- (ncol(a) - 1)
  iD <- 1 + (1:nsites)
  
  # Simulation stops when one type disappears, or when maxtime is reached
  # Compute abundances of the two types
  totD <- sum(a[iD])
  totO <- K * nsites - totD
  
  if(totO > 0 && totD > 0){
    # Coexistence
    typeOutcome <- 0
  }else{
    if(totO > 0){
      # WT won
      typeOutcome <- +1
    }else{
      # Drive won
      typeOutcome <- -1
    }
  }
  # Return type of outcome and final time
  out <- c(typeOutcome, c(a[1]))
  return(out)
}

# Characterize all the simulations
# Initiate output vector
vTSN <- matrix(0, ncol = 2, nrow = nrow(prms.text))
# Loop on all simulations
for(i in seq_len(nrow(prms.text))){
  cat(i, "")
  pprm <- prms.text[i,]
  
  vTSN[i,] <- do.call(characterizeSimTSN, as.list(pprm[1:2]))
}

# Add the result to the table of parameters
vpTSN <- cbind(prms, vTSN, 1/vTSN[,2])
vpTSN.df <- as.data.frame(vpTSN)
names(vpTSN.df) <- c("s", "mig", "type", "time", "speed")

## PLOTTING
themig <- 0.5
subvp <- vpTSN.df[vpTSN.df$mig == themig,]

# Defining max for colors using quantiles (or not, 1 = 100%)
thrspeed <- quantile(abs(subvp$type*subvp$speed), probs = c(1))
thrspeed

# Compute vector of colors
# We color-code the speeds using color intensities
# and direction (left or right) using color value
cols <- rep(0, nrow(subvp))
cols[subvp$type == 0] <- "gray"
cols[subvp$type > 0] <- vapply(subvp[subvp$type > 0, "speed"]/thrspeed, FUN.VALUE = "a", FUN = makeTransparent, col = colWTwins)
cols[subvp$type < 0] <- vapply(subvp[subvp$type < 0, "speed"]/thrspeed, FUN.VALUE = "a", FUN = makeTransparent, col = colDwins)

# Plot!
par(las = 1)
plot(subvp$s, 0.2+0*subvp$s, col = cols, pch = 16, cex = 1.2, axes = FALSE, 
     xlab = "s", ylab = "", ylim = c(0, 0.3))
axis(1, pos = 0)

# Add limits for threshold drives
abline(v = 0.5)
abline(v = 0.7)

