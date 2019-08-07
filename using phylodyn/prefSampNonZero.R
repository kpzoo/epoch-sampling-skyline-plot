# Simulate coalescent process with preferential (proportional) sampling times and varying beta

# Assumptions and modifications
# - removes zeros between epochs by adjusting epoch ends
# - adjustment is to last sample in the previous epoch    
# - adapted Karcher functions to allow time-varying beta
# - uses phylodyn package of Karcher 2016 et al
# - several modified trajectories with _traj2 _traj3 etc

# Clean the workspace and console
closeAllConnections()
rm(list=ls())
cat("\014")  
graphics.off()

# Packages for phylodyn
library("sp")
library("INLA")
library("spam")
library("ape")
library("devtools")
library("phylodyn")

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Function to write simple csv files to correct path
tableWrite <- function(val, name, pathname) {
  # Add path to name
  str0 <- paste(c(pathname, name), collapse = "")
  # Write table
  write.table(val, str0, row.names=FALSE, col.names=FALSE, sep=",")
}

# Modified demographic functions --------------------------------------------------------------------

# Define a logistic trajectory with larger N
N = 10^3
N0 = 0.01*N
logistic_traj2 <- function (t, offset = 0, a = 2) 
{
  t = t + offset
  result = rep(0, length(t))
  result[(t%%12) <= 6] = N0 + N/(1 + exp((3 - (t[(t%%12) <= 6]%%12)) * a))
  result[(t%%12) > 6] = N0 + N/(1 + exp(((t[(t%%12) > 6]%%12) - 12 + 3) * a))
  return(result)
}
# Define a boom-bust with a later changepoint
boombust_traj2 <- function (t, bust = 20, scale = 5000) 
{
  result = rep(0, length(t))
  result[t <= bust] = scale*exp(t[t <= bust] - bust)
  result[t > bust] = scale*exp(bust - t[t > bust])
  return(result)
}
# Define a boom-bust with a later changepoint and an offset
boombust_traj3 <- function (t, bust = 10, scale = 1000, offset = 100) 
{
  result = rep(0, length(t))
  result[t <= bust] = scale*exp(t[t <= bust] - bust) + offset
  result[t > bust] = scale*exp(bust - t[t > bust]) + offset
  return(result)
}

# Define larger exponential
exp_traj2 <- function (t, scale = 10000, rate = 1) 
{
  return(scale*exp(-t*rate))
}
# Define a middling bottleneck
bottleneck_traj2 <- function (t) 
{
  result = rep(0, length(t))
  result[t <= 15] <- 500
  result[t > 15 & t < 40] <- 20
  result[t >= 40] <- 500
  return(result)
}
# Low and high constant populations
unif_traj_low <- function (t, level = 50) 
{
  n = length(t)
  return(rep(level, n))
}
unif_traj_high <- function (t, level = 5000) 
{
  n = length(t)
  return(rep(level, n))
}

# Main code for preferential simulations --------------------------------------------------------------------

# Possible trajectories
trajNames = c('logis', 'exp', 'steep', 'unif_low', 'unif_high', 'boom', 'cyc', 'bottle', 'mesa')
expTrajs = c(2, 6) # exp type trajectories

# Set population true trajectory
type = 1
trajName = trajNames[type]

# Choose trajectory type
trajType = switch(type,
  "1"= logistic_traj2,
  "2"= exp_traj2,
  "3"= steep_cyc_traj,
  "4"= unif_traj_low,
  "5"= unif_traj_high,
  "6"= boombust_traj3,
  "7"= cyclic_traj,
  "8"= bottleneck_traj2,
  "9"= mesa_traj
)
traj = trajType

# Create folder for traj specific results
trajNameSplit = paste(c(trajName, '_test'), collapse = '')
dir.create(file.path(this.dir, trajNameSplit))

# Set sampling interval end and no. samples
if(!is.element(type, expTrajs)){
  # Non-exp trajs so do not decay quickly
  all_samp_end = 48
  nsamps = 2000; ndivs = 100
} else{
  # Exp trajs finish quickly
  if(type == 2){
    all_samp_end = 5
  }
  if(type == 6){
    all_samp_end = 25 
  }
  nsamps = 1000; ndivs = 50
}
#ndivs = ndivs/2

# Uniformly distribute epochs
samps_split = rep(nsamps/ndivs, ndivs)
nsplit = length(samps_split)
tsplit = all_samp_end/nsplit

# Prop constants and sample times for expected sample split
Cprop = rep(0, nsplit)
samp_end = 0; id2 = 0
samp_mx = matrix(list(), nsplit, 1)
len_split = rep(0, nsplit)
samp_times = list()
sampEnd1 = rep(0, nsplit); sampEnd2 = sampEnd1

# Sampling in each epoch, with adjustment of ends to sample times
for(i in 1:nsplit){
  # Split start and end times
  if(i > 1){
    samp_start = sampEnd2[i-1] 
  } else{
    samp_start = 0
  }
  samp_end = samp_start + tsplit
  
  # Different proportionality constant in each split
  Cprop[i] = samps_split[i]/integrate(traj_beta, samp_start, samp_end, traj=traj, beta=1)$value
  
  # Sample times for this split
  sampTemp = pref_sample(traj, c=Cprop[i], lim=c(samp_start, samp_end), beta=1)
  samp_mx[[i]] = sampTemp
  # Combine all sample times
  len_split[i] = length(sampTemp)
  id1 = id2+1; id2 = id1 + len_split[i] - 1
  samp_times[id1:id2] = sampTemp
  # Sample end times (epochs)
  if(i > 1){
    # Start time is last sample in previos epoch
    sampEnd1[i] = sampEnd2[i-1] 
  } else{
    # Or the very first sample time in first case
    sampEnd1[i] = sampTemp[1]
  }
  # End of epoch is last sample within split
  sampEnd2[i] = sampTemp[length(sampTemp)] 
}

# Introduce single samples
samp_times = as.numeric(samp_times)
samps = rep(1, length(samp_times))
            
# Simulate genealogy and get all times
gene = coalsim(samp_times = samp_times, n_sampled = samps, traj = traj, lower_bound = 10, method = "thin")
coal_times = gene$coal_times
coalLin = gene$lineages

# Obtain true trajectory across time
tmax = max(c(samp_end, coal_times))
t = seq(0, tmax, length=40000)
y = traj(t)
# Plot trajectory
#quartz()
#plot(t, y)

# Export data for Matlab
pathf = paste(c(this.dir, '/', trajNameSplit, '/'), collapse = "")
tableWrite(coal_times, 'coaltimes.csv', pathf)
tableWrite(samp_times, 'samptimes.csv', pathf)
tableWrite(coalLin, 'coalLin.csv', pathf)
tableWrite(y, 'trajy.csv', pathf)
tableWrite(t, 'trajt.csv', pathf)
tableWrite(samps, 'sampIntro.csv', pathf)
tableWrite(Cprop, 'beta.csv', pathf)
tableWrite(sampEnd1, 'sampEnd1.csv', pathf)
tableWrite(sampEnd2, 'sampEnd2.csv', pathf)
tableWrite(len_split, 'lensplit.csv', pathf)

# Plot and write tree
tree <-generate_newick(gene)
quartz()
plot(tree$newick, show.tip.label = F)
currDir = this.dir
setwd(file.path(this.dir, trajNameSplit))
write.tree(tree$newick, file="tree.txt")
setwd(currDir)
