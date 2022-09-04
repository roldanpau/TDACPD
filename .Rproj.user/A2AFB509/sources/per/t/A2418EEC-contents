# mixture_of_normals.R
#
# Generate random samples from a Mixture of Normals distribution.
# By definition, a mixture of normals has density function 
#    f(x) = p*g1(x) + (1-p)*g2(x)
# where g1 and g2 are normal and p is the mixture parameter.
#
# For simplicity we assume that g1 and g2 are symmetrically located 
# with respect to the origin, and that they have the same variance.
#
# In the code below, I have chosen the means mu1, mu2 
# and the standard deviations sd1, sd2 in such a way that the variance 
# of the mixture distribution is precisely equal to 1.5, 
# to match previous experiments by Elisheva.

# Working directory is where this script is at.
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Load required packages
library(R.utils)  # for printf
library(TDA)      # for ripsDiag

# Parameters of probability distributions.
# g1 follows a normal with mean mu1 and SD sd1, while
# g2 follows a normal with mean mu2 and SD sd2.
mu1=sqrt(3 - 1/8)
mu2=-sqrt(3 - 1/8)
sd1=sqrt(1/8)
sd2=sqrt(1/8)

# Parameter of mixture of normals
# f(x) = p*g1(x) + (1-p)*g2(x)
p=0.5

N=10000 # Number of samples from the mixture distribution (time series length)

#Sample N random uniforms U
U =runif(N)

#Variable to store the samples from the mixture distribution                                             
x = rep(NA,N)

#Sampling from the mixture
for(i in 1:N){
  if(U[i]<p){
    x[i] = rnorm(1,mu1,sd1)
  }else{
    x[i] = rnorm(1,mu2,sd2)
  }
}

#Density plot of the random samples
plot(density(x),main="Density Estimate of the Mixture Model")

# Check empirically that sample variance of mixture is 1.5.
print(var(x))

##### Find persistence diagrams #####

# max dimension of the homological features to be computed. 
# (e.g. 0 for connected components, 1 for connected components and loops, ...)
maxdimension = 0  

max_scale = 10  # maximum value of the rips filtration.

win_size = 100            # size of windows
n_windows = N-win_size+1  # number of windows

Diags.rips.ls = list()    # list of persistence diagrams (one for each window)

step <- 10
spot <- seq(from=1, to=n_windows, by=step)

print("Computing persistence diagrams...")
pb <- txtProgressBar(min = 0, max = length(spot), style = 3)
for(i in 1:length(spot)){
  Diags.rips.ls[[i]] = 
    ripsDiag(X = x[spot[i]:(spot[i]+win_size-1)], maxdimension, max_scale, library = "Dionysus")
  
  # update progress bar
  setTxtProgressBar(pb,i)
}
#plot(Diags.rips.ls[[1]]$diagram)
print(Diags.rips.ls[[1]]$diagram)


##### Calculate the max norm of persistence diagrams #####

norms = vector(length=length(Diags.rips.ls)) # collects values of norm at a particular window

print("Computing the max norm of persistence landscapes...")

for (i in 1:length(Diags.rips.ls)){
  diagram <- Diags.rips.ls[[i]]$diagram
  
  maxim=0
  # Iterate over all diagram features, excluding the first one (which corresponds to max_scale)
  for(j in 1:nrow(diagram)) {   
    if(diagram[j,1]==0 && diagram[j,3]!=max_scale) { # 0-dimensional feature
      maxim = max(maxim,diagram[j,3])
    }
  }
  #print(minim)
  #print(maxim)
  
  if(maxim==0)  # There is no 0-dimensional homology feature in this window
  {
    norms[i]=0
  }
  else 
  {
    norms[i]=maxim
  }
}
plot(norms)
mean(norms)
var(norms)
