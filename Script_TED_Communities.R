#Install flowCore via bioconductor 
##install.packages("flowCore", repos="http://bioconductor.org/biocLite.R")
library(flowCore)

#Packages for plots 
library(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggcyto")
library(ggcyto)

#Packages for TED
library(geozoo)
library(flexmix)
library(geometry)


#Install flowCore via bioconductor - cf other script for library
#Import fcs file

##find the work directory 
getwd()
##path to the file
setwd("C:/Users/Elif Ka/Documents/DataFlowCytometer/LabExp4th")
##print all the files in this directory
files <- dir()
dir()

#summarizing the data
totes <- list()
files <- files[1]
frames <- list() #contains the data, frames[[i]] contains all the measuemts of one species
for (i in c(1:72)){
  frame <- lapply(files, read.FCS,emptyValue = F, dataset = i, alter.names = T) #make a frame with several FCS
  frames[[i]] <- as(frame, "flowSet") #loading all the datas
  totes[[i]] <-frame
}
raw.data <- as(unlist(totes), "flowSet")
cols <- colnames(raw.data)

#get rid of negative values in sample (=measure failures)
mat<-matrix(rep(c(0.5,Inf),11),ncol = 11, dimnames = list(c("min", "max"), cols))
noneg <- rectangleGate(filterId="positive", .gate=mat)
noneg.fs <- Subset(raw.data, noneg)

#log transform the data set
log.trans <- logTransform()
myTrans <- transformList(cols, log.trans) #log transform all parameters
log.fs <- transform(noneg.fs, myTrans)

#Plots to more precisely define the edges - cf other script for library
#autoplot (log.fs[10:18], "RED.B.HLin", "YEL.B.HLin") #for species sample 1 to 9 
#autoplot (log.fs[10:18], "SSC.HLin", "SSC.ALin") #for singlets
#autoplot (log.fs[10:18], "FSC.HLin", "SSC.HLin") #for debris

#get rid of debris 
##Here: -0.5, -0.5, 1.5, 1.5, 3, -0.5, -0.5, 3
#autoplot (log.fs[1:18], "FSC.HLin", "SSC.HLin") #for debris
edges.deb <- matrix(c(2,2,3,3,5,1,1,5), ncol = 2, nrow = 4)
colnames(edges.deb) <- c("FSC.HLin", "SSC.HLin")
debris <- polygonGate(filterId="Debris", .gate=edges.deb)
##show without debris 
without.debris.fs <- Subset(log.fs[1:18], !debris)
#autoplot (without.debris.fs[1:18], "FSC.HLin", "SSC.HLin")

#only take singlets
#autoplot (log.fs[1:18], "SSC.HLin", "SSC.ALin") #for singlets
edges.sing <- matrix(c(1, 2, 5, 4, 2, 1, 4, 5), ncol = 2, nrow = 4)
colnames(edges.sing) <- c("SSC.HLin", "SSC.ALin")
singlets <- polygonGate(filterId="Singlets", .gate=edges.sing)

# Sp filtering - Subsetting the different sp of algae - /w kmeans function
rawsp.filter <- !debris & singlets 
raw.sp <- Subset (log.fs, rawsp.filter)
#autoplot (raw.sp[10:18], "SSC.HLin", "SSC.ALin")

#purify more the species
#autoplot(raw.sp[10:18], "RED.B.HLin", "YEL.B.HLin") 
#autoplot(raw.sp[28:36], "RED.B.HLin", "YEL.B.HLin")
#autoplot(raw.sp[46:54], "RED.B.HLin", "YEL.B.HLin")
#autoplot(raw.sp[64:72], "RED.B.HLin", "YEL.B.HLin")
edges.species <- matrix(c(1.5, 1.5, 4, 4, 5, 0, 0, 5), ncol = 2, nrow = 4)
colnames(edges.species) <- c("RED.B.HLin", "YEL.B.HLin")
species <- polygonGate(filterId="species", .gate=edges.species)

## apply
sp.filter <- !debris & singlets & species
sp <- Subset (log.fs, sp.filter)
#autoplot(sp[10:18], "RED.B.HLin", "YEL.B.HLin") 
#autoplot(sp[28:36], "RED.B.HLin", "YEL.B.HLin")
#autoplot(sp[46:54], "RED.B.HLin", "YEL.B.HLin")
#autoplot(sp[64:72], "RED.B.HLin", "YEL.B.HLin")
# I take off the param 9 and 11 (-> 9 and 10)
sp <- sp [, -9]
sp <- sp [, -10]
print(sp)

#autoplot(raw.sp[10:18], "RED.B.HLin", "YEL.B.HLin")

species.filter<- kmeansFilter("RED.B.HLin"=c("Species1","Species2"),
                              filterId="Species")

#############################################
#FOR C0, COMBI sample1                      #
#############################################
spkF <- split(sp[[10]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)


#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample2                      #
#############################################
spkF <- split(sp[[11]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample3                      #
#############################################
spkF <- split(sp[[12]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample4                      #
#############################################
spkF <- split(sp[[13]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample5                      #
#############################################
spkF <- split(sp[[14]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample6                      #
#############################################
spkF <- split(sp[[15]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample7                      #
#############################################
spkF <- split(sp[[16]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample8                      #
#############################################
spkF <- split(sp[[17]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample                      #
#############################################
spkF <- split(sp[[18]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample9                      #
#############################################
spkF <- split(sp[[28]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample9                      #
#############################################
spkF <- split(sp[[29]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C0, COMBI sample9                      #
#############################################
spkF <- split(sp[[30]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[31]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[32]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[33]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[34]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[35]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[36]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[46]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[47]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[48]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[49]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[50]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[51]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[52]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[53]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[54]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[64]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[65]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[66]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[67]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[68]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[69]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[70]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)

#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[71]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)


#############################################
#FOR C1, COMBI sample9                      #
#############################################
spkF <- split(sp[[72]], species.filter)

Species1 <- exprs(spkF$Species1)
Species2 <- exprs(spkF$Species2)
Species <- rbind2(Species1, Species2)

#####################
##    TED index    ## sample 1
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

#############################
# List of possible REFERENCES
#############################
# Open biggest sample
# Change filename to name of file with maximum points
traits.max <- Species
# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
nrow(traits.max)
ncol(traits.max)

max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traits.max)