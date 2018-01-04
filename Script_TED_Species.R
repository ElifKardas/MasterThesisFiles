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
files <- files[4]
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
#autoplot (without.debris.fs[1:18], "RED.B.HLin", "YEL.B.HLin")

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



species.filter<- kmeansFilter("RED.B.HLin"=c("Species1","Species2"),
                              filterId="Species")

#--------------------------------------FONCTION----------------------------------------------

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



stockAS1<-array(dim=c(12,2,2))
stockAS2<-array(dim=c(12,2,2))
stockAS3<-array(dim=c(12,2,2))

#--------------------------------------LOOP1sp1----------------------------------------------

samples<-c(10:12, 28:30, 46:48, 64:66)

j=1

for (s in samples){
  
  spkF <- split(sp[[s]], species.filter)
  Species1 <- exprs(spkF$Species1)
  
  traits.max <- Species1
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(paste("nrow=",nrow(traits.max)))
  
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
  stockAS1[j, 1, 1]<-s
  stockAS1[j, 2, 1]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS1))) #changed stockC in stockAS1#elif
  
}

#--------------------------------------LOOP1sp2----------------------------------------------


samples<-c(10:12, 28:30, 46:48, 64:66)

j=1

for (s in samples){
  
  spkF <- split(sp[[s]], species.filter)
  Species2 <- exprs(spkF$Species2)
  
  traits.max <- Species2
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(nrow(traits.max))
  print(ncol(traits.max))
  
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
  stockAS1[j, 1, 2]<-s
  stockAS1[j, 2, 2]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS1))) #Idem
  
}


#--------------------------------------FILTRE ASS2----------------------------------------------


## Here the kmeans function does not work. We have to delimite the gates by hand. 
# To delimite CHlamydomonas: 
##purify more the species
edges.ChlamAss2 <- matrix(c(2.5, 2.5, 4, 4, 2.8, 1.75, 1.75, 2.8), ncol = 2, nrow = 4)
colnames(edges.ChlamAss2) <- c("RED.B.HLin", "YEL.B.HLin")
ChlamAss2 <- polygonGate(filterId="ChlamAss2", .gate=edges.ChlamAss2)
## apply
ChlamAss2.filter <- !debris & singlets & ChlamAss2
ChlamAss2 <- Subset (log.fs, ChlamAss2.filter)
## Take of the unnecessary param
ChlamAss2 <- ChlamAss2 [, -9]
ChlamAss2 <- ChlamAss2 [, -10]
print(ChlamAss2)

# To delimite Rhodomonas: 
##purify more the species
edges.RhodAss2 <- matrix(c(3, 3, 4, 4, 3.75, 2.81, 2.81, 3.75), ncol = 2, nrow = 4)
colnames(edges.RhodAss2) <- c("RED.B.HLin", "YEL.B.HLin")
RhodAss2 <- polygonGate(filterId="RhodAss2", .gate=edges.RhodAss2)
## apply
RhodAss2.filter <- !debris & singlets & RhodAss2
RhodAss2 <- Subset (log.fs, RhodAss2.filter)
## Take of the unnecessary param
RhodAss2 <- RhodAss2 [, -9]
RhodAss2 <- RhodAss2 [, -10]
print(RhodAss2)


#--------------------------------------LOOP2Chlam2----------------------------------------------


samples<-c(13:15, 31:33, 49:51, 67:69)

j=1

for (s in samples){
  
  Species1 <- exprs(ChlamAss2[[s]])
  
  traits.max <- Species1
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(nrow(traits.max))
  print(ncol(traits.max))
  
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
  stockAS2[j, 1, 1]<-s
  stockAS2[j, 2, 1]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS2)))
  
}

#--------------------------------------LOOP2Rhod2----------------------------------------------


samples<-c(13:15, 31:33, 49:51, 67:69)

j=1

for (s in samples){
  
  Species2 <- exprs(RhodAss2[[s]])
  
  traits.max <- Species2
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(nrow(traits.max))
  print(ncol(traits.max))
  
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
  stockAS2[j, 1, 2]<-s
  stockAS2[j, 2, 2]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS2)))
  
}

#--------------------------------------FILTRE ASS3----------------------------------------------

## Here the kmeans function does not work. We have to delimite the gates by hand. 
# To delimite Chlorella: 
##purify more the species
edges.ChlorAss3 <- matrix(c(2, 2, 3.5, 3.5, 2.4, 0.5, 0.5, 2.4), ncol = 2, nrow = 4)
colnames(edges.ChlorAss3) <- c("RED.B.HLin", "YEL.B.HLin")
ChlorAss3 <- polygonGate(filterId="ChlorAss3", .gate=edges.ChlorAss3)
## apply
ChlorAss3.filter <- !debris & singlets & ChlorAss3
ChlorAss3 <- Subset (log.fs, ChlorAss3.filter)
## Take of the unnecessary param
ChlorAss3 <- ChlorAss3 [, -9]
ChlorAss3 <- ChlorAss3 [, -10]
print(ChlorAss3)

# To delimite Rhodomonas: 
##purify more the species
edges.RhodAss3 <- matrix(c(3, 3, 3.75, 3.75, 4, 2.5, 2.5, 4), ncol = 2, nrow = 4)
colnames(edges.RhodAss3) <- c("RED.B.HLin", "YEL.B.HLin")
RhodAss3 <- polygonGate(filterId="RhodAss3", .gate=edges.RhodAss3)
## apply
RhodAss3.filter <- !debris & singlets & RhodAss3 
RhodAss3  <- Subset (log.fs, RhodAss3.filter)
## Take of the unnecessary param
RhodAss3 <- RhodAss3 [, -9]
RhodAss3 <- RhodAss3 [, -10]
print(RhodAss3)

#--------------------------------------LOOP3Chlor3----------------------------------------------


samples<-c(16:18, 34:36, 52:54, 70:72)

j=1

for (s in samples){
  
  Species1 <- exprs(ChlorAss3[[s]])
  
  traits.max <- Species1
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(nrow(traits.max))
  print(ncol(traits.max))
  
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
  stockAS3[j, 1, 1]<-s
  stockAS3[j, 2, 1]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS3)))
  
}

#--------------------------------------LOOP3Rhod3----------------------------------------------


samples<-c(16:18, 34:36, 52:54, 70:72)
stockAS3<-matrix(NA, ncol=2, nrow=length(samples))

j=1

for (s in samples){
  
  Species2 <- exprs(RhodAss3[[s]])
  
  traits.max <- Species2
  
  # Define maximum number of points (max1) and number of traits under consideration (dim1)
  # Alternatively, it is possible to manually define max1 and dim1!!
  print(nrow(traits.max))
  print(ncol(traits.max))
  
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
  stockAS3[j, 1, 2]<-s
  stockAS3[j, 2, 2]<-TED.index(traits.max)
  
  j=j+1
  print(paste(j, "/", nrow(stockAS3)))
  
}



#---------------------------Saving------------------------
setwd("C:/Users/Elif Ka/OneDrive/Mémoire")

write.csv2(stockAS1[,,,1], "A1S1.csv")
write.csv2(stockAS1[,,,2], "A1S2.csv")

write.csv2(stockAS2[,,,1], "A2S1.csv")
write.csv2(stockAS2[,,,2], "A2S2.csv")

write.csv2(stockAS3[,,,1], "A3S1.csv")
write.csv2(stockAS3[,,,2], "A3S2.csv")