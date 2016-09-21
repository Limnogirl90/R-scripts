###Some Ordination Methods

rm(list=ls())
bumpus<-read.csv("Bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of linear measurements
bumpus.data<-bumpus.data[,-3]  #remove wt

sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))

########################
#	PCA: Principal Components Analysis

#### PCA "by hand" #####
vcv.bumpus<-var(bumpus.data)	#Calculate PC axes
pc.bumpus<-eigen(vcv.bumpus) 
pc.eval<-pc.bumpus$values
pc.evec<-pc.bumpus$vectors
sum(pc.eval)
sum(diag(vcv.bumpus))	#Same answer. Check!
bumpus.centered<-as.matrix(scale(bumpus.data,scale=F))    #mean-center
pc.scores<-bumpus.centered%*%pc.evec	#Projection
plot(pc.scores,xlab="PC I", ylab="PC II",asp=1)

###########
### PCA easier way
pca.bumpus<-prcomp(bumpus.data)	#uses SVD for computations (mean-centered is default)
summary(pca.bumpus)
PC.scores<-pca.bumpus$x #This gives PCA scores
plot(PC.scores,xlab="PC I", ylab="PC II",asp=1)

#########PC plot by group:  make color vector for groups to use
colo<-rep("blue",nrow(bumpus.data)); colo[which(sex == 'f')]<-"red"
   plot(PC.scores,pch=21,bg=colo,cex=1.5,asp=1)

## PC plot by sex*surv groups
symb<-rep(21,nrow(bumpus.data));symb[which(surv == 'FALSE')]<-24
   plot(PC.scores,pch=symb,bg=colo,cex=1.5,asp=1) #cex is size of symbols


#############Biplot with PCA
biplot(pca.bumpus)

###########################
# 	PCoA
library(vegan)
library(MASS)
bumpus.dist<-dist(bumpus.data) #make distance matrix use what best for your data
PCoA<-cmdscale(bumpus.dist)
plot(PCoA,pch=symb,bg=colo,col="black",cex=1.5,asp=1)

#############################
#	NMDS
bumpus.dist<-dist(bumpus.data)
bumpus.nmds<-isoMDS(bumpus.dist, k = 2) #NOTE: metaMDS can also be used
ordiplot(bumpus.nmds)

bumpus.new<-as.matrix(bumpus.nmds$points)
plot(bumpus.new,pch=symb,bg=colo,col="black",cex=1.5,asp=1)

###########################
### CVA
lda.bumpus<-lda(bumpus.data,SexBySurv) #lda is linear discriment analysis
cva.bumpus<-predict(lda.bumpus,bumpus.data)
cv.scores<-cva.bumpus$x
plot(cv.scores,pch=symb,bg=colo,col="black",cex=1.5,asp=1)
	
