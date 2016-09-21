# Multivariate exploration and MANOVA

rm(list=ls())
bumpus<-read.csv("bumpus.csv",header=T)
bumpus.data<-as.matrix(bumpus[,(5:13)]) # matrix of linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))
 bumpus.data<-log(bumpus.data)
Y<-bumpus.data[,-1] # Tells to take complete data minus first column
TotalLength<-bumpus.data[,1]

colo<-rep("blue",nrow(bumpus.data)); colo[which(SexBySurv == 'f TRUE')]<-"red";
	colo[which(SexBySurv == 'f FALSE')]<-"pink"; colo[which(SexBySurv == 'm FALSE')]<-"lightblue";

#__________________________________________________________________________#
#Describe your data
var(bumpus.data)
  (100*apply(bumpus.data,2,sd)) /apply(bumpus.data,2,mean)  # CV: Coefficient of Variation

cor.bumpus<-cor(bumpus.data)
cor.bumpus
pairs(bumpus.data)
vcv.bumpus<-var(bumpus.data)
vcv.bumpus
var(scale(bumpus.data))
dist(bumpus.data)

#__________________________________________________________________________#
#single factor MANOVA
model1<-lm(bumpus.data~sex)
summary(model1)	#yields a set of univariate analyses

summary(manova(model1))	#does multivariate test
	#using Wilks' lambda as test-statistic
summary(manova(model1),test="Wilks")	#does multivariate test

#__________________________________________________________________________#
#Factorial MANOVA
model2<-lm(bumpus.data~sex*surv)
summary(manova(model2))

# VARIANCE COMPONENTS
  summary(manova(model2))$SS[[1]]	#SSCP matrices for each term in model
  summary(manova(model2))$SS[[2]]
  summary(manova(model2))$SS[[3]]
  summary(manova(model2))$SS[[4]]

###  Advanced note: One can also obtain SSCP via matrix computations in lecture.
	# Recall that the SSCP of a model is really its SSCP.error
  
#__________________________________________________________________________#
#Pairwise comparisons based on D_Euclidean among LS.means

yhat.full<-predict(manova(model2))
lsmeans.obs <- rowsum(yhat.full, SexBySurv)/as.vector(table(SexBySurv))    #group means
plot(lsmeans.obs[,1],lsmeans.obs[,2])

D.obs<-as.matrix(dist(lsmeans.obs)) # Euclidean D between pairs of groups
D.obs

permute<-999
PDist<-array(1, dim=c(dim(lsmeans.obs)[1],dim(lsmeans.obs)[1]))
for(i in 1:permute){
	bumpus.rand<-bumpus.data[sample(nrow(bumpus.data)),]		#shuffle specimens 
  	yhat.rand<-predict(manova(bumpus.rand~sex*surv))	
  	lsmeans.rand <- rowsum(yhat.rand, SexBySurv)/as.vector(table(SexBySurv))    #group means/Shuffled
	D.rand<-as.matrix(dist(lsmeans.rand))		#Euclidean distance
	PDist<-ifelse(D.rand>=D.obs, PDist+1, PDist)
}  #end permute

PDist<-PDist/(permute+1)

D.obs
PDist

#__________________________________________________________________________#
### Multivariate Regression
summary(manova(Y~TotalLength))

### Visualizing multivariate regression with summary variables
# 1: Visualize Regression PC1
PC1<-prcomp(Y)$x[,1]
plot(TotalLength, PC1, pch=21, cex=1.25, bg="black")

# 2: Visualize with Regression scores: Drake and Klingenberg 2008
B<-(coef(lm(Y~TotalLength)))
s<-Y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) #slope portion only
plot(TotalLength, s, pch=21, cex=1.25, bg="black")

# 3: Stylized Regression lines: Adams and Nistri 2010
yhat.mancova<-predict(lm(Y~TotalLength))
P<-prcomp(yhat.mancova)$x[,1] #1st dim of predicteds: Adams and Nistri 2010
##NOTE: Adams and Nistri had multiple groups and used: shape~group*size)
plot(TotalLength, P, pch=21, cex=1.5, bg="black")

#__________________________________________________________________________#
#MANCOVA	RUN with Totlength as covariate
summary(manova(lm(Y~TotalLength*sex*surv)))
  # FIT COMMON SLOPE
summary(manova(lm(Y~TotalLength+sex*surv)))

## COMPARE PREDICTED VALUES
  yhat.full2<-predict(manova(lm(Y~TotalLength*sex+surv)))
lsmeans.obs2 <- rowsum(yhat.full2, SexBySurv)/as.vector(table(SexBySurv))    #group means
lsmeans.obs2
lsmeans.obs[,-1]  #ls means from before  (note the values change due to including covariate)

### Visualizing MANCOVA
# 1: Visualize Regression PC1
PC1<-prcomp(Y)$x[,1]
plot(TotalLength, PC1, pch=21, cex=1.25, bg=colo)

# 2: Visualize with Regression scores: Drake and Klingenberg 2008
B<-(coef(lm(Y~TotalLength+SexBySurv)))
s<-Y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) #slope portion only
plot(TotalLength, s, pch=21, cex=1.25, bg=colo)

# 3: Stylized Regression lines: Adams and Nistri 2010
yhat.mancova<-predict(lm(Y~TotalLength+SexBySurv))
P<-prcomp(yhat.mancova)$x[,1] #1st dim of predicteds: Adams and Nistri 2010
##NOTE: Adams and Nistri had multiple groups and used: shape~group*size)
plot(TotalLength, P, pch=21, cex=1.5, bg=colo)

#__________________________________________________________________________#
### Distances from data
	#There are MANY similarity and distance coefficients for different data types. 
	#Must use appropriate measure/metric for your data.

#Some useful functions are:  'dist' [base], 'vegdist' [vegan]

#__________________________________________________________________________#
### PERMUTATIONAL MANOVA: based on distance matrices
#install.packages("vegan")
library(vegan)
summary(manova(lm(bumpus.data~sex*surv)))
bumpus.dist<-dist(bumpus.data)   #generate distance matrix
adonis(bumpus.dist~sex*surv)


