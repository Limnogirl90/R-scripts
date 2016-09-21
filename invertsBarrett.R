library(mvabund)

inverts = read.csv("invertsBarrett.csv",header=TRUE)
str(inverts)
X = inverts[,1:3]
X$Year = factor(X$Year)
# Analysing the ordinal (?) data, i.e. excluding Pulmonata/Prosobranchia
inverts = mvabund(inverts[,4:24])

# fit a main effects model
invertFit = manyglm(inverts~Habitat*Month*Year,data=X, offset = log(apply(inverts,1,mean)))
# Year entered into hte model last to look at effects of year after accounting for season and habitat
# an offset for log row sum included to (approximately) account for differences in sampling intensity across samples

# residual plot
plot(invertFit)
# looks good

# do an anova to get univariate results
# this will take about a minute to run.  For publication should actually use nBoot=999 and leave to run for 10 minutes.
an.inverts = anova(invertFit,p.uni="adjusted",nBoot=99)
an.inverts
# everything is significant, but it seems like there are really strong patterns here with year. Most (about 60%?) of the change in deviance is due to Year. I assume it is year that is of most interest.
# note from the univariate tests, that even after adjusting for multiple testing, there are significant year effects in ALL taxa. This is unusual.

# have more of a look at where the action is, from univariate test statistics:
an.snails$uni.test
# the big stats are mostly for year (fourth row)
sort(an.snails$uni.test[4,],decreasing=TRUE)
# the taxa that seem to have the strongest change over time are Elimia, Dreissena polymorpha, Acari, Gammarus faciatus
# although as before there are large and significant effects pretty much all over the place.

# plot data to look at year effect
plot.mvformula(log(inverts+1)~X$Year,pch=as.numeric(X$Habitat))
# this isn't the prettiest plot. Different colours for different years, different symbols for Habitat.
# you can see:
# - that there is a lot of variation across reps within site (some is due to habitat, e.g. strong effect in ..polymorpha, ...bugensis and ...gammarus)
# - many taxa have dropped to zero abundance at 2014, including the two most with strongest year effect (Elimia and D poly).
# - A few have started at zero in 1983 and moved up from there (including Acari)
# Might be worth using a line graph to visualise effects for each species?

# try a line plot of the means:
source("default.meanvar.plot.R") # to use the updated function that makes it easier to extract means
source("meanvar.plot.R")

mvYear=meanvar.plot(inverts~X$Year,table=TRUE)

# construct a line plot for all taxa and save as a png file
png(file="linePlot.png",res=500,width=3000,height=1500)
{
par(mar=c(3,3,0.5,0.5),mgp=c(1.75,0.75,0))
matplot(levels(X$Year),log(mvYear$mean+1),type="l",col=rainbow(21)[1:21],lty=1,xaxt="n",xlab="Year",ylab="Mean count",yaxt="n")
axis(1,c(levels(X$Year)))
yax = c(0,1,10,100,1000)
axis(2,at=log(yax+1),labels=yax)
legend("topright",colnames(mvYear$mean),col=rainbow(21)[1:21],lty=1,cex=0.5,y.intersp=0.75,ncol=3)
}
dev.off()

# could do something similar for year:Habitat or year:Month:Habitat if desired
mv3Int=meanvar.plot(inverts~interaction(X$Year,X$Month,X$Habitat),table=TRUE)
mvYearHab=meanvar.plot(inverts~interaction(X$Year,X$Habitat),table=TRUE)
# the idea there would probably be to do separate line plots for different habitats.



# haven't done an ordination - not sure if one is really needed here??? The effects are quite striking on plots of means already, you couldn't expect to learn much more or an ordination
