###########################################################
# 1. Load the necessary packages for the Cluster analysis #
###########################################################

library(car)
library(cluster)
library(corrplot)
library(fpc)
library(Hmisc)
library(HSAUR)
library(HSAUR2)
library(maptree)
library(mclust)
library(MVA)
library(lattice)
library(proxy)
library(psych)
library(reshape)
library(useful)

##################################
# 2. Load the survey data into R #
##################################

load("apphappyData.RData")
numdata <- apphappy.3.num.frame
labsdata <- apphappy.3.labs.frame

############################################
# 3. Create subsets for easier exploration #
############################################

labssub.demographic <- subset(labsdata, select = c("caseID","q1","q48","q49","q50r1","q50r2","q50r3",
                                                   "q50r4","q50r5","q54","q55","q56","q57"))

labssub.behavioral <- subset(labsdata, select = c("caseID","q2r1","q2r2","q2r3","q2r4","q2r5","q2r6",
                                                  "q2r7","q2r8","q2r9","q2r10","q4r1","q4r2",
                                                  "q4r3","q4r4","q4r5","q4r6","q4r7","q4r8",
                                                  "q4r9","q4r10","q4r11","q11","q12","q13r1",
                                                  "q13r1","q13r2","q13r3","q13r4","q13r5",
                                                  "q13r6","q13r7","q13r8","q13r9","q13r10",
                                                  "q13r11","q13r12"))

numsub.attitudinal <- subset(numdata, select=c("caseID","q24r1","q24r2","q24r3","q24r4","q24r5","q24r6",
                                               "q24r7","q24r8","q24r9","q24r10","q24r11","q24r12",
                                               "q25r1","q25r2","q25r3","q25r4","q25r5","q25r6",
                                               "q25r7","q25r8","q25r9","q25r10","q25r11","q25r12",
                                               "q26r3","q26r4","q26r5","q26r6","q26r7","q26r8",
                                               "q26r9","q26r10","q26r11","q26r12","q26r13","q26r14",
                                               "q26r15","q26r16","q26r17","q26r18"))

numsub <- numsub.attitudinal

### Get Rid of Cases with all the Same Attitudinal Responses
attach(numsub)
numsub.drop <- numsub
numsub.drop$mean <- (q24r1+q24r2+q24r3+q24r4+q24r5+q24r6+
                         q24r7+q24r8+q24r9+q24r10+q24r11+q24r12+
                         q25r1+q25r2+q25r3+q25r4+q25r5+q25r6+
                         q25r7+q25r8+q25r9+q25r10+q25r11+q25r12+
                         q26r3+q26r4+q26r5+q26r6+q26r7+q26r8+q26r9+
                         q26r10+q26r11+q26r12+q26r13+q26r14+q26r15+
                         q26r16+q26r17+q26r18)/40

numsub.drop$var <- 0
for (i in 1:nrow(numsub.drop)) {
    numsub.drop$var[i] <- var(as.numeric(numsub.drop[i,2:41]))
}
numsub.drop <- numsub.drop[numsub.drop$var != 0,1:41]

labssub.behavioral.drop <- labssub.behavioral[labssub.behavioral$caseID %in% numsub.drop$caseID ,]
labssub.demographic.drop <- labssub.demographic[labssub.demographic$caseID %in% numsub.drop$caseID ,]
numdata.drop <- numdata[numdata$caseID %in% numsub.drop$caseID ,]

##################################
# 4. Exploratory data analysis 1 #
##################################

str(numdata)
head(numdata)
tail(numdata)
summary(numdata)
# caseId is non-numeric
# NAs should be of concern for q5r1, q12, q57

str(numsub.attitudinal)
head(numsub.attitudinal)
tail(numsub.attitudinal)
summary(numsub.attitudinal)
# no missing values; all between 1 and 6, as expected

str(numsub.drop)
head(numsub.drop)
tail(numsub.drop)
summary(numsub.drop)

## Univariate Plots for Attitudinal Variables
par(mfrow = c(4,5))
for (i in 2:21) {
    hist(numsub.drop[,i], xlab = names(numsub.drop[i]), 
         main = paste("Histogram", names(numsub.drop[i])), ylim = c(0,1000), col="grey")}
for (i in 22:41) {
    hist(numsub.drop[,i], xlab = names(numsub.drop[i]), 
         main = paste("Histogram", names(numsub.drop[i])), ylim = c(0,1000), col="grey")}

## Univatriate Behavioral Variable Bar Plots
par(mfrow = c(3,3))
for (i in 2:10) {
    barplot(table(labssub.behavioral.drop[,i]), xlab = names(labssub.behavioral.drop[i]), 
            main = paste("Histogram", names(labssub.behavioral.drop[i])), ylim = c(0,1800))}
for (i in 11:19) {
    barplot(table(labssub.behavioral.drop[,i]), xlab = names(labssub.behavioral.drop[i]), 
            main = paste("Histogram", names(labssub.behavioral.drop[i])), ylim = c(0,1800))}
for (i in 20:28) {
    barplot(table(labssub.behavioral.drop[,i]), xlab = names(labssub.behavioral.drop[i]), 
            main = paste("Histogram", names(labssub.behavioral.drop[i])), ylim = c(0,1800))}
for (i in 29:37) {
    barplot(table(labssub.behavioral.drop[,i]), xlab = names(labssub.behavioral.drop[i]), 
            main = paste("Histogram", names(labssub.behavioral.drop[i])), ylim = c(0,1800))}

## Univariate Demographic Variable Bar Plots
par(mfrow = c(3,4))
for (i in 2:13) {
    barplot(table(labssub.demographic.drop[,i]), xlab = names(labssub.demographic.drop[i]), 
            main = paste("Histogram", names(labssub.demographic.drop[i])), ylim = c(0,1800))}

## Correlation Plot for Attitudinal Variables
par(mfrow = c(1,1))
rcorr(as.matrix(numsub.drop[,2:41]), type="pearson")
numsub.drop.correlation <- cor(numsub.drop[,2:41])
corrplot(numsub.drop.correlation, method="shade", shade.col=NA, tl.col="black", cex.lab = 0.25)

## Factor Analysis for Variable Reduction
fit <- principal(numsub.attitudinal[,2:41], nfactors=12, rotate="varimax")
summary(fit)
fit # print results

#######################
# 5. Transformation 1 #
#######################

#####################################
## DR. SRINIVASAN'S BASIS VARIABLES #
#####################################

attach(numsub.drop)

numsub.drop$tech.positive <- (q24r1+q24r2+q24r3+q24r5+q24r6)/5 ## POSITIVE TO TECH
numsub.drop$media.positive <- (q24r7+q24r8)/2 ## POSTIVE TO MEDIA
numsub.drop$internet.comms <- (q24r10+q24r11)/2 ## INTERNET/COMMS
numsub.drop$negative.to.tech <- (q24r4+q24r9+q24r12)/3 ## NEGATIVE TO TECH

numsub.drop$opinion.leader <- (q25r1+q25r2+q25r3+q25r4+q25r5)/5 ## OPINION LEADER
numsub.drop$control <- (q25r7+q25r8)/2 ## CONTROL
numsub.drop$drive <- (q25r9+q25r10+q25r11)/3 ## DRIVE
numsub.drop$negative <- (q25r6+q25r12)/2 ##NEGATIVE

numsub.drop$bargain <- (q26r3+q26r4+q26r5+q26r6+q26r7)/5 ## BARGAIN SEEKER
numsub.drop$show.off <- (q26r8+q26r9+q26r10)/3 ## SHOW OFF
numsub.drop$children <- q26r11 ## CHILDREN
numsub.drop$hot <- (q26r12+q26r13+q26r14)/3 ## HOT
numsub.drop$brand <- (q26r15+q26r16+q26r17+q26r18)/4 ## BRAND

numsub2a <- subset(numsub.drop, select=c("tech.positive","media.positive","internet.comms",
                                        "opinion.leader","control","drive",
                                        "bargain","show.off","children","hot","brand"))

numsub2b <- subset(numsub.drop, select=c("tech.positive","media.positive","internet.comms",
                                         "opinion.leader","control","drive",
                                         "bargain","show.off","children","hot","brand",
                                         "negative.to.tech","negative"))

## Univariate Plots for Attitudinal Variables after Variable Reduction
par(mfrow = c(4,4))
for (i in 1:13) {
    hist(numsub2b[,i], xlab = names(numsub2b[i]), 
         main = paste("Histogram", names(numsub2b[i])), ylim = c(0,1000), col="grey")}

## Correlation Plot for Transformed Attitudinal Variables
par(mfrow = c(1,1))
rcorr(as.matrix(numsub2b), type="pearson")
numsub.attitudinal.correlation <- cor(numsub2b)
corrplot(numsub.attitudinal.correlation, method="shade", shade.col=NA, tl.col="black", cex.lab = 0.25)

###############################
## CLASSMATE'S BASIS APPROACH #
###############################

attach(numsub.drop)

numsub.drop$active.influencer <- (q25r1+q25r3+q25r4+q25r5+q25r10)/5 ## ACTIVE INFLUENCER
numsub.drop$active.physical <-q25r11 ## ACTIVE PHYSICAL
numsub.drop$active.shopper <-(q24r3+q26r3+q26r4+q26r5+q26r6+q26r7+q26r12+q26r13+q26r16+q26r18)/10 ## ACTIVE SHOPPER
numsub.drop$not.tech.savvy <-(q24r3+q24r9)/2 ## TECH NOT SAVVY
numsub.drop$social.media <-(q24r10+q24r11+q26r14)/3 ## ACTIVE SOCIAL MEDIA
numsub.drop$active.style <-q26r15 ## ACTIVE STYLE
numsub.drop$considers.self <- (q25r2+q25r8+q25r9+q25r10+q25r12)/5 ## CONSIDERS SELF
numsub.drop$control2 <-(q24r12+q24r7)/2 ## CONTROL
numsub.drop$passive.influencer <-q24r2 ## PASSIVE INFLUENCER
numsub.drop$tech.savvy <-(q24r1+q24r5+q24r6+q24r8+q26r8)/5 ## TECH SAVVY

numsub3 <- subset(numsub.drop, select=c("active.influencer","active.physical","active.shopper",
                                        "social.media","active.style","considers.self",
                                        "control2","passive.influencer","tech.savvy"))

## Univariate Plots for Attitudinal Variables after Variable Reduction
par(mfrow = c(3,4))
for (i in 1:10) {
    hist(numsub3[,i], xlab = names(numsub3[i]), 
         main = paste("Histogram", names(numsub3[i])), ylim = c(0,1000), col="grey")}

## Correlation Plot for Transformed Attitudinal Variables
par(mfrow = c(1,1))
rcorr(as.matrix(numsub3), type="pearson")
numsub.attitudinal.correlation <- cor(numsub3)
corrplot(numsub.attitudinal.correlation, method="shade", shade.col=NA, tl.col="black", cex.lab = 0.25)

##############
## MY DESIGN #
##############
attach(numsub.drop)
numsub.drop$shopper <- (q26r4+q26r6+q26r16+q26r12+q26r13+q26r8)/6 ## SHOPPER	
numsub.drop$active <- q24r3 ## ACTIVE
numsub.drop$advisor <- (q25r1+q25r3+q24r2)/3 ## ADVISOR
numsub.drop$control.oriented <- (q25r7+q25r4+q24r12+q24r5)/5 ## CONTROL ORIENTED
numsub.drop$cost.oriented <- q26r3 ## COST ORIENTED
numsub.drop$driven <- (q25r10+q25r11)/2 ## DRIVEN
numsub.drop$entertainment.oriented <- (q26r17+q24r7+q24r8)/3 ## ENTERTAINMENT ORIENTED	
numsub.drop$family.oriented <- q26r11 ## FAMILY ORIENTED	
numsub.drop$image.oriented <- (q26r14+q26r9+q26r10+q26r7+q26r15+q26r18)/6 ## IMAGE ORIENTED
numsub.drop$socially.oriented <- (q24r10+q24r11)/2 ## SOCIALLY ORIENTED
numsub.drop$tech.averse <- (q24r4+q24r9)/2 ## TECH AVERSE
numsub.drop$tech.follower <- q24r1 ## TECH FOLLOWER
numsub.drop$time.oriented <- (q26r5+q25r12+q24r6)/3 ## TIME ORIENTED
numsub.drop$trail.blazer <- (q25r8+q25r9+q25r2+q25r5)/4 ## TRAIL BLAZER

numsub4 <- subset(numsub.drop, select=c("shopper","active","advisor","control.oriented",
                                        "cost.oriented","driven","entertainment.oriented",
                                        "family.oriented","image.oriented","socially.oriented",
                                        "tech.follower","time.oriented","trail.blazer"))

## Univariate Plots for Attitudinal Variables after Variable Reduction
par(mfrow = c(4,4))
for (i in 1:15) {
    hist(numsub4[,i], xlab = names(numsub4[i]), 
         main = paste("Histogram", names(numsub4[i])), ylim = c(0,1000), col="grey")}

## Correlation Plot for Transformed Attitudinal Variables
par(mfrow = c(1,1))
rcorr(as.matrix(numsub4), type="pearson")
numsub.attitudinal.correlation <- cor(numsub4)
corrplot(numsub.attitudinal.correlation, method="shade", shade.col=NA, tl.col="black", cex.lab = 0.25)

#####################################
## DR. SRINIVASAN'S BASIS VARIABLES #
#####################################

attach(numsub.drop)

numsub.drop$media.positive2 <- (q24r7+q24r8+q26r17)/3 ## POSTIVE TO MEDIA
numsub.drop$bargain2 <- (q26r3+q26r5)/2 ## BARGAIN SEEKER
numsub.drop$app.fan <- (q26r8+q26r9+q26r10+q26r12)/4 ## APP FAN
numsub.drop$spender <- (q26r4+q26r6+q26r13+q26r16)/4 ## SPENDER
numsub.drop$brand2 <- (q26r7+q26r14+q26r15+q26r18)/4 ## BRAND

numsub5 <- subset(numsub.drop, select=c("tech.positive","media.positive2","internet.comms",
                                         "opinion.leader","control","drive",
                                         "bargain2","app.fan","children","spender","brand2"))
summary(numsub5)

## Univariate Plots for Attitudinal Variables after Variable Reduction
par(mfrow = c(3,4))
for (i in 1:11) {
    hist(numsub5[,i], xlab = names(numsub5[i]), 
         main = paste("Histogram", names(numsub5[i])), ylim = c(0,1000), col="grey")}

## Correlation Plot for Transformed Attitudinal Variables
par(mfrow = c(1,1))
rcorr(as.matrix(numsub5), type="pearson")
numsub.attitudinal.correlation <- cor(numsub5)
corrplot(numsub.attitudinal.correlation, method="shade", shade.col=NA, tl.col="black", cex.lab = 0.25)

## PCA
pca <- princomp(numsub5)
summary(pca)
par(mfrow=c(1,1))
plot(pca$scores[,1],pca$scores[,2]) 
plot(pca$scores[,1],pca$scores[,2]) 

loadings(pca)

## Outliers with all 1's
## sort(pca$scores[,1])
## numsub.attitudinal["431",]           
## etc. 


##################################
# 5. K Means Clustering Approach #
##################################

# Create a 'scree' plot to determine the num of clusters

wssplot <- function(numsub, nc=20, seed=1234) {
  wss <- (nrow(numsub)-1)*sum(apply(numsub,2,var))
  for (i in 2:nc) {
    set.seed(seed)
    wss[i] <- sum(kmeans(numsub, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")} 

par(mfrow=c(1,1))
wssplot(numsub5) # Scree plot suggests 4, 5, or 6 clusters would be ideal

############
############
############

clusterresults <- kmeans(numsub5,5)
##names(clusterresults)
##clusterresults$withinss
##clusterresults$tot.withinss
##clusterresults$totss
##clusterresults$betweenss
clusterresults$size
rsquare <- clusterresults$betweenss/clusterresults$totss
rsquare

# Create a PC (Principal Component plot)

plot(clusterresults, data=numsub5)
clusterresults$centers
head(clusterresults$cluster)

dev.off()
dissE <- daisy(numsub5)
## names(dissE)
dE2   <- dissE^2
sk2   <- silhouette(clusterresults$cluster, dE2)
str(sk2)
plot(sk2)

###############################################
# Create new data set with Cluster Assignment #
###############################################

newdf <- as.data.frame(clusterresults$cluster)
write.csv(newdf, file = "clusterresults.csv")
write.csv(numsub.attitudinal, file = "numsub.csv")

newdf <- read.csv("clusterresults.csv")
combdata <- cbind(numsub.drop,newdf,numsub.drop$caseID)
combdata2 <- cbind(numdata.drop,newdf,numdata.drop$caseID)
combdata3 <- cbind(labssub.behavioral.drop,newdf,labssub.behavioral.drop$caseID)
combdata4 <- cbind(labssub.demographic.drop,newdf,labssub.demographic.drop$caseID)

combdata <- rename(combdata, c(clusterresults.cluster="cluster"))
combdata2 <- rename(combdata2, c(clusterresults.cluster="cluster"))
combdata3 <- rename(combdata3, c(clusterresults.cluster="cluster"))
combdata4 <- rename(combdata4, c(clusterresults.cluster="cluster"))

attach(combdata2)
combdata2$device.count <- q2r1+q2r2+q2r3+q2r4+q2r5+q2r6+q2r7+q2r8+q2r9
combdata2$child.at.home <- q50r2+q50r3+q50r4
detach(combdata2)

aggregate(combdata,by=list(byvar=combdata$cluster), mean)
aggregate(combdata2,by=list(byvar=combdata$cluster), mean)



par(mfrow=c(4,4))
for (i in 2:17) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
    }

for (i in 18:33) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
}

for (i in 34:41) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
}

for (i in 42:57) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
}

for (i in 58:73) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
}

for (i in 74:89) {
    barplot(prop.table(table(combdata2[,i], combdata2$cluster), 2),
            main=names(combdata2[i]),
            xlab="Segment")
}

par(mfrow=c(4,4))
for (i in 42:57) {
    boxplot(combdata[,i]~combdata$cluster,
            main=names(combdata[i]),
            ylab="Segment",horizontal=TRUE)
}

for (i in 58:73) {
    boxplot(combdata[,i]~combdata$cluster,
            main=names(combdata[i]),
            ylab="Segment",horizontal=TRUE)
}

for (i in 74:80) {
    boxplot(combdata[,i]~combdata2$cluster,
            main=names(combdata[i]),
            ylab="Segment",horizontal=TRUE)
}

###########################
# Hierarchical Clustering #
###########################

numsub.dist = dist(numsub5)

hclustmodel <- hclust(dist(numsub5), method = "complete")
plot(hclustmodel)

cut.5 <- cutree(hclustmodel, k=5)
plot(silhouette(cut.5,numsub.dist))
head(cut.5)

write.csv(cut.5, file = "cut5results.csv")
########################################
# For hclust how to calculate BSS & TSS
######################################
numsubmat <- as.matrix(numsub5)
overallmean <- matrix(apply(numsubmat,2,mean),nrow=1)
overallmean
TSS <- sum(dist(numsubmat,overallmean)^2)
TSS

###################################
# Compute WSS based on 5 clusters #
###################################
combcutdata <- cbind(numsub5,cut.5)
head(combcutdata)

combcutdata <- rename(combcutdata, c(cut.5="cluster"))
head(combcutdata)

clust1 <- subset(combcutdata, cluster == 1)
clust1 <- subset(clust1, select=c("tech.positive","media.positive2","internet.comms",
                                  "opinion.leader","control","drive",
                                  "bargain2","app.fan","children","spender","brand2"))
clust1 <- as.matrix(clust1,rowby=T)
clust1mean <- matrix(apply(clust1,2,mean),nrow=1)
dis1 <- sum(dist(clust1mean,clust1)^2)

clust2 <- subset(combcutdata, cluster == 2)
clust2 <- subset(clust2, select=c("tech.positive","media.positive2","internet.comms",
                                  "opinion.leader","control","drive",
                                  "bargain2","app.fan","children","spender","brand2"))
clust2 <- as.matrix(clust2,rowby=T)
clust2mean <- matrix(apply(clust2,2,mean),nrow=1)
dis2 <- sum(dist(clust2mean,clust2)^2)

clust3 <- subset(combcutdata, cluster == 3)
clust3 <- subset(clust3, select=c("tech.positive","media.positive2","internet.comms",
                                  "opinion.leader","control","drive",
                                  "bargain2","app.fan","children","spender","brand2"))
clust3 <- as.matrix(clust3,rowby=T)
clust3mean <- matrix(apply(clust3,2,mean),nrow=1)
dis3 <- sum(dist(clust3mean,clust3)^2)

clust4 <- subset(combcutdata, cluster == 4)
clust4 <- subset(clust4, select=c("tech.positive","media.positive2","internet.comms",
                                  "opinion.leader","control","drive",
                                  "bargain2","app.fan","children","spender","brand2"))
clust4 <- as.matrix(clust4,rowby=T)
clust4mean <- matrix(apply(clust4,2,mean),nrow=1)
dis4 <- sum(dist(clust4mean,clust4)^2)

clust5 <- subset(combcutdata, cluster == 5)
clust5 <- subset(clust5, select=c("tech.positive","media.positive2","internet.comms",
                                  "opinion.leader","control","drive",
                                  "bargain2","app.fan","children","spender","brand2"))
clust5 <- as.matrix(clust5,rowby=T)
clust5mean <- matrix(apply(clust5,2,mean),nrow=1)
dis5 <- sum(dist(clust5mean,clust5)^2)

WSS <- sum(dis1,dis2,dis3,dis4,dis5)
WSS

BSS <- TSS - WSS
BSS
## calculating the % of Between SS/ Total SS
rsquare <- BSS/TSS
rsquare

#######################################################
### A little function to calculate the average silhouette width
### for a variety of choices of k for PAM method:
###########################################################
my.k.choices <- 2:8
avg.sil.width <- rep(0, times=length(my.k.choices))
for (ii in (1:length(my.k.choices)) ){
  avg.sil.width[ii] <- pam(numsub, k=my.k.choices[ii])$silinfo$avg.width
}
print( cbind(my.k.choices, avg.sil.width) )

#################
# 9. PAM method #
#################
clusterresultsPAM <-pam(numsub5,5)
summary(clusterresultsPAM)
plot(clusterresultsPAM, which.plots=1)
plot(clusterresultsPAM, which.plots=2)

##############################
# 10. Model based clustering #
##############################

fit <- Mclust(numsub5, 5)
plot(fit, data=numsub5, what="density") # plot results
summary(fit) # display the best model
clusplot(numsub.drop, mclust5$class)
