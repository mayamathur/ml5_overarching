###########Many Labs 5 Overarching Analysis###############
#####Charlie Ebersole#####
######Begin February 10, 2020######

###Many Labs 5 will investigate whether replications that have undergone formal expert review produce effects more consistent with original studies compared to replications that have not undergone extra review. We are investigating 10 previous effects that were all replicated during the Reproducibility Project: Psychology (RP:P; OSC, 2015). For each effect, we are conducting two sets of replications: one using the methods used in RP:P and a second using a protocol revised through expert review. 
###We investigate our primary question in two ways.

####Approach #1####
#Reading data
#setwd("C:/Users/Charlie/Dropbox/Charlie_Harddrive/Desktop/Ongoing Projects/Many Labs 5/Overarching Paper/Results Included")
setwd("~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie")
Data<-read.csv(file="Summary Effect Sizes - Single Effects.csv",header=TRUE)
head(Data)
str(Data)

###Explanation of variables###
#Study - Which effect is being replicated? 
#Site - Where were the data collected? 
#Version - Is the effect size from the RP:P or Revised protocol?
#ES.Metric - what is the effect size for this estimate?
#Value - the value of the effect size, with more positive numbers indicating stronger effects in the same direction as the original study
#N - The sample size for that estimate
#Reverse - Does the effect size need to be reversed? (for variance explained effect sizes)
#Label - combination of Study and Version for labeling figures

require(metafor)
require(psych)

#Convert Effect Sizes to a common metric (correlation coefficient)

#Variance explained measures
list(levels(Data$ES.Metric))
R2<-subset(Data,ES.Metric == "pseudo R2" | ES.Metric == "eta sq")
R2$Value<-sqrt(R2$Value)

#d or g
D<-subset(Data,ES.Metric == "d" | ES.Metric == "g")
D$Value<-d2r(D$Value)

#r or partial r
R<-subset(Data,ES.Metric == "r" | ES.Metric == "rp")

###Merging data back together
Data<-rbind(R2,D,R)

#Reversing Effect Sizes where revised version showed weaker/reversal effect

NoReverse<-subset(Data,Reverse == "no")
Reverse<-subset(Data,Reverse == "yes")
Reverse$Value<-Reverse$Value*-1

###Merging back together
Data<-rbind(NoReverse,Reverse)


####Analyses####
ModelStat<-escalc(ri=Value,ni=N,data=Data,measure="COR",vtype="LS")
ModelStat
str(ModelStat$Version)
ModelStat$Version<-relevel(ModelStat$Version,"RP:P")

#We will begin with random intercepts of Site nested within Study. If the model fails to converge, we will drop the random intercept of Site. 
#Our main outcome of interest is testing for moderation by Version. That is, are effects from the revised protocols reliably different from those using the RP:P methods? If so, how large is the difference?

# MM: *** MAIN APPROACH #1 MODEL
ModelMV<-rma.mv(yi,vi,mods = ~ Version,random = ~ 1|Study/Site,data=ModelStat)
summary(ModelMV)
forest.rma(ModelMV,xlab="Focal Replication Effect (Pearson's r)",slab=ModelStat$Label,cex=.6)

#Breaking apart protocols for visualization
RPP1<-subset(ModelStat,Version=="RP:P")
Revised1<-subset(ModelStat,Version=="Revised")

ModelMVRPP<-rma.mv(yi,vi,random = ~ 1|Study/Site,data=RPP1)
summary(ModelMVRPP)
forest.rma(ModelMVRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP1$Label,cex=.6)

ModelMVRev<-rma.mv(yi,vi,random = ~ 1|Study/Site,data=Revised1)
summary(ModelMVRev)
forest.rma(ModelMVRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised1$Label,cex=.6)


####Approach 2: Analyzing effect of Version from each Study####
#For this approach, we will estimate the effect of Version within each Study (e.g., how much stronger is the replication effect in the Revised Version compared to the RP:P Version). We will then meta-analyze those estimates in order to estimate an overall effect of Version across all Studies. More positive effect sizes are indicative of stronger effects from the Revised Version compared to the RP:P Version.

#Reading in the data
Data2<-read.csv(file="Summary Effect Sizes - Interaction Effects.csv",header=TRUE)
Data2
str(Data2)

#Convert Effect Sizes to a common metric (correlation coefficient)

#Variance explained measures
list(levels(Data2$ES.Metric))
R2<-subset(Data2,ES.Metric == "pseudo R2" | ES.Metric == "partial eta2")
R2$Value<-sqrt(R2$Value)

#g
G<-subset(Data2,ES.Metric == "g")
G$Value<-d2r(G$Value)

###Merging data back together
Data2<-rbind(R2,G)

#Reversing Effect Sizes where revised version showed weaker/reversal effect

NoReverse<-subset(Data2,Reverse == "no")
Reverse<-subset(Data2,Reverse == "yes")
Reverse$Value<-Reverse$Value*-1

###Merging back together
Data2<-rbind(NoReverse,Reverse)




####Analyses####
ModelStat2<-escalc(ri=Value,ni=N,data=Data2,measure="COR",vtype="LS")
ModelStat2

# MM: *** MAIN ANALYSIS MODEL #2
Model2<-rma.uni(yi,vi,data=ModelStat2, test = "knha")
summary(Model2)
# compare to the LMM version:
summary(ModelMV)

forest.rma(Model2,xlab="Effect of Protocol",slab=ModelStat2$Study)


####Breaking apart by study and version for summary purposes####
head(ModelStat)
list(levels(ModelStat$Study))




###Albarracin S5
StudyData<-subset(ModelStat, Study == "Albarracin S5")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

###Albarracin S7
StudyData<-subset(ModelStat, Study == "Albarracin S7")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

###Crosby
StudyData<-subset(ModelStat, Study == "Crosby")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#Forster
StudyData<-subset(ModelStat, Study == "Forster")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#LoBue
StudyData<-subset(ModelStat, Study == "LoBue")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#Payne
StudyData<-subset(ModelStat, Study == "Payne")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#Risen
StudyData<-subset(ModelStat, Study == "Risen")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#Shnabel
StudyData<-subset(ModelStat, Study == "Shnabel")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#van Dijk
StudyData<-subset(ModelStat, Study == "van Dijk")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

#Vohs
StudyData<-subset(ModelStat, Study == "Vohs")

RPP<-subset(StudyData, Version == "RP:P")
StudyMetaRPP<-rma(yi,vi,data=RPP)
summary(StudyMetaRPP)
forest.rma(StudyMetaRPP,xlab="Focal Replication Effect (Pearson's r)",slab=RPP$Site,cex=.6)

Revised<-subset(StudyData, Version == "Revised")
StudyMetaRev<-rma(yi,vi,data=Revised)
summary(StudyMetaRev)
forest.rma(StudyMetaRev,xlab="Focal Replication Effect (Pearson's r)",slab=Revised$Site,cex=.6)

####Summary Statistics####
require(doBy)
head(Data)
str(Data)
summaryBy(Value~Version, FUN = c(mean,sd,median), data = Data)
sum(Data[Data$Study == "Albarracin 5"])
length(Data$Value)

###Calculating Average Sample Sizes for Studies and Versions###
Alb5RPP<-subset(Data, Study == "Albarracin S5" & Version == "RP:P")
Alb7RPP<-subset(Data, Study == "Albarracin S7" & Version == "RP:P")
CrosbyRPP<-subset(Data, Study == "Crosby" & Version == "RP:P")
VohsRPP<-subset(Data, Study == "Vohs" & Version == "RP:P")
PayneRPP<-subset(Data, Study == "Payne" & Version == "RP:P")
RisenRPP<-subset(Data, Study == "Risen" & Version == "RP:P")
ShnabelRPP<-subset(Data, Study == "Shnabel" & Version == "RP:P")
LobueRPP<-subset(Data, Study == "LoBue" & Version == "RP:P")
vanDijkRPP<-subset(Data, Study == "van Dijk" & Version == "RP:P")
ForsterRPP<-subset(Data, Study == "Forster" & Version == "RP:P")

#Median N RP:P
median(sum(Alb5RPP$N),sum(Alb7RPP$N),sum(CrosbyRPP$N),sum(VohsRPP$N),
       sum(PayneRPP$N),sum(RisenRPP$N),sum(ShnabelRPP$N),sum(LobueRPP$N),
       sum(vanDijkRPP$N),sum(ForsterRPP$N))


Alb5REV<-subset(Data, Study == "Albarracin S5" & Version == "Revised")
Alb7REV<-subset(Data, Study == "Albarracin S7" & Version == "Revised")
CrosbyREV<-subset(Data, Study == "Crosby" & Version == "Revised")
VohsREV<-subset(Data, Study == "Vohs" & Version == "Revised")
PayneREV<-subset(Data, Study == "Payne" & Version == "Revised")
RisenREV<-subset(Data, Study == "Risen" & Version == "Revised")
ShnabelREV<-subset(Data, Study == "Shnabel" & Version == "Revised")
LobueREV<-subset(Data, Study == "LoBue" & Version == "Revised")
vanDijkREV<-subset(Data, Study == "van Dijk" & Version == "Revised")
ForsterREV<-subset(Data, Study == "Forster" & Version == "Revised")

#Median N revised
median(sum(Alb5REV$N),sum(Alb7REV$N),sum(CrosbyREV$N),sum(VohsREV$N),
       sum(PayneREV$N),sum(RisenREV$N),sum(ShnabelREV$N),sum(LobueREV$N),
       sum(vanDijkREV$N),sum(ForsterREV$N))

#Median N total

Alb5<-subset(Data, Study == "Albarracin S5")
Alb7<-subset(Data, Study == "Albarracin S7")
Crosby<-subset(Data, Study == "Crosby")
Vohs<-subset(Data, Study == "Vohs")
Payne<-subset(Data, Study == "Payne")
Risen<-subset(Data, Study == "Risen")
Shnabel<-subset(Data, Study == "Shnabel")
Lobue<-subset(Data, Study == "LoBue")
vanDijk<-subset(Data, Study == "van Dijk")
Forster<-subset(Data, Study == "Forster")

median(sum(Alb5$N),sum(Alb7$N),sum(Crosby$N),sum(Vohs$N),
       sum(Payne$N),sum(Risen$N),sum(Shnabel$N),sum(Lobue$N),
       sum(vanDijk$N),sum(Forster$N))

range(sum(Alb5$N),sum(Alb7$N),sum(Crosby$N),sum(Vohs$N),
       sum(Payne$N),sum(Risen$N),sum(Shnabel$N),sum(Lobue$N),
       sum(vanDijk$N),sum(Forster$N))

