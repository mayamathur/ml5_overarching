
prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"
# the rest is verbatim from Charlie's code, "ML5 Overarching Analyses.R"

############################# VERBATIM FROM CHARLIE'S SCRIPT ############################# 

#Reading data
#setwd("C:/Users/Charlie/Dropbox/Charlie_Harddrive/Desktop/Ongoing Projects/Many Labs 5/Overarching Paper/Results Included")
setwd(raw.data.dir)
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

# calculate Fisher's z
Data = escalc(ri=Value,ni=N,data=Data,measure="COR",vtype="LS")
names(Data)[ names(Data) == "yi" ] = "yi.f"
names(Data)[ names(Data) == "vi" ] = "vi.f"

##### Write Prepped Data #####
setwd(prepped.data.dir)
write.csv(Data, "summary_focal_prepped.csv", row.names = FALSE)



############################# PREP INTERACTION EFFECT DATA #############################

#Reading in the data
setwd(raw.data.dir)
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


##### Write Prepped Data #####
setwd(prepped.data.dir)
write.csv(Data2, "summary_interaction_prepped.csv", row.names = FALSE)


