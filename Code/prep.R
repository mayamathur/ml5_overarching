
prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"

library(metafor)
library(psych)
library(MetaUtility)
library(dplyr)

############################# PREP CHARLIE'S DATA ############################# 

# verbatim from his "ML5 Overarching Analyses.R" script
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


############################# PREP RPP DATA ############################# 

setwd(raw.data.dir)
rpp = read.csv("RPP Contested Reps.csv")

# correlation and sample size of original study
# look at RPP variables
rpp$T_r..O.
rpp$T_N..O.; rpp$N..O. # should be the same
rpp$Authors..O.

rpp$N..O. = as.numeric( as.character(rpp$N..O.) )


# recode author variable to make merger variable
# some contested replications in RPP weren't redone in ML5, so no need to recode those
rpp$Study = rep(NA, nrow(rpp))
rpp$Study[ rpp$Authors..O. == "D Albarrac\x90n, IM Handley, K Noguchi, KC McCulloch, H Li, J Leeper, RD Brown, A Earl, WP Hart" &
             rpp$Replicated.study.number..R. == 7 ] = "Albarracin S7"

rpp$Study[ rpp$Authors..O. == "D Albarrac\x90n, IM Handley, K Noguchi, KC McCulloch, H Li, J Leeper, RD Brown, A Earl, WP Hart" &
             rpp$Replicated.study.number..R. == 5 ] = "Albarracin S5"

rpp$Study[ rpp$Authors..O. == "JR Crosby, B Monin, D Richardson" ] = "Crosby"

rpp$Study[ rpp$Authors..O. == "BK Payne, MA Burkley, MB Stokes" ] = "Payne"

rpp$Study[ rpp$Authors..O. == "E van Dijk, GA van Kleef, W Steinel, I van Beest" ] = "van Dijk"

rpp$Study[ rpp$Authors..O. == "J F_rster, N Liberman, S Kuschel" ] = "Forster"

rpp$Study[ rpp$Authors..O. == "JL Risen, T Gilovich" ] = "Risen"

rpp$Study[ rpp$Authors..O. == "N Shnabel, A Nadler" ] = "Shnabel"

rpp$Study[ rpp$Authors..O. == "KD Vohs, JW Schooler" ] = "Vohs"

rpp$Study[ rpp$Authors..O. == "V Lobue, JS DeLoache" ] = "LoBue"


# should be none
unique( Data$Study[ !Data$Study %in% rpp$Study ] )

# remove unneccessary rows
rpp = rpp[ rpp$Study %in% Data$Study, ]
nrow(rpp)  # should be 10

# make my own point estimate and variance variables
rpp$yio.f = r_to_z( rpp$T_r..O. )
rpp$vio.f = 1 / (rpp$N..O. - 3)

# remove unwanted columns
rpp = rpp %>% select( Study, 
                      T_r..O.,
                      N..O.,
                      yio.f,
                      vio.f )

# merge them
d = merge( Data,
           rpp,
           by = "Study",
           all.x = TRUE )

nrow(d)  # should be 101, just like Data


##### Write Prepped Data #####
setwd(prepped.data.dir)
write.csv(d, "summary_focal_prepped.csv", row.names = FALSE)
