
###################################### PRELIMINARIES ######################################

prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"
code.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Code"
results.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Results from R"

# helper fns
setwd(code.dir)
source("helper.R")

# for debugging discrepancies
setwd(prepped.data.dir)
# focal effects at each site
df = read.csv("summary_focal_prepped.csv")
# note that per prep code, "Value" is always on Pearson's r scale, NOT R^2

# for checking results
setwd(results.dir)
res = read.csv("results_table_full.csv")

library(dplyr)
library(ggplot2)
library(MetaUtility)
library(Replicate)
library(robumeta)
library(boot)
library(metafor)
library(psych)
library(testthat)

###################################### MANUALLY REPRODUCE ONE META ######################################

# try different subsets and ES.metrics


##### Start from Raw Datasets #####
setwd(raw.data.dir)
rpp = read.csv("RPP Contested Reps.csv")

setwd(raw.data.dir)
Data<-read.csv(file="Summary Effect Sizes - Single Effects.csv",header=TRUE)

# CHOOSE WHICH SUBSET TO ANALYZE HERE:
table(Data$Study)
table(rpp$Authors..O.)
charlie.meta.name = "van Dijk"
rpp.meta.name = "E van Dijk, GA van Kleef, W Steinel, I van Beest"
subset = "RP:P"

# some that I have already checked:
# charlie.meta.name = "Risen"
# rpp.meta.name = "JL Risen, T Gilovich"
# subset = "Revised"

# charlie.meta.name = "Payne"
# rpp.meta.name = "BK Payne, MA Burkley, MB Stokes"
# subset = "all"

# subset the datasets and results
if ( subset == "all" ) Data = Data %>% filter( Study == charlie.meta.name )
if ( subset == "Revised" ) Data = Data %>% filter( Study == charlie.meta.name & Version == "Revised" )
if ( subset == "RP:P" ) Data = Data %>% filter( Study == charlie.meta.name & Version == "RP:P" )

( k = nrow(Data) )

res = res[ res$Meta == paste( charlie.meta.name, ", ", subset, sep = "" ), ]
res

# the prepped dataset for debugging discrepancies
if ( subset == "all" ) df = df %>% filter( Study == charlie.meta.name )
if ( subset == "Revised" ) df = df %>% filter( Study == charlie.meta.name & Version == "Revised" )
if ( subset == "RP:P" ) df = df %>% filter( Study == charlie.meta.name & Version == "RP:P" )


##### Estimates on Fisher's Z Scale #####

# convert to Fisher's z
Data$ES.Metric
if ( Data$ES.Metric[1] %in% c("pseudo R2", "eta sq") ) Data$yif = r_to_z(sqrt(Data$Value))
if ( Data$ES.Metric[1] %in% c("d", "g") ) Data$yif = r_to_z(psych::d2r(Data$Value))
if ( Data$ES.Metric[1] %in% c("r", "rp") ) Data$yif = r_to_z(Data$Value)

# variance of replication estimates
Data$vif = 1 / (Data$N - 3)

# *sanity check
expect_equal( sort(Data$yif), sort(df$yi.f) )

# get estimate and variance of original study
( yio.f = r_to_z( rpp$T_r..O.[ rpp$Authors..O. == rpp.meta.name] ) )
( vio.f = 1 / (rpp$T_N..O.[ rpp$Authors..O. == rpp.meta.name] - 3) )

# *sanity check
expect_equal( yio.f, df$yio.f[1] )
expect_equal( vio.f, df$vio.f[1] )


##### Original and Replication Estimates #####
# meta-analyze replications
meta = robu( yif ~ 1,
             data = Data,
             studynum = 1:k,  # assume no clustering
             var.eff.size = vif,
             small = TRUE)

est = meta$b.r
t2 = meta$mod_info$tau.sq
mu.lo = meta$reg_table$CI.L
mu.hi = meta$reg_table$CI.U
mu.se = meta$reg_table$SE
mu.pval = meta$reg_table$prob


# **check Est, Pval, Tau in results
res$Est; round( z_to_r(est),2 ); round( z_to_r(mu.lo),2 ); round( z_to_r(mu.hi),2 )
res$Pval; mu.pval 
res$Tau; round( sqrt(t2), 2 )


##### Porig #####
Porig = p_orig( orig.y = yio.f,
                  orig.vy = vio.f,
                  yr = est,
                  t2 = t2,
                  vyr = mu.se^2 )

# *sanity check
expect_equal( as.numeric(as.character(res$Porig)), round(Porig, 3) )


##### P(signif agree) #####
Psignif.agree = prob_signif_agree( orig.y = yio.f,  # all entries of this variable will be the same
                                   orig.vy = vio.f,
                                   rep.vy = mu.se^2,
                                   t2 = t2 )

# *sanity check
expect_equal( as.numeric(as.character(res$Psignif.agree)), round(Psignif.agree, 2) )


##### Phat(>0) #####

if( t2 > 0 ){
  
  Phat0 = prop_stronger(q = 0,
                        M = est,
                        t2 = t2,
                        tail = "above",
                        dat = Data,
                        yi.name = "yif",
                        vi.name = "vif")
  
  Phat0.1 = prop_stronger(q = r_to_z(.1),
                        M = est,
                        t2 = t2,
                        tail = "above",
                        dat = Data,
                        yi.name = "yif",
                        vi.name = "vif")
  
  Phat0.2 = prop_stronger(q = r_to_z(.2),
                        M = est,
                        t2 = t2,
                        tail = "above",
                        dat = Data,
                        yi.name = "yif",
                        vi.name = "vif")
  
  
} else {
  Phat0 = est > 0
  Phat0.1 = est > r_to_z(.1)
  Phat0.1 = est > r_to_z(.2)
}


# * sanity check
# point estimates should agree exactly, but CIs may differ due to bootstrapping
res$Percent.above.0; round(100*Phat0$est,1); round(100*Phat0$lo,1); round(100*Phat0$hi,1)
res$Percent.above.0.1; round(100*Phat0.1$est,1); round(100*Phat0.1$lo,1); round(100*Phat0.1$hi,1)
res$Percent.above.0.2; round(100*Phat0.2$est,1); round(100*Phat0.2$lo,1); round(100*Phat0.2$hi,1)



# in case we need to debug the function
# analyze_one_meta( dat = df,  # subset to analyze
#                   yi.name = "yi.f",  # variable name of point estimate in replications (on Fisher's z scale)
#                   vi.name = "vi.f", # variable name of variance estimate in replications (on Fisher's z scale)
#                   meta.name = charlie.meta.name,  # name of meta-analysis/project
#                   
#                   ql = ql,  # on Fisher's z scale
#                   z.to.r = z.to.r,  # should we transform back to r for reporting?
#                   boot.reps = 2000,
#                   digits = 2 )




