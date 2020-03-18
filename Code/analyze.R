
prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"
code.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Code"

setwd(code.dir)
source("helper.R")

setwd(prepped.data.dir)
# interactions of protocol with target effect
di = read.csv("summary_interactions_prepped.csv")

# focal effects at each site
df = read.csv("summary_focal_prepped.csv")
# note that per prep code, "Value" is always on Pearson's r scale, NOT R^2



library(dplyr)
library(ggplot2)
library(MetaUtility)
library(Replicate)
library(robumeta)

###################################### EXPLORE ######################################

# codebook verbatim from Charlie's code:
#Study - Which effect is being replicated? 
#Site - Where were the data collected? 
#Version - Is the effect size from the RP:P or Revised protocol?
#ES.Metric - what is the effect size for this estimate?
#Value - the value of the effect size, with more positive numbers indicating stronger effects in the same direction as the original study
#N - The sample size for that estimate
#Reverse - Does the effect size need to be reversed? (for variance explained effect sizes)
#Label - combination of Study and Version for labeling figures

# number of sites per replication
df %>% group_by(Study) %>%
  summarise( n() )

# and number of sites by revised vs. RPP status
df %>% group_by(Study, Version) %>%
  summarise( n() )


rm("resE")
# make all the subsets to analyze
studies = unique(df$Study)
ql = c(0, MetaUtility::r_to_z(.10) )
boot.reps = 500  # ~~~ INCREASE LATER
digits = 2
z.to.r = TRUE

# 
for( i in studies ) {
  # all replications (RPP and Revised)
  analyze_one_meta( dat = df %>% filter( Study == i ),
                    yi.name = "yi.f",
                    vi.name = "vi.f",
                    meta.name = paste( i, ", all", sep = ""),
                    ql = ql,
                    boot.reps = boot.reps,
                    n.tests = 1,
                    digits = digits,
                    z.to.r = z.to.r )
  # automatically writes to resE
  
  # RPP only
  analyze_one_meta( dat = df %>% filter( Study == i & Version == "RP:P" ),
                    yi.name = "yi.f",
                    vi.name = "vi.f",
                    meta.name = paste( i, ", RP:P", sep = ""),
                    ql = ql,
                    boot.reps = boot.reps,
                    n.tests = 1,
                    digits = digits,
                    z.to.r = z.to.r )
  
  # Revised only
  analyze_one_meta( dat = df %>% filter( Study == i & Version == "Revised" ),
                    yi.name = "yi.f",
                    vi.name = "vi.f",
                    meta.name = paste( i, ", Revised", sep = ""),
                    ql = ql,
                    boot.reps = boot.reps,
                    n.tests = 1,
                    digits = digits,
                    z.to.r = z.to.r )
}

View(resE)


