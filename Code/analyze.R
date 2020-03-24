

###################################### PRELIMINARIES ######################################

prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"
code.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Code"
results.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Results from R"

setwd(code.dir)
source("helper.R")

setwd(prepped.data.dir)
# # interactions of protocol with target effect
# di = read.csv("summary_interactions_prepped.csv")

# focal effects at each site
df = read.csv("summary_focal_prepped.csv")
# note that per prep code, "Value" is always on Pearson's r scale, NOT R^2

library(dplyr)
library(ggplot2)
library(MetaUtility)
library(Replicate)
library(robumeta)
library(boot)

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
ql = c(0,
       MetaUtility::r_to_z(.10),
       MetaUtility::r_to_z(.20) )
boot.reps = 2000 # ~~~ increase
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

# # debugging
# analyze_one_meta( dat = df %>% filter( Study == "Albarracin S5" & Version == "RP:P" ),
#                   yi.name = "yi.f",
#                   vi.name = "vi.f",
#                   meta.name = paste( i, ", Revised", sep = ""),
#                   ql = ql,
#                   boot.reps = boot.reps,
#                   n.tests = 1,
#                   digits = digits,
#                   z.to.r = z.to.r )

# add variable for which group the estimate is in
resE$subset = unlist( lapply( X = resE$Meta,
                              FUN = function(x) strsplit( as.character(x), ", " )[[1]][[2]] ) )


# save results
setwd(results.dir)
write.csv(x = resE,
          file = "results_table_full.csv",
          row.names = FALSE)

# pretty version without extra columns
temp = resE[ , which( grepl( pattern = "unrounded", names(resE) ) == FALSE ) ]
temp = temp %>% select(-c(robu.error))
write.csv(x = temp,
          file = "results_table_pretty.csv",
          row.names = FALSE)


# subset to analyses of all replications combined
write.csv(x = temp %>% filter(subset == "all"),
          file = "results_table_pretty_combined.csv",
          row.names = FALSE)


# and to analyses of RP:P replications only
write.csv(x = temp %>% filter(subset == "RP:P"),
          file = "results_table_pretty_rpp.csv",
          row.names = FALSE)

# and to analyses of Revised only
write.csv(x = temp %>% filter(subset == "Revised"),
          file = "results_table_pretty_revised.csv",
          row.names = FALSE)


###################################### SUMMARY STATS ######################################

# summary table by subset
t = resE %>% group_by(subset) %>%
  summarise( k.reps = sum(k),
             n.studies = n(),
             Tau.mn = mean(Tau),
             Perc.Tau.GT0 = round( 100 * mean(Tau>0), 0 ),
             Phat0.mn = round( mean(`Percent above 0 unrounded`), 0 ),
             Phat10.mn = round( mean(`Percent above 0.1 unrounded`), 0 ),
             Phat20.mn = round( mean(`Percent above 0.2 unrounded`), 0 ),
             
             Porig.mn = round( mean(Porig.unrounded), digits ) )
View(t)
setwd(results.dir)

write.csv(x = t,
          file = "results_aggregated_by_subset.csv",
          row.names = FALSE)


