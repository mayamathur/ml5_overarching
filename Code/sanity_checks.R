
# check how datasets are merged

###################################### PRELIMINARIES ######################################

prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Prepped data"
raw.data.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Data/Raw data from Charlie"
code.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Code"
results.dir = "~/Dropbox/Personal computer/Independent studies/Many Labs 5 (ML5)/Charlie's overarching manuscript/MM analyses for ML5 overarching (git)/Results from R"

# helper fns
setwd(code.dir)
source("helper.R")

setwd(prepped.data.dir)

# focal effects at each site
df = read.csv("summary_focal_prepped.csv")
# note that per prep code, "Value" is always on Pearson's r scale, NOT R^2

library(dplyr)
library(ggplot2)
library(MetaUtility)
library(Replicate)
library(robumeta)
library(boot)

###################################### MANUALLY REPRODUCE ONE META ######################################

dat = df %>% filter( Study == "Crosby")  
yi.name = "yi.f"
vi.name = "vi.f"
ql = c(0, MetaUtility::r_to_z(.10) )
boot.reps = 500
digits = 2
z.to.r = TRUE

# standardize variable names
dat$yi = dat[[yi.name]]
dat$vyi = dat[[vi.name]]

k = nrow(dat)

##### Meta-Analyze the Replications #####
library(robumeta)
meta = robu( yi ~ 1,
             data = dat,
             studynum = 1:nrow(dat),  # assume no clustering
             var.eff.size = vyi,
             small = TRUE)

est = meta$b.r
t2 = meta$mod_info$tau.sq
mu.lo = meta$reg_table$CI.L
mu.hi = meta$reg_table$CI.U
mu.se = meta$reg_table$SE
mu.pval = meta$reg_table$prob
  

  # use rma.uni instead (parametric)
  meta <<- rma.uni( yi = yi,
                    vi = vyi,
                    data = dat,
                    test = "knha" )
  est <<- meta$b
  t2 <<- meta$tau2
  mu.lo <<- meta$ci.lb
  mu.hi <<- meta$ci.ub
  mu.se <<- meta$se
  mu.pval <<- meta$pval


##### Porig #####
Porig = p_orig( orig.y = dat$yio.f[1],  # they will all be the same
                orig.vy = dat$vio.f[1],
                yr = est,
                t2 = t2,
                vyr = mu.se^2 )

##### Probability of Significance Agreement #####
# P(significance agreement) for pooled replication estimate
Psignif.agree = prob_signif_agree( orig.y = dat$yio.f[1],  # all entries of this variable will be the same
                                   orig.vy = dat$vio.f[1],
                                   rep.vy = mu.se^2,
                                   t2 = t2 )

##### NPPhat #####
# skip this if k=1 or if there is no heterogeneity
if ( t2 > 0 ) {
  Phat.l = lapply( ql,
                   FUN = function(q) {
                     
                     # get new ensemble estimates for this subset
                     # yi and vyi aren't using yi.name and vi.name intentionally 
                     #  since these are newly created variables
                     ens = MetaUtility::calib_ests( yi = dat$yi, 
                                                    sei = sqrt(dat$vyi) )
                     
                     # set tail based on sign of q
                     if (q >= 0) tail = "above" else tail = "below"
                     if ( tail == "above" ) Phat.NP.ens = sum(ens > c(q)) / length(ens)
                     if ( tail == "below" ) Phat.NP.ens = sum(ens < c(q)) / length(ens)
                     
                     Note = NA
                     boot.lo.ens = NA  # new
                     boot.hi.ens = NA
                     tryCatch({
                       boot.res.ens = boot( data = dat, 
                                            parallel = "multicore",
                                            R = boot.reps, 
                                            statistic = function(original, indices) {
                                              
                                              b = original[indices,]
                                              
                                              ens.b = MetaUtility::calib_ests( yi = b$yi, 
                                                                               sei = sqrt(b$vyi) )
                                              if ( tail == "above" ) return( sum(ens.b > c(q)) / length(ens.b) )
                                              if ( tail == "below" ) return( sum(ens.b < c(q)) / length(ens.b) )
                                            }
                       )
                       
                       bootCIs.ens = boot.ci(boot.res.ens, type="bca")
                       boot.lo.ens = bootCIs.ens$bca[4]
                       boot.hi.ens = bootCIs.ens$bca[5]
                       
                     }, error = function(err){
                       boot.lo.ens <<- NA
                       boot.hi.ens <<- NA
                       print( paste(meta.name, ": ", err$message, sep = " ") )
                       Note <<- err$message
                       
                     }, warning = function(w) {
                       # catch "extreme order statistics used as endpoints"
                       boot.lo.ens <<- NA
                       boot.hi.ens <<- NA
                       print( paste(meta.name, ": ", w$message, sep = " ") )
                       Note <<- w$message
                     }
                     
                     )  # end tryCatch
                     
                     return( data.frame( Est = Phat.NP.ens,
                                         lo = boot.lo.ens,
                                         hi = boot.hi.ens,
                                         boot.note = Note ) )
                   } )  # end lapply
  
  
  Phat.df = do.call( rbind, 
                     Phat.l )
  Phat.df$string = paste( round( 100*Phat.df$Est,
                                 digits = 0 ),
                          format_CI( 100*Phat.df$lo, 
                                     100*Phat.df$hi,
                                     digits = 0 ),
                          sep = " " )
  # omit CI if it was NA
  Phat.df$string[ is.na(Phat.df$lo) ] = round( 100*Phat.df$Est[ is.na(Phat.df$lo) ],
                                               digits = 0 )
  # if t2 = 0 exactly: 
} else {
  
  # just check if Phat = 0 or 1 but omit CI
  Phat.df = data.frame( Est = 100 * (c(est) > ql),
                        lo = rep( NA, length(ql) ),
                        hi = rep( NA, length(ql) ),
                        boot.note = rep( "t2 = 0 exactly, so no Phat CI", length(ql) ) )
  
  Phat.df$string = Phat.df$Est # omit the CI
  # Phat.df$string = paste( Phat.df$Est, 
  #                         " [NA, NA]", 
  #                         sep = "")
}
  
  
  
Porig
Psignif.agree
Phat.df
