
################################ FOR ANALYSIS ################################

# analyze a subset 
analyze_one_meta = function( dat,  # subset to analyze
                             yi.name,  # variable name of point estimate in replications (on Fisher's z scale)
                             vi.name, # variable name of variance estimate in replications (on Fisher's z scale)
                             meta.name,  # name of meta-analysis/project
                             
                             ql,  # on Fisher's z scale
                             z.to.r,  # should we transform back to r for reporting?
                             boot.reps = 2000,
                             digits = 2 ) {  # digits for rounding
  
  # #~~~ TEST ONLY:
  # dat = df %>% filter( Study == "Albarracin S5" & Version == "RP:P" )  # k = 1 case
  # dat = df %>% filter( Study == "Albarracin S7" )  # k = 14
  # 
  # dat = df %>% filter( Study == "Albarracin S7" & Version == "Revised" ) 
  # yi.name = "yi.f"
  # vi.name = "vi.f"
  # meta.name = "Fake"
  # ql = c(0, MetaUtility::r_to_z(.10) )
  # boot.reps = 500
  # digits = 2
  # z.to.r = TRUE
  # # ~~~

  dat$yi = dat[[yi.name]]
  dat$vyi = dat[[vi.name]]
  
  k = nrow(dat)
  
  ##### Meta-Analyze the Replications #####
  # prepare to catch errors with robumeta
  robu.error = NA

  # fine to proceed with meta-analysis even if k=1; it just returns the
  # study's yi.f and vi.f and sets t2 = 0
  tryCatch({
    # sometimes hits errors
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
    
  }, error = function(err){
    robu.error <<- err$message
    
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
  })
  

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
  

  ##### Porig #####
  Porig = p_orig( orig.y = dat$yio.f[1],  # they will all be the same
          orig.vy = dat$vio.f[1],
          yr = est,
          t2 = t2,
          vyr = mu.se^2 )
  
  ##### Put Results in Dataframe #####
  # transform back to r if needed
  if (z.to.r == TRUE) {
    est = z_to_r(est)
    lo = z_to_r(mu.lo)
    hi = z_to_r(mu.hi)
  }
  
  est.string = paste( round( est, digits ),
                      format_CI( lo, 
                                 hi,
                                 digits),
                      sep = " " )
  
  tau.string = round( sqrt(t2), digits)
  
  
  new.row = data.frame( Meta = meta.name,
                        k = k,
                        Est = est.string,
                        Pval = format_stat(mu.pval, cutoffs = c(.1, .0001) ),
                        Tau = tau.string,
                        Porig = round( Porig, digits = digits ) )

  # add Phat results to dataframe
  # tail is now just for the purpose of creating the column name
  tail = rep("above", length(unlist(ql)))
  tail[ unlist(ql) < 0 ] = "below"
  # transform q back to r, if desired, for the column name
  if (z.to.r == TRUE) q.vec = z_to_r(unlist(ql)) else q.vec = unlist(ql)
  Phat.names = paste( "Percent ", tail, " ", q.vec, sep = "" )
  new.row[ , Phat.names ] = Phat.df$string

  # add unrounded results to facilitate reporting overall stats
  Phat.names.2 = paste( Phat.names,
                        "unrounded",
                        sep = " ")
  new.row[ , Phat.names.2 ] = Phat.df$Est

  # details at the end of df, to easily lop off for prettiness
  new.row$Porig.unrounded = Porig
  new.row$robu.error = robu.error
  
  # resE is a global variable
  if ( !exists("resE") ){
    resE <<- new.row
  } else {
    library(plyr)
    resE <<- rbind.fill(resE, new.row)
    detach("package:plyr", unload=TRUE)
  }
} 




################################ MISCELLANEOUS ################################

# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
update_result_csv = function( name,
                              section = NA,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    res = read.csv( "stats_for_paper.csv",
                    stringsAsFactors = FALSE,
                    colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% res$name) ) res[ res$name %in% name, ] = new.rows
    else res = rbind(res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    res = new.rows
  }
  
  write.csv( res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  setwd(overleaf.dir)
  write.csv( res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  if ( print == TRUE ) {
    View(res)
  }
}