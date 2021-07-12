library(doParallel)
library(survival)

gen_subs = function(cfg, parameters, nsub){


 # Creating subject IDs and the weight covariate
 WT = rnorm(n = nsub,
                mean = 79.7,
                sd   = 20.3)

 

 SUBS = data.frame(SIMINT_ID   = 1:nsub,
                   SIMINT_TIME = 0,
                   WT          = WT) %>% 
  mutate(#Cap lower bound of WT at 30
         #To maintain positive values for WT
         WT = case_when(WT < 30 ~ 30,
                        TRUE ~ WT))

  # Now we call simulate_subjects to setting ponly to TRUE to generate the
  # subject level parameters based on the specified IIV

  cfg = system_set_option(cfg, group="stochastic", option="nsub",    value=nsub)
# cfg = system_set_option(cfg, group="stochastic", option="ci",      value=95 )
# cfg = system_set_option(cfg, group="stochastic", option="seed",    value=8675309)
  cfg = system_set_option(cfg, group="stochastic", option="ponly",   value=TRUE)
# cfg = system_set_option(cfg, group="stochastic", option="states",  value=list())
# cfg = system_set_option(cfg, group="stochastic", option="outputs", value=c("OP1", "OP2"))

  som   = simulate_subjects(parameters, cfg)


  # Putting a cap on the body weight used for dosing
  BWdose               = WT
  BWdose[BWdose > 100] = 100
  # Now we append the subject parameters to the SUBS data.frame
  SUBS = cbind(SUBS, as.data.frame(som$subjects$parameters)) %>%
      mutate(BWdose = BWdose)  # This replaces the default weight for dose
                               # calculations with the actual subject weight
                               # or 100 which ever is less.


SUBS}


#' @param nsub   number of subject to sample for each trial
#' @param ntrial number of trials to simulate
#' @param dsname name of the dataset to sample from
#' @param multicore Boolean variable indicating if you want to use multiple cores (TRUE) or a single core (FALSE)
#' @param scenario  character description of the dosing regimen, parametric perturbations, etc
#' @param tfinal final time to simulate PN (days)
#' @param tinclude output times (days) to include in the simulated output the default value of \code{NULL} will not include any user specified times
#' @param dtimes dosing times (days)
#' @param dvals  corresponding dose values (mg/kg)
#' @param inf_dur infusion duration (hours)
#' @param cycle_dur  duration of treatment cycle (days)
#' @param ncycles  number of cycles (#)
sim_pn      = function(cfg, 
                        nsub      = 100,
                        ntrial    = 2, 
                        dsname    = "SUBS",
                        multicore = FALSE,
                        scenario  = "Regimen Description",
                        tfinal    = 3*10*7,
                        tinclude  = NULL, 
                        dtimes    = 0,
                        dvals     = 2.0,
                        inf_dur   = 0.5,
                        cycle_dur = 28,
                        ncycles   = 1.0){

  times_hr_keep    = seq(24,tfinal*24,24)
  times_hr_keep    = c(0.001, times_hr_keep)
  out_keep         = c("DVADC", "SUR")


  # If the user specified times to include we add them here:
  if(!is.null(tinclude)){
    # First we remove any times beyond tfinal here:
    tinclude = tinclude[tinclude < tfinal]
    if(length(tinclude) > 0){
      tinclude_hr = tinclude*24
      times_hr_keep = sort(unique(c(times_hr_keep, tinclude_hr)))
    }
  }

  # Checking some user input:
  isgood = TRUE 

  if( cycle_dur*(ncycles) > tfinal){
    vp(cfg, paste("The number of days in the dosing cycles you've specified (",cycle_dur*(ncycles) , ")", sep=""))
    vp(cfg, paste("is greater than the duration of the simulation (",tfinal , " days)", sep=""))
  }
  
  if( length(dtimes) != length(dvals)){
    isgood = FALSE
    vp(cfg, paste("The number of dosing times (",length(dtimes), ") must be equal to the number of dosing values (", length(dvals), ")", sep=""))
  
  }

  
  
  
  
  if(isgood){


    # Making sure we have infusion durations for each dosing time
    if(length(inf_dur) == 1){
       inf_dur = rep(inf_dur, times=length(dtimes))
    }

    # expanding the single cycle dosing out to multiple cycles
    if(ncycles > 1){
       dtimes_tot  = c()
       dvals_tot   = c()
       inf_dur_tot = c()

       for(cycle in 1:ncycles){
         toffset = cycle_dur*(cycle - 1)
         dtimes_tot  = c(dtimes_tot,   dtimes+toffset)
         dvals_tot   = c(dvals_tot,    dvals)
         inf_dur_tot = c(inf_dur_tot,  inf_dur)
       }
    
    } else {

       dtimes_tot  = dtimes
       dvals_tot   = dvals
       inf_dur_tot = inf_dur
    }


    # Saving final dosing times and amounts to return below
    dosing_summary = tibble(times    = dtimes_tot,
                            amounts  = dvals_tot,
                            DREG     = scenario)
    
    
    #-----------------------------------------------------
    # General simulation optoins
    # Output times to include
    cfg = system_set_option(cfg, group  = "simulation", 
                                 option = "output_times", 
                                 c(0,times_hr_keep))

    # Applying the infusion rate dosing:
    cfg = pn_set_doses(cfg, dtimes   = dtimes_tot,
                            dvals    = dvals_tot,
                            inf_dur  = inf_dur_tot)
    #-----------------------------------------------------
    # Population simulation options
    # Defining the dataset to sample from
    cfg=system_set_option(cfg, group  = "stochastic",
                               option = "sub_file",
                               value  = dsname)
    
    # Sampling with replacement
    cfg=system_set_option(cfg, group  = "stochastic",
                               option = "sub_file_sample",
                               value  = "with replacement")
    
    # Defining the states and outputs to keep from population simulations
    cfg = system_set_option(cfg, group="stochastic", option="ssp",     value=list())
    cfg = system_set_option(cfg, group="stochastic", option="states",  value=list())
    cfg = system_set_option(cfg, group="stochastic", option="outputs", value=out_keep)
    
    # Number of subjects to sample 
    cfg = system_set_option(cfg, group="stochastic", option="nsub",    value=nsub)
    
    #-----------------------------------------------------
    # Parallel options
    if(multicore){
      cfg=system_set_option(cfg, group  = "simulation",
                                 option = "parallel",    
                                 value  = "multicore")
      
      cfg=system_set_option(cfg, group  = "simulation",
                                 option = "compute_cores", 
                                 value  = detectCores() - 1)
    }
    
    #-----------------------------------------------------
    # Parameters will be pulled from the defined dataset, but the function
    # requires a parametric input, so we just pull the default parameter set
    # and send it it. 
    parameters = system_fetch_parameters(cfg)
    
    # Running population simulations:
    times_day_keep = times_hr_keep/24
    
    simall = matrix(nrow = length(times_hr_keep)*ntrial*nsub,
                    ncol = 4 + length(out_keep))
    
    colnames(simall) = c("ID", "Trial", "time", "time_days", out_keep)
    


    set.seed(8675309)

    # Creating a separate seed for each trial
    myseeds = runif(ntrial)

    start_time = as.numeric(as.POSIXct(Sys.time()))
    for(tidx in 1:ntrial){
       # Setting a trial specific seed:
       cfg = system_set_option(cfg, group="stochastic", option="seed",    value=myseeds[tidx])
    
       # For details on the simulation output format see the population section of
       # the Simulation vignetted:
       # vignette("Simulation")
       # ODE trial simulation
       vp(cfg, paste("Trial ", tidx, " of ", ntrial, " (", date(), ")",  sep=""))
       som   = simulate_subjects(parameters, cfg)
    
    
       odata = list()
       # Resampling at the observed output times:
       for(oname  in out_keep){
         # pulling out the subjects for the current output
         # By default there is a column for each time and a row for each subject,
         # we're going to transpose t() that here:
         odata[[oname]] = t(som[["outputs"]][[oname]])
        
         # Keeping the output times that are important
         odata[[oname]] =  odata[[oname]][ som$times$ts.hours %in%  times_hr_keep, ]
        
       }
    
       for(sidx in 1:nsub){
         row_offset = (tidx-1)*length(times_hr_keep)*nsub + length(times_hr_keep)*(sidx -1) + 1
         cr = row_offset:(row_offset+length(times_hr_keep)-1)
    
    
         # Calculating a unique subject ID based on the trial and the nominal id
         tsid = tidx*10000 + sidx
         
         simall[cr,1] = tsid 
         simall[cr,2] = tidx
         simall[cr,3] = times_hr_keep
         simall[cr,4] = times_day_keep
    
         for(oidx in 1:length(out_keep)){
           oname = out_keep[oidx]
           simall[cr,(oidx+4)] = odata[[oname]][,sidx]
         }
       }
    }
    stop_time = as.numeric(as.POSIXct(Sys.time()))
    #-----------------------------------------------------
    # Processing the survival data
    
    # This contains the timecourse
    simall = as_tibble(simall) %>% 
        mutate(DREG = scenario) # Appending the scenario information 

    
    
    tte = list()
    
    # We do this for each trial
    for(tidx in unique(simall[["Trial"]])){
    
      outsub = simall %>% filter(Trial == tidx)
      
      times = unique(outsub$time)
    
      # Get individual survival curves
      outsplt <- split(outsub$SUR, outsub$ID)
      
      
      preds <- data.frame(do.call(cbind,outsplt))
      
      # Simulate from uniform distribution for survival values
      Svals = runif(unique(outsub$ID))
      
      # Calculate corresponding survival times from the model
      sims <- sapply(1:length(Svals),function(i) {
        tryCatch(approx(x=c(preds[,i]),y=c(times),xout=Svals[i])$y, error=function(e) NA)
      })
    
      
      # Censor at maximum observation time
      status <- ifelse(is.na(sims),0,1)
      sims[is.na(sims)] <- max(times)
      
      simout <- data.frame(ID=sort(unique(outsub$ID)),
                           sims=sims,
                           status=status,
                           irep=unique(outsub$Trial)) %>%
       dplyr::left_join(outsub %>% dplyr::distinct(ID,DREG))
    
      # Return data frame of status and observation time
      # Stacking them all together
      tte[[tidx]]  = simout
    }

    
    pn_tc =  make_plotssum(tte, 
               output_times = c(0,times_hr_keep),
               endTIME      = tfinal*24)

    # Adding scenario columns to the outputs
    pn_tc[["DREG"]] = scenario

    # Creating the median and  prediction interval timecourse dataset
    plo = 0.025
    phi = 0.975
    pi_SUR_tc = simall %>% 
            group_by(time)  %>%
            dplyr::summarize(pred_qlo=quantile(SUR,prob=plo),
                             pred_med=median(SUR),
                             pred_qhi=quantile(SUR,prob=phi)) %>%
            mutate(name="Survival") %>%
            mutate(time_days=time/24) %>%
            mutate(DREG=scenario) %>%
            ungroup()
    pi_ADC_tc = simall %>% 
            group_by(time)  %>%
            dplyr::summarize(pred_qlo=quantile(DVADC,prob=plo),
                             pred_med=median(DVADC),
                             pred_qhi=quantile(DVADC,prob=phi)) %>%
            mutate(name="ADC") %>%
            mutate(time_days=time/24) %>%
            mutate(DREG=scenario) %>%
            ungroup()
     
    pi_tc = rbind(pi_SUR_tc, pi_ADC_tc)

    res = list(pi_tc          = pi_tc,
               tte            = tte,
               pn_tc          = pn_tc,
               dosing_summary = dosing_summary,
               system_time    = stop_time - start_time,
               isgood         = isgood)

  } else {
    res = list(isgood = isgood)
  }

}



#' @export
#' @title 
#' @description 
#'
#' @param cfg ubiquity system object  
#' @param dtimes dosing times (days)
#' @param dvals  corresponding dose values (mg/kg)
#' @param inf_dur infusion duration (hours)
#'
#'@return vector of numbers from \code{a} to \code{b} with
#'\code{n} linearly spaced apart
pn_set_doses =  function(cfg, 
                         dtimes    = 0,
                         dvals     = 2.0,
                         inf_dur   = 0.5){

  # Converting dose times in days to hours
  dtimes_hr   = dtimes*24


  if(length(inf_dur) == 1){
    inf_dur = rep(inf_dur, length(dtimes))
  }

  rate_times  = c()
  rate_levels = c()

  # Now we conver the absolute doses in mg/kg to rates (mg/kg/hr) over the
  # infusion duration 
  for(ridx in 1:length(dtimes)){
    rate_times  = c(rate_times,  c(dtimes_hr[ridx],           (dtimes_hr[ridx] + inf_dur[ridx])))
    rate_levels = c(rate_levels, c(dvals[ridx]/inf_dur[ridx],   0                              ))
  }

  # Clear out any current dosing and apply the rates above
  cfg = system_zero_inputs(cfg)
  cfg = system_set_rate(cfg, rate   = "Dinf", 
                             times  = rate_times,    # hours  
                             levels = rate_levels)   # mg/hr 
  
cfg}


find_resources = function(){

  res = NULL 
  if(dir.exists("resources")){
    res = "resources" }

  if(dir.exists(file.path("..", "resources"))){
    res = file.path("..", "resources") }

res}

QC_load = function(cfg){
 
  mg_to_nmoles = (1/145455) * 3.585 * (1000*1000)

  file_data   = file.path(find_resources(), "metrum","01PN-app", "lab-notebook", "testres2.RData")
  file_app    = file.path(find_resources(), "metrum","01PN-app", "apps","TTE", "app.R")
  file_model  = file.path(find_resources(), "metrum","01PN-app", "modeling","polaPN", "1.cpp")
  file_rmd    = file.path(find_resources(), "metrum","01PN-app", "lab-notebook", "exploration2.Rmd")
  file_input  = file.path(find_resources(), "metrum","01PN-app", "lab-notebook", "input_data.RData")
  file_pn     = file.path(find_resources(), "metrum","01PN-app", "scripts", "simobj", "plotsumVarCL.RDS")

  isgood = TRUE  


  data_pn_tc = NULL
  # Reading in the Metrum TTE PN timecourse predictions
  if(file.exists(file_pn)){
    data_pn_tc  = readRDS(file_pn) %>%
       data.frame() %>%
       mutate(group = stringr::str_replace(.$group,'DREG=','')) %>%
       mutate(med  = 1 - med) %>%
       mutate(qlo  = 1 - qlo) %>%
       mutate(qhi  = 1 - qhi) %>%
       rename(DREG=group)
    
       data_pn_tc[data_pn_tc$DREG == "2.4 mg/kg Q3W (21 day cycle): CL*1", ]$DREG = "2.4 mg/kg Q3W (21 day cycle): CL*1.0"
  } else {
    isgood = FALSE
    vp(cfg, "Unable to find peripherial neuropathy timecourse data:")
    vp(cfg, paste(" -->", file_pn))
    vp(cfg, "pn_tc field will be NULL")
  }

  res_all     = NULL
  dose_all    = NULL
  subs_all    = NULL

  if(file.exists(file_data)){
    load(file_data)
    for(trial in 1:length(testres2)){
      
      subs      = as_tibble(testres2[[trial]][[4]]$subs)
      subs      = subs %>% mutate(trial=trial)
      subs_all  = rbind(subs_all,  subs)
      
      dose    =           testres2[[trial]][[4]]$SimData
      if(!is.null(dose)){
        dose    = dose %>% mutate(trial=trial)
        dose_all  = rbind(dose_all, dose)
      }
    
      res     =           testres2[[trial]][[4]]$res
      if(!is.null(res)){
        res     = res %>% mutate(trial=trial)
        res_all   = rbind(res_all,  res)
      }
    }
    
    res_all     = res_all   %>% mutate(unique_ID = paste("Trial: ", trial, ", ID: ", ID, sep=""))
    dose_all    = dose_all  %>% mutate(unique_ID = paste("Trial: ", trial, ", ID: ", ID, sep=""))
    subs_all    = subs_all  %>% mutate(unique_ID = paste("Trial: ", trial, ", ID: ", ID, sep=""))
    
    parameters = system_fetch_parameters(cfg)
    CL_BL = parameters$CL
    
    # Adding columns to the input data to make sure I understand the relationships correctly:
    # Conversion defined in: packages/internal/R/data-gen.R
      dose_all    = dose_all                       %>%
            mutate(CL_TV   = THETA1)               %>% # Mappign TV of CL to THETA1
            mutate(CL_mult = CL_TV/CL_BL)          %>% # Calculating the CL multiplier
            mutate(DOSE_mpk = amt/WT/mg_to_nmoles) %>% # Converting the AMT in nmoles to mg/kg
            mutate(INF_DUR_hr = amt/rate)              # Duration of infusion in hours
    
       subs_all = subs_all %>% 
            mutate(CL_mult = CL_TV/CL_BL)              # Calculating the CL multiplier
       
  } else {
    isgood = FALSE
    vp(cfg, "Unable to find timecourse simulations data:")
    vp(cfg, paste(" -->", file_data))
    vp(cfg, "sims, dosing, and subjects fields will be NULL")
  }


  if(!isgood){
    vp(cfg, "QC_load()")
    vp(cfg, "Unable to find some or all data, see above.")
  }
  
  res = list(mg_to_nmoles = mg_to_nmoles,
             data     = list( sims     = res_all ,
                              dosing   = dose_all,
                              pn_tc    = data_pn_tc  , 
                              subjects = subs_all)   ,
             files    = list( data     = file_data   ,
                              app      = file_app    ,
                              model    = file_model  ,
                              rmd      = file_rmd    ,
                              pn_tc    = file_pn     ,
                              input    = file_input )) 
res}




#---------------------------------------------------------------------------
# The following functions were taken from the Metrum code 
#---------------------------------------------------------------------------


#' make plots summary
#' @param simout list of individual simulations
#' @param output_times prediction times in hours
#' @param endTIME - max time to simulate survival probability out to
#' @export
make_plotssum <- function(
  simout,
  output_times,
  endTIME = 31*7*24
) {
  
  
  #times_for_prediction = sort(c(seq(0,endTIME,length=200)))
  times_for_prediction = sort(unique(output_times))
  
  model <- 'Surv(TIME,DV) ~ DREG'

  simpred <- purrr::map(model, ~ osPredictiveDistribution(simData = simout,
                                                          SurvFit=.x,
                                                          pred.times=times_for_prediction))
  simpred[[1]]$group = as.character(simpred[[1]]$group)
  
  # Calculate quantiles from simulation
  plo=0.025
  phi = 0.975
  
  quants3 <- simpred[[1]] %>% dplyr::group_by(group,time) %>%
    dplyr::summarize(qlo=1-quantile(surv,prob=plo),
                     med=1-median(surv),
                     qhi=1-quantile(surv,prob=phi))
  
  return(quants3)
  
}
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
# @param simData data frame of simulated data from simOSdata function
# @param origSurvFit  K-M formula to follow
# @param pred.times numeric vector of times at which to obtain predicted survival
# @return dataframe with kaplain meir from 0 to max time survival info by replicate
osPredictiveDistribution <- function(simData,
                                     SurvFit,
                                     pred.times) {
  
  #Model for observed K-M
  mod.km <- formula(SurvFit)
  
  #Simulation model (this expression should have the same r.h.s as the KM model just prior)
  mod.sim <- update(mod.km, Surv(sims,status) ~ .)
  
  # Define function for getting predictions at a grid of time points
  predfun <- function(simobj,predtimes) {
    # Get K-M estimate from simulated data
    simfit <- survival::survfit(mod.sim,data=simobj)
    
    # Obtain predictied survival at the grid of times
    est <- numeric()
    ifelse(is.null(simfit$strata),
           strata.tot <- length(simfit$surv),
           strata.tot <- cumsum(simfit$strata))
    
    istart <- c(1,strata.tot+1)
    iend <- strata.tot
    nstrat <- ifelse(is.null(simfit$strata),1,length(simfit$strata))
    
    for (i in 1:nstrat) {
      est <- c(est,approx(x=c(0,simfit$time[istart[i]:iend[i]]),
                          y=c(1,simfit$surv[istart[i]:iend[i]]),
                          xout=pred.times,rule=2)$y)
    }
    
    ifelse(is.null(simfit$strata),
           bootsim3 <- data.frame(time=predtimes,
                                  group=1,
                                  surv=est,
                                  irep = simobj$irep[1]),
           
           bootsim3 <- data.frame(time=rep(predtimes,times=length(simfit$strata)),
                                  group=rep(rep(names(simfit$strata),
                                                each=length(predtimes))),
                                  surv=est,
                                  irep = simobj$irep[1]
           )
    )
    
    return(bootsim3)
  }
  
  
  simPredDist <- map_df(simData, predfun, predtimes=pred.times)
  return(simPredDist)
}
#---------------------------------------------------------------------------

