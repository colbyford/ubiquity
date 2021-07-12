library(azuremlsdk)
library(ubiquity)
library(tidyverse)
library(here)
library(optparse)


## Parameters --------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
cat(args)

data_folder <- as.character(args[2])

scen <- as.character(args[4])
log_metric_to_run("scenario", scen)

print(paste("Scenario: ", scen))


cfg = build_system(file.path("analysis", "system.txt"))

# This pulls out the default paramters, you can overwrite 
# them in the scnenarios below
parameters = system_fetch_parameters(cfg)


#-----------------------------------
# Change analysis options here:
analysis_name = "analysis_short_name"

# Run simulations for scenarios not currently available 
runsims       = TRUE 
# If runsims is TRUE and overwrite is TRUE then previously 
# run simulations will be overwritten. If overwrite is FALSE the
# previously run simulations will be skipped
overwrite     = FALSE

# Switches between 1 and n-1 cores. 
multicore     = TRUE

# This is the final simulation time in days. It can be set for individual
# scenarios below if needed
tfinal        = 3*10*7

# Simulated output will contain predefined times to make smooth profiles, but
# you may want to force specific times to be included. 
# Put the output times you want in the simulated output here. 
# Set it to NULL to ignore this option
tinclude      = c(10.5, 30000)

# Number of subjects to generate. These are the population 
# that will be sampled from when performing the Monte 
# Carlo simulations
nsub_gen  = 500

# Each trial will consist of nsub subjects sampled from 
# those generated above
nsub      = 25 

# This is the number of trials to run:
ntrial    = 25 

# For Metrums analysis
# nsub_gen  = 500
# nsub      = 250
# ntrial    = 250 

scens = list()

scens = list(
  S01= list(DREG      = "0.8 mg/kg QW (21 day cycle)",    # Unique description for the current scenario
            pvals     = parameters,                       # Typical parameter values
            nsub_gen  = nsub_gen,                         # Total number of subjects to generate
            nsub      = nsub,                             # Number of subjects to sample per trial
            ntrial    = ntrial,                           # Number of trials to simulate
            cycle_dur = 21,                               # Cycle duration (days)
            ncycles   = 10,                               # Number of treatment cycles (#)
            tfinal    = tfinal,                           # Final simulation time in days
            dtimes    = c(  0,   7,  14),                 # Dose times (days) relative to start of cycle
            dvals     = c(0.8, 0.8, 0.8)),                # Corresponding dose values (mg/kg)
  S02= list(DREG      = "1.2 mg/kg QW (21 day cycle)",
            pvals     = parameters, 
            nsub_gen  = nsub_gen,
            nsub      = nsub,
            ntrial    = ntrial,
            cycle_dur = 21,
            ncycles   = 10,
            tfinal    = tfinal,  
            dtimes    = c(  0,   7,  14),
            dvals     = c(1.2, 1.2, 1.2)))


source(file.path("resources", "pn_funcs.R"))

pnall = NULL
piall = NULL
dsall = NULL

data_file = paste(analysis_name, "_pn_scens-", scen, ".RData", sep="")
data_file_path = file.path("outputs", data_file)

pn = NULL
# If we're running simulations then we run them here
if(runsims){
  if(!file.exists(data_file_path) | overwrite){
    vp(cfg, paste("Generating subjects for scneario", scen))
    # First we generate subjects for the trial
    SUBS = gen_subs(cfg         = cfg, 
                    parameters  = scens[[scen]]$pvals,
                    nsub        = scens[[scen]]$nsub_gen)
    
    # Loading the dataset into the system object
    cfg = system_load_data(cfg, dsname     = "SUBS", 
                           data_file  =  SUBS)
    pn = sim_pn(cfg, 
                nsub      =  scens[[scen]]$nsub,       #  number of subject to sample for each trial 
                ntrial    =  scens[[scen]]$ntrial,     #  number of trials to simulate
                dsname    = "SUBS",                    #  name of the dataset to sample from
                multicore = multicore,                 #  Boolean variable indicating if you want to use multiple cores (TRUE) or a single core (FALSE)
                scenario  = scens[[scen]]$DREG,        #  character description of the dosing regimen, parametric perturbations, etc
                tfinal    = scens[[scen]]$tfinal,      #  final time to simulate PN (days)
                tinclude  = tinclude,                  #  output times to include in the simulated output
                dtimes    = scens[[scen]]$dtimes,      #  dosing times (days)
                dvals     = scens[[scen]]$dvals,       #  corresponding dose values (mg/kg)
                inf_dur   = 0.5,                       #  infusion duration (hours)
                cycle_dur = scens[[scen]]$cycle_dur,   #  duration of treatment cycle (days)
                ncycles   = scens[[scen]]$ncycles)     #  number of cycles (#)
    
    save(pn, SUBS, file=data_file_path)
  }
} 

# If we're not running simulations we'll try to load previously stored data
if(file.exists(data_file_path)){
  vp(cfg, paste("loading file:", data_file_path))
  load(file=data_file_path)
} else {
  vp(cfg, paste("No file:", data_file_path))
  vp(cfg, paste("Skipping scenario:", scen)) 
}

# Compile the results into a dataset
pnall = rbind(pnall, pn$pn_tc)
piall = rbind(piall, pn$pi_tc)
dsall = rbind(dsall, pn$dosing_summary)


log_metric_to_run("scen_output", 1)

# log_table_to_run("pnall",
#                  as.list(pnall))
# log_table_to_run("piall",
#                  as.list(piall))
# log_table_to_run("dsall",
#                  as.list(dsall))