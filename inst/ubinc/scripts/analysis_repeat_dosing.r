#clearing the workspace
rm(list=ls())
graphics.off()
options(show.error.locations = TRUE)

# If we are in a stand alone ubiquity distribution we run 
# from there otherwise we try to load the package
if(file.exists(file.path('library', 'r_general', 'ubiquity.R'))){
  source(file.path('library', 'r_general', 'ubiquity.R'))
} else { 
  library(ubiquity) }


# Creating the mAb system file
system_new(file_name="system.txt", system_file="mab_pk", overwrite = TRUE)

# For documentation explaining how to modify the commands below
# See the "R Workflow" section at the link below:
# http://presentation.ubiquity.grok.tv

# Rebuilding the system (R scripts and compiling C code)
cfg = build_system()

# set name                  | Description
# -------------------------------------------------------
# default                   | TMDD: Membrane bound target

cfg = system_select_set(cfg, "default")

# fetching the parameter values
parameters = system_fetch_parameters(cfg)


# The following applies to both individual and stochastic simulations:
# Define the solver to use
cfg=system_set_option(cfg,group = "simulation", option = "solver", value = "lsoda")
# Specify the output times 
cfg=system_set_option(cfg, group  = "simulation", 
                           option = "output_times", 
                           seq(0,10*7,1))


# -------------------------------------------------------------------------
# Fixed dosing
cfg = system_zero_inputs(cfg) 
cfg = system_set_bolus(cfg, state   = "Cc", 
                            times   = c(0, 14, 28, 42, 56),         #  days  
                            values  = c(500, 500, 500, 500, 500))   #  mg 


som_fix = run_simulation_ubiquity(parameters, cfg)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Same dosing using titration
# Enabling titration 
cfg=system_set_option(cfg,
                      group       = "titration",
                      option      = "titrate",   
                      value       = TRUE)

# Creating a titration event            
cfg=system_new_tt_rule(cfg,   
                       name       = "ivdose",
                       times      = c(0, 2, 4, 6, 8),   
                       timescale  = "weeks")
                                  
cfg=system_set_tt_cond(cfg,       
                       name       = "ivdose",
                       cond       = "TRUE",   
                       action     = "SI_TT_BOLUS[state='Cc', values=500, times=0]",
                       value      = "1")
                           

som_tt = run_simulation_titrate(parameters, cfg)
# -------------------------------------------------------------------------

 myfig = ggplot() + 
         geom_line(data=som_fix$simout, aes(x=ts.days,   y=Cc, color="Fixed Dosing"), linetype=1) +
         geom_line(data=som_tt$simout,  aes(x=ts.days,   y=Cc, color="Titration"   ), linetype=2)  +
         scale_colour_manual(values=c("Fixed Dosing"="darkblue", "Titration"="firebrick3"))  +
         theme(legend.title = element_blank()) +
         theme(legend.position = 'bottom')     

 myfig = prepare_figure('present', myfig) 

 plot(myfig)

ggsave(sprintf('output%srepeat_dosing.png', .Platform$file.sep), width=8, height=3.4, plot=myfig)

