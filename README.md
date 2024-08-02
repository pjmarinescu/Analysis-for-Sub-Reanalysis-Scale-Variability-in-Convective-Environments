# env_reanalysis

This repository has code that calculates environmental parameters and the variability of these parameters within reanalysis-sized grid boxes within LES simulations

1) run_dom_subdom_prof.py --> calculates mean profiles of atmospheric variables within reanalysis sized-boxes and subboxes from INCUS LES simulations
2) calc_env_stats.py --> uses mean profiles from (1) to calculate specific convective environment parameters
3) plot_env_stats.py --> creates and saves statistics for the convective environment parameters that are used in the various figure codes 
