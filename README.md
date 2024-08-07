# env_reanalysis

This repository has code that calculates environmental parameters and the variability of these parameters within reanalysis-sized grid boxes within LES simulations

1) run_dom_subdom_prof.py --> calculates mean profiles of atmospheric variables within reanalysis sized-boxes and subboxes from INCUS LES simulations
2) calc_env_stats.py --> uses mean profiles from (1) to calculate specific convective environment parameters. This code uses MetPy (May et al., 2024) for some of the calculations and test figures.
3) plot_env_stats.py --> creates and saves statistics for the convective environment parameters calculated in (2) that are then used in the various notebooks for analysis and figures 
4) Figures directory has individual python notebooks that are used to generate analysis and figures for this work. These notebooks use the data that are generated from 1), 2) and 3).
5) Proc_Data includes the statistics of SUBGRID convective environment parameter data created in plot_env_stats.py in pickle files and used in creation of Figures 2, 3, 4 and 5. For example the files look like this: Var_Stats_c0.0001_all_100.p, where c0.0001 represents the condensate threshold and 100 represents the number of grid boxes to ignore at the boundaries (100=10km). Several sensitisivity tests were run in terms of the condensate threshold to use and in terms of shifting the grid (_150, _200). The files that have _120 represent the files run at the end of the hour-long analysis period and the files that include p15 represent the statistics that were calculated using the 16KM_AREAs, as opposed to the 27KM_AREAs. 
6) Proc_Data all includes Reanalysis Proxy values, i.e. mean values averaged over reanalysis-sized grid boxes used for the creation of Figure 5. These are saved in casename directories (i.e., AUS1.1-R/all_100/AUS1.1-R_Environments_ERA5_all_100.p, where AUS1.1-R is the Australia case simulation and all_100 represents the data at the initial time that excludes 10km near the boundaries for the 27KM_AREAs)

References:
May, R. M., Arms, S. C., Marsh, P., Bruning, E., Leeman, J. R., Goebbert, K., Thielen, J. E.,
    Bruick, Z., and Camron, M. D., 2024: MetPy: A Python Package for Meteorological Data.
    Unidata, Unidata/MetPy, doi:10.5065/D6WW7G29.
