# env_reanalysis

This repository has code that calculates environmental parameters and the variability of these parameters within reanalysis-sized grid boxes within LES simulations

1) run_dom_subdom_prof.py --> calculates mean profiles of atmospheric variables within reanalysis sized-boxes and subboxes from INCUS LES simulations
2) calc_env_stats.py --> uses mean profiles from (1) to calculate specific convective environment parameters. This code uses MetPy (May et al., 2024) for some of the calculations and test figures.
3) plot_env_stats.py --> creates and saves statistics for the convective environment parameters calculated in (2) that are then used in the various notebooks for analysis and figures 
4) Figures directory has individual python notebooks that are used to generate analysis and figures for this work. These notebooks use the data that are generated from 1), 2) and 3).

References:
May, R. M., Arms, S. C., Marsh, P., Bruning, E., Leeman, J. R., Goebbert, K., Thielen, J. E.,
    Bruick, Z., and Camron, M. D., 2024: MetPy: A Python Package for Meteorological Data.
    Unidata, Unidata/MetPy, doi:10.5065/D6WW7G29.
