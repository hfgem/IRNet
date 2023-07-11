# plotting folder
This folder contains functions related to plotting.

## plot_burst_stats.m
This function plots the statistics of burst events stored in the structure "bursts", output by the function burst_calculator.m. The figure shows burst lengths, burst sizes, and inter-burst-intervals both over the course of a simulation and averaged across the simulation by some set number of bursts (say 5 bursts at a time).

## plot_bursts.m
This function plots individual burst rasters from the structure "bursts".

## plot_conns.m
This function plots connection strength changes from start to finish of a simulation as found in the 3D matrix "conns".

## plot_param_test_results.m
This function plots the results of parameter testing where pairs of varied parameters are plotted against each other using average values as well as binary "good parameter" delineations based on burst statistics.

## plot_rast.m
This function plots the full simulation raster found in "spikes_bin".