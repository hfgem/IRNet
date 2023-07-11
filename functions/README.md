# functions folder
This folder contains functions related to parameter setting, network creation, simulation calculation, and various analyses. It contains a sub-folder of functions related to plotting.

## burst_calculator.m
This function takes in a spike raster and simulation parameters / information and calculates when bursts occur (based on burst parameters), the number of neurons in a burst, the length (in seconds) of a burst, and the inter-burst-interval from the last burst. A structure containing this information is the output.

## create_clusters.m
This function takes in parameters to create a network structure with specified probabilities of excitatory and inhibitory neurons and cluster number / assignment / overlap. A structure containing this information is the output.

## ir_net_calculator_light.m
This function performs the calculations of membrane potentials, conductances, etc... given certain inputs and simulation parameters. Only the spike times (timestep index) of each neuron and the final connectivity matrix are output. Spike times are output as a structure with times per neuron index.

## ir_net_calculator.m
This function performs the calculations of membrane potentials, conductances, etc... given certain inputs and simulation parameters. It keeps track of changes in membrane potential, conductances, and connectivity matrix changes throughout the simulation and outputs matrices containing these values. IT USES A LOT OF MEMORY AND IS BEST USED FOR VERY SHORT SIMULATIONS WITH SMALLER NETWORK SIZES.

## parallelize_parameter_tests.m
This function is to be used in conjunction with misc/param_testing_code.m as it is optimized to be called from a parfor loop. All individual simulation results are saved in a structure called "test_burst_var" and summarized averages for network simulations are stored in the structure "net_burst_results". The specific averages stored are: (1) average number of neurons per burst event, (2) average length of burst events, and (3) average inter-burst-interval.

## set_dependent_parameters.m
This function uses the parameters structure defined in IR_Net.m or param_testing_code.m or any simulation codes in the misc folder. The parameters are set in those files and this function is called to define other parameters or values that are calculated using the defined parameters.