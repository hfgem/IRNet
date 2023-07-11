# misc folder
This folder contains miscellaneous simulation-running code, including parameter testing and simulations of aspects of the code base.

## mock_sim_decay_type.m
This set of code allows the toggling of different connection strength decay modalities (linear, exponential, gaussian) (as well as other parameters) to show how the decay affects connection strengths. Rather than running a full simulation with inputs, this mock simulation has specifically set times when an excitatory neuron fires before an inhibitory neuron, or vice versa, so both sets of immediate connection changes can be seen followed by the specified decay type and decay parameters.

## param_testing_code.m
This set of code allows for parameter testing of any number of parameters at once. The code is written to parallelize model simulations and store outputs for further visualization compared to "success criteria". This code primarily calls on the function "parallelize_parameter_tests.m", which re-packages outputs of simulations in a way that cooperates with parfor loop parallelization standards. All individual simulation results are saved in a structure called "testresults" and summarized averages for network simulations are stored in the structure "netresults". The specific averages stored are: (1) average number of neurons per burst event, (2) average length of burst events, and (3) average inter-burst-interval.
