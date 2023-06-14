%Mock Simulation Code 2
%This code is used to create a mock simulation of I-E spike pairs and
%visualize how the connection strength changes as a function of time

clear all

%% Get file directory and add functions to path

current_file = matlab.desktop.editor.getActiveFilename;
split_dir = split(current_file,'/');
joined_dir = join(split_dir(1:length(split_dir)-2),'/');
current_dir = strcat(joined_dir{1},'/');

cd(current_dir)

addpath('functions')

clear current_file split_dir joined_dir current_dir

%% Simulation options and save path

parameters.saveFlag = 1; % 1 to save simulation results
parameters.selectPath = 1; % 1 to select save destination, 0 to save in current dir
parameters.plotResults = 1; % 1 to plot basic simulation results

%If uploading a parameter file, uncomment the next load line and skip next 2
%sections that deal with parameters:
% load(strcat(save_path,'/parameters.mat'))

%% Parameter Set Up

% Run light simulator
parameters.run_light = 1; %If 1 = only V_m stored; 0 = all parameters stored.

% Network structure parameters
parameters.n = 200; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

% Time
parameters.dt = 1*10^(-3); %timestep (s)
parameters.t_max = 600*parameters.dt; %maximum amount of time (s)

% Basic model parameters
% tau_E ~= 10 ms from direct data, DOI: 10.1126/science.aaf1836
parameters.tau_syn_E = 10*10^(-3); % Exc. synaptic decay time constant (s) PF19=50ms, HF18=10ms for figs 7-8 and longer for earlier figs
% tau_I ~= 1.2-8 ms from direct data, https://doi.org/10.1073/pnas.192233099
parameters.tau_syn_I = 2*10^(-3);  % Inh. synaptic decay time constant (s) PF19=5ms,  HF18=10ms for figs 7-8 and for earlier figs
parameters.tau_stdp = 5*10^(-3); %STDP time constant (s)                 
parameters.E_K = -75*10^(-3); %potassium reversal potential (V) %-75 or -80 mV
parameters.E_L = -70*10^(-3); %leak reversal potential (V) %-60 - -70 mV range
parameters.G_L = 25*10^(-9); %leak conductance (S) %10 - 30 nS range
parameters.C_m = 0.2*10^(-9); %total membrane capacitance (F) %Huge range from 0.1 - 100 pF
parameters.V_m_noise = 10^(-4); % 10^(-4); %magnitude of noise to use in membrane potential simulation (V)
parameters.V_th = -50*10^(-3); %threshold membrane potential (V)
parameters.V_reset = -80*10^(-3); %reset membrane potential (V)
parameters.V_syn_E = 0; %synaptic reversal potential (excitatory)
parameters.V_syn_I = -80*10^(-3); %synaptic reversal potential (inhibitory) %generally -70 pr -80 mV

% Recurrent connection strengths
parameters.del_G_syn_E_E = 5*10^(-10); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_I = 0; %1.4*del_G_syn_E_E; %synaptic conductance step following spike (S)
parameters.del_G_syn_E_I = 5*10^(-10); %synaptic conductance step following spike (S)
parameters.del_G_syn_I_E = 10*10^(-10); %synaptic conductance step following spike (S)

% SRA parameters
parameters.del_G_sra = 200e-09; %spike rate adaptation conductance step following spike %ranges from 1-200 *10^(-9) (S)
parameters.tau_sra = 50*10^(-3); %spike rate adaptation time constant (s)

% STDP parameters
parameters.potentiation_gain_E = 5*10^(-4); %amount to increase connectivity (percent)
parameters.potentiation_gain_I = 1*10^(-4); %amount to increase connectivity (percent)
parameters.depression_gain = 10*10^(-4); %amount to decrease connectivity (percent)
parameters.I_max = 2; %Maximum inhibitory connection strength multiplier
parameters.E_max = 1.5; %Maximum excitatory connection strength multiplier

% Inhibitory gain parameters
parameters.I_E_delay = 10; %Time to delay the sigmoid curve for increased I-E connections (s)

% Input parameters:
% Poisson input
parameters.inputType = 2; % 0 = randn(), 1 = poisson, 2 = theta + randn()
if parameters.inputType == 0
    % Conductance input
    parameters.G_std = 18*10^-9; % STD of the input conductance G_in, if using randn()
    parameters.G_mean = 0*10^-12; % mean of the input conductance G_in, if using randn()
elseif parameters.inputType == 1
    % Poisson input
    parameters.rG = 10; % input spiking rate, if using poisson inputs
    parameters.W_gin = 8*10^-9; % increase in conductance, if using poisson inputs
elseif parameters.inputType == 2
    % Theta input
    parameters.t_freq = 3.5; %Theta frequency (3.5 - 7.5 Hz)
    parameters.t_amp = 9.75*10^(-9); %Theta wave amplitude (tbd a good range)
    parameters.G_std = 0.7*10^(-9); %Theta wave noise std (tbd a good range)
end

% Network connection parameters
parameters.conn_prob = 0.08; %set a total desired connection probability
parameters.p_E = 0.75; %probability of an excitatory neuron
parameters.include_all = 1; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.global_inhib = 0; % if 1, I-cells are not clustered and have connection probability p_I
parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

%Set burst selection parameters
parameters.burst_t_min = 10*10^(-3); %Seconds that must pass without activity to separate bursts

disp("Set Base Parameters")

%Update dependent parameters
parameters = set_dependent_parameters(parameters);

disp("Set Dependent Parameters")
%% Run Simulation + Plot

seed = 1;
network = create_clusters(parameters, 'seed', seed, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);

%Select an I-E pair
I_ind = network.E_indices(randi(length(network.I_indices)));
E_ind = network.E_indices(randi(length(network.E_indices)));

%Copy connectivity matrix
conns = network.conns; %separately update a connectivity matrix

%Store changes
conns_change = zeros(1,parameters.t_steps+1); %Store the I-E connection strength
del_conn_change = zeros(1,parameters.t_steps+1);

%Set I-E connection to 1 and remove reciprocal connection
conns(I_ind,E_ind) = 1;
conns(E_ind,I_ind) = 0;

%Create a maximum connection strength matrix
conns_max = conns;
conns_max(network.E_indices,:) = conns_max(network.E_indices,:)*parameters.E_max;
conns_max(network.I_indices,:) = conns_max(network.I_indices,:)*parameters.I_max;

%Set spike times
V_m = parameters.V_reset*ones(parameters.n,parameters.t_steps+1); %set all neurons to reset value at all times
V_m(I_ind,40) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,45) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,80) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,85) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,120) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,125) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,160) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,165) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,200) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,205) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,240) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,245) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,280) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,285) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,320) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,325) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,360) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,365) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,400) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,405) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,440) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,445) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,480) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,485) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,520) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,525) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike
V_m(I_ind,560) = parameters.V_th + 1*10^(-3); %set inhibitory neuron first spike
V_m(E_ind,565) = parameters.V_th + 1*10^(-3); %set excitatory neuron first spike


%Binary indices of excitatory and inhibitory neurons
E_bin = zeros(parameters.n,1);
E_bin(network.E_indices) = 1;
I_bin = zeros(parameters.n,1);
I_bin(network.I_indices) = 1;
temp_spike_bin = zeros(parameters.n,1); %to use for 

%Variables for STDP
t_spike = zeros(parameters.n,1); %vector to store the time of each neuron's last spike, for use in STDP
t_stdp = parameters.tau_stdp;
pre_strength_E = parameters.potentiation_gain_E;
pre_strength_I = parameters.potentiation_gain_I;
post_strength = parameters.depression_gain;

%Run through each timestep and calculate
for t = 1:parameters.t_steps
    t_curr = t*parameters.dt;
    %check for spiking neurons and postsynaptic and separate into E and I
    spikers_bin = V_m(:,t) >= parameters.V_th; %binary vector
    t_spike(spikers_bin) = t_curr; %update spike times
    spikers_I = spikers_bin.*I_bin; %binary column vector of spiking inhibitory neurons
    spikers_E = spikers_bin.*E_bin; %binary column vector of spiking excitatory neurons
    %______________________________________
    %Update which neurons spiked
    ever_spiked = t_spike > 0;
    E_ever_spiked = ever_spiked.*E_bin; %Binary vector of excitatory ever spiked
    I_ever_spiked = ever_spiked.*I_bin; %Binary vector of inhibitory ever spiked
    pre_syn_n_E = (conns.*(E_ever_spiked*spikers_bin')) > 0; %binary matrix pre-E-n x spikers
    pre_syn_n_I = (conns.*(I_ever_spiked*spikers_bin')) > 0; %binary matrix pre-I-n x spikers
    post_syn_n = (conns.*(spikers_bin*ever_spiked')) > 0; %binary matrix spikers x post-n

    pre_syn_t_E = t_spike.*pre_syn_n_E; %spike times of pre-synaptic excitatory neurons (0 where not pre-synaptic neuron)
    pre_syn_t_I = t_spike.*pre_syn_n_I; %spike times of pre-synaptic inhibitory neurons
    post_syn_t = post_syn_n.*t_spike'; %spike times of all post-synaptic neurons

    t_diff_pre_E = t_curr*pre_syn_n_E - pre_syn_t_E; %time diff between pre-synaptic excitatory neurons and current
    t_diff_pre_I = t_curr*pre_syn_n_I - pre_syn_t_I; %time diff between pre-synaptic inhibitory neurons and current
    t_diff_post = t_curr*post_syn_n - post_syn_t; %time diff between post-synaptic and current
    
    conn_strength_E_pre = (conns_max.*pre_syn_n_E - conns.*pre_syn_n_E)./(conns_max.*pre_syn_n_E);
    conn_strength_E_pre(isnan(conn_strength_E_pre)) = 0;
    conn_strength_I_pre = (conns_max.*pre_syn_n_I - conns.*pre_syn_n_I)./(conns_max.*pre_syn_n_I);
    conn_strength_I_pre(isnan(conn_strength_I_pre)) = 0;
    conn_strength_post = (conns_max.*post_syn_n - conns.*post_syn_n)./(conns_max.*post_syn_n);
    conn_strength_post(isnan(conn_strength_post)) = 0;
    
    del_conn_pre_E = (pre_strength_E/100)*(conn_strength_E_pre).*exp(-t_diff_pre_E/t_stdp).*pre_syn_n_E;
    del_conn_pre_I = (pre_strength_I/100)*(conn_strength_I_pre).*exp(-t_diff_pre_I/t_stdp).*pre_syn_n_I;
    del_conn_post = (post_strength/100)*(conn_strength_post).*exp(-t_diff_post/t_stdp).*post_syn_n;
    
    conns_pre = conns;
    
    conns = conns + del_conn_pre_E + del_conn_pre_I - del_conn_post; %update connections
    
    conns_post = conns;
    
    %Store changes
    conns_change(t) = conns(I_ind,E_ind);
    del_conn_change(t) = conns_post(I_ind,E_ind) - conns_pre(I_ind,E_ind);
    
end

%Plot spikes and connection strength over time
I_spike_inds = find(V_m(I_ind,:) > parameters.V_th,2);
figure;
ax1 = subplot(3,1,1);
plot(parameters.dt*(1:parameters.t_steps+1),V_m(I_ind,:),'DisplayName','Inhibitory Neuron')
hold on
plot(parameters.dt*(1:parameters.t_steps+1),V_m(E_ind,:),'DisplayName','Excitatory Neuron')
title('Membrane Potential')
xlabel('Time (s)')
ylabel('Membrane Potential (V)')
legend()
ax2 = subplot(3,1,2);
plot(parameters.dt*(1:parameters.t_steps+1),conns_change)
hold on
title('Connection Strength')
xlabel('Time (s)')
ylabel('Connection Strength')
legend()
ax3 = subplot(3,1,3);
plot(parameters.dt*(1:parameters.t_steps+1),del_conn_change)
hold on
title('Connection Strength Change')
xlabel('Time (s)')
ylabel('Connection Strength Change')
legend()
