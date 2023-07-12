%IR_Net Project: a model of CTA, latent inhibition, and latent enhancement
%using replay of taste memories and increased inhibition.

clear all

%% Get file directory and add functions to path

current_file = matlab.desktop.editor.getActiveFilename;
split_dir = split(current_file,'/');
joined_dir = join(split_dir(1:length(split_dir)-1),'/');
current_dir = strcat(joined_dir{1},'/');

cd(current_dir)

addpath('functions')
addpath('functions/plotting')

clear current_file split_dir joined_dir current_dir

%% Simulation options and save path

parameters.saveFlag = 1; % 1 to save simulation results
parameters.selectPath = 1; % 1 to select save destination, 0 to save in current dir
parameters.plotResults = 1; % 1 to plot basic simulation results


if parameters.saveFlag & parameters.selectPath
    disp("Select Save Path for Results and Parameters")
    save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored
else
    save_path = [pwd, '/results/irnet'];
end

%If uploading a parameter file, uncomment the next load line and skip next 2
%sections that deal with parameters:
% load(strcat(save_path,'/parameters.mat'))

%% Parameter Set Up

% Run light simulator
parameters.run_light = 1; %If 1 = only spike times and final conns stored; 0 = all parameters stored.

% Network structure parameters
parameters.n = 200; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

% Time
parameters.dt = 1*10^(-3); %timestep (s)
parameters.init_period = 30; %initialization time (s)
parameters.sim_period = 270; %simulation time (s)

% Basic model parameters
% tau_E ~= 10 ms from direct data, DOI: 10.1126/science.aaf1836
parameters.tau_syn_E = 10*10^(-3); % Exc. synaptic decay time constant (s) PF19=50ms, HF18=10ms for figs 7-8 and longer for earlier figs
% tau_I ~= 1.2-8 ms from direct data, https://doi.org/10.1073/pnas.192233099
parameters.tau_syn_I = 2*10^(-3);  % Inh. synaptic decay time constant (s) PF19=5ms,  HF18=10ms for figs 7-8 and for earlier figs
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
parameters.stdp_type = 'decay'; %'growth'; %decay = fast excitatory decay; growth = slow inhibitory growth;
if strcmp(parameters.stdp_type,'decay')
    %Options: 'linear','exponential','gaussian'
    parameters.decay_type = 'linear';
end
parameters.tau_stdp = 5*10^(-3); %STDP time constant (s)                 
if strcmp(parameters.stdp_type,'growth')
    parameters.I_max = 4.5; %Maximum inhibitory connection strength multiplier
    parameters.E_max = 1.5; %Maximum excitatory connection strength multiplier
    parameters.potentiation_gain_E = 5*10^(-2); %amount to increase connectivity (percent)
    parameters.potentiation_gain_I = 1*10^(-2);  %amount to increase connectivity (percent)
elseif strcmp(parameters.stdp_type,'decay')
    parameters.tau_E_decay = -30/log(0.5); %E connection strength decay timescale (s)
    parameters.tau_I_decay = -180/log(0.5); %I connection strength decay timescale (s)
    parameters.potentiation_gain_E = 5*10^(-2); %amount to increase connectivity (percent)
    parameters.potentiation_gain_I = 3*10^(-2); %amount to increase connectivity (percent)
end
parameters.depression_gain = 1*10^(-2); %amount to decrease connectivity (percent)

% Input parameters for conductance input
parameters.inputType = 2; % 0 = randn(), 1 = poisson, 2 = theta + randn()
if parameters.inputType == 0
    % Random noise input - pink noise
    parameters.N_amp = 0.7*10^(-9); %Pink noise max amplitude
elseif parameters.inputType == 1
    % Poisson input
    parameters.rG = 1; % input spiking rate, if using poisson inputs
    parameters.W_gin = 18*10^-9; % increase in conductance, if using poisson inputs
elseif parameters.inputType == 2
    % Theta input
    %   Nitzan et al. 2022 had a criterion of looking at SWR separated by at least
    %   500 ms, which we use for the 2 Hz base oscillation criterion and
    %   randomized added pink noise on top.
    parameters.t_freq = 2; %SWR frequency - Nitzan et al. 2022
    parameters.t_amp = 9.75*10^(-9); %Wave amplitude (tbd a good range)
    parameters.N_amp = 0.7*10^(-9); %Pink noise max amplitude
end

% Network connection parameters
parameters.conn_prob = 0.08; %set a total desired connection probability
parameters.p_E = 0.75; %probability of an excitatory neuron
parameters.include_all = 1; % if a neuron is not in any cluster, take cluster membership from a highly connected neuron
parameters.global_inhib = 0; % if 1, I-cells are not clustered and have connection probability p_I
parameters.p_I = 0.5; % probability of an I cell connecting to any other cell

%Set burst selection parameters
parameters.burst_n_min = 0.05; %Fraction of neurons that must be active in a burst for it to count
parameters.burst_t_min = 10*10^(-3); %Seconds that must pass without activity to separate bursts
parameters.num_burst_avg = 5; %Number of bursts to average together for visualizations
parameters.num_rast_to_plot = 25; %Number of burst rasters to plot

%Set conductance test parameters
parameters.pre_spike = 1; %Seconds before spike to visualize conductance inputs
parameters.post_spike = 0.5; %Seconds after spike to visualize conductance inputs

% Number of trials per net to run
parameters.nTrials = 1; % How many tests of different initializations to run
parameters.nNets = 1; % How many networks to run

disp("Set Base Parameters")

% Update dependent parameters and save
parameters = set_dependent_parameters(parameters);

disp("Set Dependent Parameters")

%Save to computer
if parameters.saveFlag & ~isfolder(save_path)
    mkdir(save_path);
end
    
if parameters.saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters','-v7.3'); 
    disp("Saved Parameters")
end

%% Run simulation

for ithNet = 1:parameters.nNets
    disp('Network ' + string(ithNet))
    %CREATE NETWORK SAVE PATH
    net_save_path = strcat(save_path,'/network_',string(ithNet));
    if parameters.saveFlag & ~isfolder(net_save_path)
        mkdir(net_save_path);
    end
    
    % Generete ith network structure
    seed = ithNet;
    network = create_clusters(parameters, 'seed', seed, 'include_all', parameters.include_all, 'global_inhib', parameters.global_inhib);
    if parameters.saveFlag
        save(strcat(net_save_path,'/network.mat'),'network');
    end
    
    %RUN MODEL AND CALCULATE
    %Run through every cluster initialization and store relevant data and
    %calculations
    burst_var = struct;
    if parameters.run_light ~= 1
        V_m_var = struct;
        G_var = struct;
        I_var = struct;
        %conns_var = struct;
    end
    
    for ithTest = 1:parameters.nTrials       
        disp('Simulation ' + string(ithNet))
        %Create input conductance variable
        if parameters.inputType == 0
            % Pink noise input
            try %With installed audio package
                pink_noise = pinknoise(parameters.t_steps+1,parameters.n);
                pink_noise_scaled = parameters.N_amp*(pink_noise./max(pink_noise));
            catch %For older MATLAB versions or uninstalled audio package
                pink_noise = dsp.ColoredNoise('Color','pink','NumChannels',parameters.n,'SamplesPerFrame',parameters.t_steps+1);
                pink_noise_scaled = parameters.N_amp*(pink_noise()./max(pink_noise()));
            end
            G_in = pink_noise_scaled;
            G_in(G_in<0) = 0;
        elseif parameters.inputType == 1
            % Poisson input
            G_in = zeros(parameters.n, parameters.t_steps+1);
            for k = 2:(parameters.t_steps+1)
                G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
                G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
            end
        elseif parameters.inputType == 2
            % Theta input + pink noise
            try %With installed audio package
                pink_noise = pinknoise(parameters.t_steps+1,parameters.n);
                pink_noise_scaled = parameters.N_amp*(pink_noise./max(pink_noise));
            catch %For older MATLAB versions or uninstalled audio package
                pink_noise = dsp.ColoredNoise('Color','pink','NumChannels',parameters.n,'SamplesPerFrame',parameters.t_steps+1);
                pink_noise_scaled = parameters.N_amp*(pink_noise()./max(pink_noise()));
            end
            theta_wave = (parameters.t_amp/2)*sin(2*pi*parameters.t_freq*(parameters.dt*ones(parameters.n,parameters.t_steps+1).*(1:parameters.t_steps+1))) + parameters.t_amp/2;
            G_in = theta_wave + (pink_noise_scaled)';
            G_in(G_in<0) = 0;
        end
        parameters.('G_in') = G_in;
        
        seed = ithTest;
        
        %Run model
        if parameters.run_light == 1
            V_m = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
            [spike_struct, last_conns] = ir_net_calculator_light(parameters, seed, network, V_m);
            spikes_bin = zeros(parameters.n, parameters.t_steps+1);
            for s_i = 1:parameters.n
               spikes_bin(s_i,spike_struct(s_i).times) = 1; 
            end
            conns = zeros(parameters.n,parameters.n,2); %save connectivity changes
            conns(:,:,1) = network.conns;
            conns(:,:,2) = last_conns;
            clear last_conns
        else
            %Create Storage Variables
            V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
            V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
            [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = ir_net_calculator(...
                parameters, seed, network, V_m);
            V_m_var(ithTest).V_m = V_m;
            G_var(ithTest).G_in = G_in;
            G_var(ithTest).G_sra = G_sra;
            G_var(ithTest).G_syn_I_E = G_syn_I_E;
            G_var(ithTest).G_syn_E_E = G_syn_E_E;
            G_var(ithTest).G_syn_I_I = G_syn_I_I;
            G_var(ithTest).G_syn_E_I = G_syn_E_I;
            clear V_m G_sra G_syn_E_E G_syn_I_E G_syn_E_I G_syn_I_I G_in
            spikes_bin = V_m > parameters.V_th;
            %conns_var(ithTest).conns = conns;
        end
        
        %Plot raster of full simulation
        if parameters.plotResults == 1
            plot_rast(parameters,spikes_bin,ithTest,net_save_path)
        end
        
        %Calculate bursts
        bursts = burst_calculator(parameters,spikes_bin,net_save_path,ithNet,ithTest);
        burst_var(ithTest).bursts = bursts;
        clear bursts
        
        %Plot connection strength changes
        if parameters.plotResults == 1
            plot_conns(parameters,conns,ithNet,ithTest,net_save_path)
        end
        
    end %end test loop
    
    %SAVE NETWORK DATA
    if parameters.saveFlag
        save(strcat(net_save_path,'/burst_var.mat'),'burst_var','-v7.3')
        if parameters.run_light ~= 1
            save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
            save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
            save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
            %save(strcat(net_save_path,'/conns_var.mat'),'conns_var','-v7.3')
        end
    end

end %end network loop
disp("Simulations Done")

%% Look at connection changes

flat_all_conns = reshape(conns(:,:,2) - conns(:,:,1),[1,parameters.n^2]);
flat_all_conns_no_0 = flat_all_conns(flat_all_conns~=0);
flat_E_conns = reshape(conns(network.E_indices,:,2) - conns(network.E_indices,:,1),[1,length(network.E_indices)*parameters.n]);
flat_E_conns_no_0 = flat_E_conns(flat_E_conns~=0);
flat_I_conns = reshape(conns(network.I_indices,:,2) - conns(network.I_indices,:,1),[1,length(network.I_indices)*parameters.n]);
flat_I_conns_no_0 = flat_I_conns(flat_I_conns~=0);

figure;
histogram(flat_all_conns_no_0,'DisplayName','All connections')
hold on
histogram(flat_E_conns_no_0,'DisplayName','E connections')
histogram(flat_I_conns_no_0,'DisplayName','I connections')
legend()
ylabel("Number of Occurrences")
xlabel("Change in Connection Strength")
title("Change in Connection Strength (Post - Pre)")

%% Visualize the effect of I-E delayed strengthening

%First find I-E pairings
I_indices = network.I_indices;
E_indices = network.E_indices;
[I_ind,E_ind] = ind2sub([length(I_indices),length(E_indices)],find(conns(I_indices,E_indices,parameters.t_steps + 1) > 0));
I_E_pairs = [I_indices(I_ind);E_indices(E_ind)]';

%Now look at the average connection strength over the course of a
%simulation
I_E_conn_changes = zeros(length(I_E_pairs),parameters.t_steps+1);
for i = 1:length(I_E_pairs)
    I_E_conn_changes(i,:) = squeeze(conns(I_E_pairs(i,1),I_E_pairs(i,2),:));
end

figure;
p_i = 1;
ax1 = subplot(2,1,1);
plot(parameters.dt*(1:parameters.t_steps+1),I_E_conn_changes(p_i,:))
ax2 = subplot(2,1,2);
plot(parameters.dt*(1:parameters.t_steps+1),(V_m(I_E_pairs(p_i,1),:) > parameters.V_th),'DisplayName','inhibitory neuron')
hold on
plot(parameters.dt*(1:parameters.t_steps+1),(V_m(I_E_pairs(p_i,2),:) > parameters.V_th),'DisplayName','excitatory neuron')
legend()
linkaxes([ax1,ax2],'x')

