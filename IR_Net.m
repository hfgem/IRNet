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
parameters.run_light = 1; %If 1 = only V_m stored; 0 = all parameters stored.

% Network structure parameters
parameters.n = 200; %number of neurons
parameters.clusters = 10; % Number of clusters in the network
parameters.mnc = 2; % mean number of clusters each neuron is a member of

% Time
parameters.t_max = 60; %maximum amount of time (s)
parameters.dt = 1*10^(-3); %timestep (s)

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
parameters.potentiation_gain = 5*10^(-4); %amount to increase connectivity
parameters.depression_gain = 10*10^(-4); %amount to decrease connectivity (at least ~2x potentiation is good)

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

% Number of trials per net to run
parameters.nTrials = 2; % How many tests of different initializations to run
parameters.nNets = 3; % How many networks to run

disp("Set Base Parameters")

%% Update dependent parameters and save

parameters = set_dependent_parameters(parameters);

disp("Set Dependent Parameters")

%Save to computer
if parameters.saveFlag & ~isfolder(save_path)
    mkdir(save_path);
end
    
if parameters.saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters'); 
    disp("Saved Parameters")
end

%% Run simulation
%close all

for ithNet = 1:parameters.nNets
    
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
    V_m_var = struct;
    if parameters.run_light ~= 1
        G_var = struct;
        I_var = struct;
        burst_var = struct;
        %conns_var = struct;
    end
    
    for ithTest = 1:parameters.nTrials       
        
        %Create input conductance variable
        if parameters.inputType == 0
            % Conductance input
            G_in = (parameters.G_std*randn(parameters.n,parameters.t_steps+1))+parameters.G_mean;
            G_in(G_in<0) = 0;
        elseif parameters.inputType == 1
            % Poisson input
            G_in = zeros(parameters.n, parameters.t_steps+1);
            for k = 2:(parameters.t_steps+1)
                G_in(:,k) = G_in(:,k-1)*exp(-parameters.dt/parameters.tau_syn_E);
                G_in(:,k) = G_in(:,k) + parameters.W_gin * [rand(parameters.n, 1) < (parameters.dt*parameters.rG)];
            end
        elseif parameters.inputType == 2
            % Theta input
            G_in = (parameters.t_amp/2)*sin(2*pi*parameters.t_freq*(parameters.dt*ones(parameters.n,parameters.t_steps+1).*(1:parameters.t_steps+1))) + (parameters.G_std*randn(parameters.n,parameters.t_steps+1)) + (parameters.t_amp/2);
            G_in(G_in<0) = 0;
        end
        parameters.('G_in') = G_in;
        
        %Create Storage Variables
        V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
        V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
        
        seed = ithTest;
        
        %Run model
        if parameters.run_light == 1
            [V_m, last_conns] = ir_net_calculator_light(parameters, seed, network, V_m);
            V_m_var(ithTest).V_m = V_m;
        else
            [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = ir_net_calculator(...
                parameters, seed, network, V_m);
            V_m_var(ithTest).V_m = V_m;
            G_var(ithTest).G_in = G_in;
            G_var(ithTest).G_sra = G_sra;
            G_var(ithTest).G_syn_I_E = G_syn_I_E;
            G_var(ithTest).G_syn_E_E = G_syn_E_E;
            G_var(ithTest).G_syn_I_I = G_syn_I_I;
            G_var(ithTest).G_syn_E_I = G_syn_E_I;
            %conns_var(ithTest).conns = conns;
        end
        
        %Calculate bursts
        bursts = burst_calculator(parameters,V_m);
        burst_var(ithTest).bursts = bursts;
        
        %Calculate average burst length as a function of time
        num_burst_avg = 5;
        num_chunks = ceil(length(bursts)/num_burst_avg);
        burst_avg_len = zeros(1,num_chunks);
        burst_avg_sz = zeros(1,num_chunks);
        for n_i = 1:num_chunks
            burst_ind = max([0,(n_i-1)*num_burst_avg + 1]);
            burst_end_ind = min([length(bursts),burst_ind + num_burst_avg]);
            burst_avg_len(n_i) = mean([bursts(burst_ind:burst_end_ind).length]);
            burst_avg_sz(n_i) = mean([bursts(burst_ind:burst_end_ind).num_neur]);
        end    
        
        %REMOVE FOLLOWING LINES BEFORE END LATER
        disp('Number of Spikes: ' + string(sum(V_m > parameters.V_th,'all')))
%         figure; 
%         imagesc(V_m > parameters.V_th)
%         xlabel('Time (in dt)')
%         ylabel('Neuron Index')
%         title('Spike Visual')
        figure;
        post_delay_ind = find([bursts.pre_delay] == 1,1);
        subplot(3,2,1)
        plot([bursts.length])
        hold on
        xline(post_delay_ind)
        title('Burst Lengths')
        subplot(3,2,3)
        plot([bursts.num_neur])
        hold on
        xline(post_delay_ind)
        title('Burst Size (Number of Neurons)')
        subplot(3,2,5)
        plot([bursts.ibi])
        hold on
        xline(post_delay_ind)
        title('Time to Next Burst')
        subplot(3,2,2)
        plot(1:num_chunks,burst_avg_len)
        hold on
        xline(post_delay_ind/num_burst_avg)
        title('Average Burst Length for ' + string(num_burst_avg) + ' Bursts')
        subplot(3,2,4)
        plot(1:num_chunks,burst_avg_sz)
        hold on
        xline(post_delay_ind/num_burst_avg)
        title('Average Burst Size for ' + string(num_burst_avg) + ' Bursts')
%         if parameters.run_light == 1
%             figure;
%             ax1 = subplot(1,2,1);
%             imagesc(network.conns)
%             title('Starting Connectivity')
%             colorbar()
%             ax2 = subplot(1,2,2);
%             imagesc(last_conns)
%             title('Ending Connectivity')
%             colorbar()
%             linkaxes([ax1,ax2])
%         else
%             figure;
%             ax1 = subplot(1,2,1);
%             imagesc(squeeze(conns(:,:,1)))
%             title('Starting Connectivity')
%             colorbar()
%             ax2 = subplot(1,2,2);
%             imagesc(squeeze(conns(:,:,parameters.t_steps + 1)))
%             title('Ending Connectivity')
%             colorbar()
%             linkaxes([ax1,ax2])
%         end
        
    end %end test loop
    
    %SAVE NETWORK DATA
    if parameters.saveFlag
        save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
        if parameters.run_light ~= 1
            save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
            save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
            save(strcat(net_save_path,'/burst_var.mat'),'burst_var','-v7.3')
            %save(strcat(net_save_path,'/conns_var.mat'),'conns_var','-v7.3')
        end
    end

end %end network loop
disp("Simulations Done")

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

