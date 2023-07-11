function [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conn_changes] = ir_net_calculator(...
    parameters, seed, network, V_m)
    %_________
    %ABOUT: This function uses the leaky integrate-and-fire model of 
    %neuronal firing to calculate the trajectory of membrane potentials,
    %currents, etc... that take place in a particular network with a 
    %particular set of parameters and initialization.
    %
    %INPUTS:
    %   parameters = a structure that contains a majority of LIF neuron
    %   initialization parameters. A few relevant listed below:
    %       n = Number of neurons in the network
    %       clusters = Number of clusters of neurons in network
    %       t_max = maximum time of simulation (s)
    %       dt = timestep of simulation (s)
    %       t_steps = The number of timesteps in the simulation
    %       V_th = The threshold membrane potential (V)
    %       V_reset = The reset membrane potential (V)
    %       connectivity_gain = Amount to increase or decrease connectivity by 
    %               with each spike (more at the range of 1.002-1.005) -
    %               keep at 1 to ensure no connectivity change
    %       G_in = input conductance matrix of size [n,t_steps+1]
    %   seed = A random number generator seed which affects membrane
    %       potential noise
    %   network = a structure that contains the following:
    %       cluster_mat = A binary [clusters x n] matrix of which neurons are
    %               in which cluster
    %       conns = An [n x n] matrix of which neurons are connected to each
    %               other, with values greater than 1 implying stronger
    %               connectivity
    %       I_indices = Vector of indices of inhibitory neurons
    %       E_indices = Vector of indices of excitatory neurons
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_in = An [n x t_steps+1] matrix of input conductance for each
    %               neuron at each timestep
    %
    %OUTPUTS:
    %   V_m = An [n x t_steps+1] matrix of membrane potential for each 
    %               neuron at each timestep
    %   G_sra = An [n x t_steps+1] matrix with refractory conductance for
    %               each neuron at each timestep (S)
    %   G_syn_E_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory to postsynaptic excitatory (S)
    %   G_syn_I_E = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory to postsynaptic excitatory (S)
    %   G_syn_E_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               excitatory to postsynaptic inhibitory (S)
    %   G_syn_I_I = An [n x t_steps+1] matrix of conductance for presynaptic 
    %               inhibitory to postsynaptic inhibitory (S)
    %   conns = The updated conns matrix (if connectivity_gain != 0)
    %
    %ASSUMPTIONS:
    %   1. LIF model with spike rate adaptation and synaptic transmission
    %   2. SRA has decay
    %   3. Synaptic transmission has decay
    %   4. There can be an externally applied current input through 
    %       I_in
    %   5. Excitatory and inhibitory connections have different reversal
    %       potentials in the postsynaptic neuron represented in vectors
    %       E_syn_E and E_syn_I
    %   6. Excitatory and inhibitory currents have different synaptic
    %       conductance steps and decay rates
    %_________
    
    %Set the random number generator seed
    rng(seed)
    
    %Create Storage Variables
    G_sra = zeros(parameters.n,parameters.t_steps+1); %refractory conductance for each neuron at each timestep (S)
    G_syn_I_E = zeros(parameters.n,parameters.t_steps+1); %conductance for pre-inhib to post-excit (S)
    G_syn_E_E = zeros(parameters.n,parameters.t_steps+1); %conductance for pre-excit to post-excit (S)
    G_syn_I_I = zeros(parameters.n,parameters.t_steps+1); %conductance for pre-inhib to post-inhib (S)
    G_syn_E_I = zeros(parameters.n,parameters.t_steps+1); %conductance for pre-excit to post-inhib (S)

    %Copy connectivity matrix in case of stdp changes
    conns = network.conns; %separately update a connectivity matrix
    conn_changes = zeros(parameters.n,parameters.n,parameters.t_steps+1); %save connectivity changes
    conn_changes(:,:,1) = conns;
    
    %Binary indices of excitatory and inhibitory neurons
    E_bin = zeros(parameters.n,1);
    E_bin(network.E_indices) = 1;
    I_bin = zeros(parameters.n,1);
    I_bin(network.I_indices) = 1;
    
    %Store a maximum connection strength matrix
    if strcmp(parameters.stdp_type,'growth')
        conns_max = conns;
        conns_max(network.E_indices,:) = conns_max(network.E_indices,:)*parameters.E_max;
        conns_max(network.I_indices,:) = conns_max(network.I_indices,:)*parameters.I_max; %#ok<*NASGU>
    end
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
        t_spike(spikers_bin) = t_curr;
        spikers_I = spikers_bin.*I_bin; %binary column vector of spiking inhibitory neurons
        spikers_E = spikers_bin.*E_bin; %binary column vector of spiking excitatory neurons
        %______________________________________
        %Adjust parameters dependent on spiking
        G_sra(spikers_bin,t) = G_sra(spikers_bin,t) + parameters.del_G_sra; %set SRA conductance values
        %Synaptic conductance is stepped for postsynaptic neurons
        %   dependent on the number of presynaptic connections, and the
        %   current will depend on the presynaptic neuron type (E_syn_I and E_syn_E)
        incoming_conn_E = sum(conns.*(spikers_E*ones(1,parameters.n)),1)'; %post-synaptic neuron E input counts
        incoming_conn_I = sum(conns.*(spikers_I*ones(1,parameters.n)),1)'; %post-synaptic neuron I input counts
        G_syn_I_E(:,t) = G_syn_I_E(:,t) + parameters.del_G_syn_I_E*incoming_conn_I.*E_bin;
        G_syn_E_E(:,t) = G_syn_E_E(:,t) + parameters.del_G_syn_E_E*incoming_conn_E.*E_bin;
        G_syn_I_I(:,t) = G_syn_I_I(:,t) + parameters.del_G_syn_I_I*incoming_conn_I.*I_bin;
        G_syn_E_I(:,t) = G_syn_E_I(:,t) + parameters.del_G_syn_E_I*incoming_conn_E.*I_bin;
        %______________________________________
        %Calculate membrane potential using integration method
        V_ss = ( parameters.G_in(:,t).*parameters.syn_E + G_syn_E_E(:,t).*parameters.syn_E + ...
            G_syn_E_I(:,t).*parameters.syn_E + G_syn_I_I(:,t).*parameters.syn_I + ...
            G_syn_I_E(:,t).*parameters.syn_I + parameters.G_L*parameters.E_L + ...
            G_sra(:,t)*parameters.E_K )./(parameters.G_L + G_sra(:,t) + G_syn_E_E(:,t) + ...
            G_syn_E_I(:,t) + G_syn_I_I(:,t) + G_syn_I_E(:,t) + parameters.G_in(:,t));
        taueff = parameters.C_m./(parameters.G_L + G_sra(:,t) + G_syn_E_E(:,t) + ...
            G_syn_E_I(:,t) + G_syn_I_I(:,t) + G_syn_I_E(:,t) + parameters.G_in(:,t));
        V_m(:,t+1) = V_ss + (V_m(:,t) - V_ss).*exp(-parameters.dt ./taueff) + ...
            randn([parameters.n,1])*parameters.V_m_noise*sqrt(parameters.dt); %the randn portion can be removed if you'd prefer no noise
        V_m(spikers_bin,t+1) = parameters.V_reset; %update those that just spiked to reset
        %______________________________________
        %Update next step conductances
        G_sra(:,t+1) = G_sra(:,t)*exp(-parameters.dt/parameters.tau_sra); %Spike rate adaptation conductance
        %Synaptic conductance updated for each postsynaptic neuron by
        %incoming connection type
        G_syn_E_E(:,t+1) = G_syn_E_E(:,t).*exp(-parameters.dt/parameters.tau_syn_E); %excitatory conductance update
        G_syn_I_E(:,t+1) = G_syn_I_E(:,t).*exp(-parameters.dt/parameters.tau_syn_I); %excitatory conductance update
        G_syn_I_I(:,t+1) = G_syn_I_I(:,t).*exp(-parameters.dt/parameters.tau_syn_I); %inhibitory conductance update
        G_syn_E_I(:,t+1) = G_syn_E_I(:,t).*exp(-parameters.dt/parameters.tau_syn_E); %inhibitory conductance update
        %______________________________________
        %Update connection strengths with decay
        if strcmp(parameters.stdp_type,'decay')
            %Options: 'linear','exponential','gaussian'
            if strcmp(parameters.decay_type,'linear')
                conns(network.E_indices,:) = conns(network.E_indices,:) + -parameters.dt/parameters.tau_E_decay;
                conns(network.I_indices,:) = conns(network.I_indices,:) + -parameters.dt/parameters.tau_I_decay;
            elseif strcmp(parameters.decay_type,'exponential')    
                conns(network.E_indices,:) = conns(network.E_indices,:)*exp(-parameters.dt/parameters.tau_E_decay);
                conns(network.I_indices,:) = conns(network.I_indices,:)*exp(-parameters.dt/parameters.tau_I_decay);
            elseif strcmp(parameters.decay_type,'gaussian')
                conns(network.E_indices,:) = conns(network.E_indices,:)*exp(-(parameters.dt)^2/parameters.tau_E_decay);
                conns(network.I_indices,:) = conns(network.I_indices,:)*exp(-(parameters.dt)^2/parameters.tau_I_decay);
            end
        end
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
        if strcmp(parameters.stdp_type,'growth')
            conn_strength_E_pre = (conns_max.*pre_syn_n_E - conns.*pre_syn_n_E)./(conns_max.*pre_syn_n_E);
            conn_strength_E_pre(isnan(conn_strength_E_pre)) = 0;
            conn_strength_I_pre = (conns_max.*pre_syn_n_I - conns.*pre_syn_n_I)./(conns_max.*pre_syn_n_I);
            conn_strength_I_pre(isnan(conn_strength_I_pre)) = 0;
            conn_strength_post = (conns_max.*post_syn_n - conns.*post_syn_n)./(conns_max.*post_syn_n);
            conn_strength_post(isnan(conn_strength_post)) = 0;
            del_conn_pre_E = (pre_strength_E/100)*(conn_strength_E_pre).*exp(-t_diff_pre_E/t_stdp).*pre_syn_n_E;
            del_conn_pre_I = (pre_strength_I/100)*(conn_strength_I_pre).*exp(-t_diff_pre_I/t_stdp).*pre_syn_n_I;
            del_conn_post = (post_strength/100)*(conn_strength_post).*exp(-t_diff_post/t_stdp).*post_syn_n;
        elseif strcmp(parameters.stdp_type,'decay')
            del_conn_pre_E = (pre_strength_E/100)*exp(-t_diff_pre_E/t_stdp).*pre_syn_n_E;
            del_conn_pre_I = (pre_strength_I/100)*exp(-t_diff_pre_I/t_stdp).*pre_syn_n_I;
            del_conn_post = (post_strength/100)*exp(-t_diff_post/t_stdp).*post_syn_n;
        end
        %Store changes
        conns = conns + del_conn_pre_E + del_conn_pre_I - del_conn_post; %update connections
 
        %Save new connectivity matrix
        conn_changes(:,:,t+1) = conns;
    end