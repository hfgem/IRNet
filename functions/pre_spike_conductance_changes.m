function [neur_conductances] = pre_spike_conductance_changes(...
    parameters, network, V_m, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I)
    %_________
    %ABOUT: This function pulls out input conductances a set amount of time 
    %prior to a spike from both intra-network activity conductances and
    %external conductance input. It calls another function that creates a
    %plot of the different types of conductances prior to a spike, as well
    %as the average waveform of conductance inputs prior to a neuron
    %spiking. In addition, it creates an output structure "change_results"
    %which contains the conductance calculations per spike.
    %
    %INPUTS:
    %   parameters: structure containing simulation parameters
    %   network: structure containing network structure information
    %   spikes_bin: binary spike matrix containing values for each neuron
    %       at each point in time in the simulation. (n x t_steps+1)
    %   G_syn_E_E: matrix containing E to E conductances (n x t_steps+1)
    %   G_syn_I_E: matrix containing I to E conductances (n x t_steps+1)
    %   G_syn_E_I: matrix containing E to I conductances (n x t_steps+1)
    %   G_syn_I_I: matrix containing I to I conductances (n x t_steps+1)
    %
    %OUTPUTS:
    %   change_results: structure containing 
    %_________
    
    %For each neuron, pull out time intervals pre- and post- spikes
    %Grab membrane potential and conductance information as well
    spikes_bin = V_m >= parameters.V_th;
    sim_ind = ceil(parameters.init_period/parameters.dt);
    spikes_bin_cropped = spikes_bin(:,sim_ind:end);
    sim_end = parameters.t_steps+1-sim_ind;
    clear spikes_bin
    pre_buffer = ceil(parameters.pre_spike/parameters.dt);
    post_buffer = ceil(parameters.post_spike/parameters.dt);
    neur_conductances = struct; %store intervals here for each neuron
    for n_i = 1:parameters.n
        spike_t = find(spikes_bin_cropped(n_i,pre_buffer:end-post_buffer));
        pre_t = spike_t;
        post_t = spike_t + pre_buffer + post_buffer;
        neur_conductances(n_i).ranges = [pre_t;post_t]';
        %Pull out the membrane potential waveform for each neuron
        V_m_i = [];
        G_syn_E_E_i = [];
        G_syn_E_I_i = [];
        G_syn_I_E_i = [];
        G_syn_I_I_i = [];
        for r_i = 1:length(spike_t)
            start_i = sim_ind+pre_t(r_i);
            end_i = sim_ind+post_t(r_i);
            V_m_i = [V_m_i;V_m(n_i,start_i:end_i)]; %#ok<*AGROW>
            G_syn_E_E_i = [G_syn_E_E_i;G_syn_E_E(n_i,start_i:end_i)];
            G_syn_E_I_i = [G_syn_E_I_i;G_syn_E_I(n_i,start_i:end_i)];
            G_syn_I_E_i = [G_syn_I_E_i;G_syn_I_E(n_i,start_i:end_i)];
            G_syn_I_I_i = [G_syn_I_I_i;G_syn_I_I(n_i,start_i:end_i)];
        end
        neur_conductances(n_i).V_m = V_m_i;
        neur_conductances(n_i).G_syn_E_E = G_syn_E_E_i;
        neur_conductances(n_i).G_syn_E_I = G_syn_E_I_i;
        neur_conductances(n_i).G_syn_I_E = G_syn_I_E_i;
        neur_conductances(n_i).G_syn_I_I = G_syn_I_I_i;
        clear V_m_i G_syn_E_E_i G_syn_E_I_i G_syn_I_E_i G_syn_I_I_i r_i
    end    
    clear spikes_bin_cropped n_i spike_t pre_t post_t
    
    %Send the pulled data to the plotting function
    if parameters.plotResults
        plot_pre_spike_conductance_changes(parameters, neur_conductances)
    end    
    
end