function plot_pre_spike_conductance_changes(parameters, network, ...
    neur_conductances)
    %_________
    %ABOUT: This function plots input conductances a set amount of time 
    %prior to a spike. The plot has different types of conductances
    %separated, as well as the average waveform of conductance inputs prior 
    %to a neuron spiking.
    %
    %INPUTS:
    %   parameters: structure containing simulation parameters
    %   network: 
    %   neur_conductances: structure containing 
    %
    %OUTPUTS: A figure with plotted conductances for different connections
    %and an average conductance profile. The figure has 6 subplots
    %arranged in a 3x2 orientation. Going left-to-right and top-to-bottom:
    %   - subplot 1 contains the average membrane potential waveform around
    %       spike time (across all neurons)
    %   - subplot 2 contains the average E-E conductance waveform around 
    %       spike time
    %   - subplot 3 contains the average E-I conductance waveform around
    %       spike time
    %   - subplot 4 contains the average I-E conductance waveform around
    %       spike time
    %   - subplot 5 contains the average I-I conductance waveform around
    %       spike time
    %   - subplot 6 contains the average conductance waveform around spike
    %       time, averaging all types of conductances
    %_________
    
    %Variables
    pre_buffer = ceil(parameters.pre_spike/parameters.dt);
    post_buffer = ceil(parameters.post_spike/parameters.dt);
    
    %Compile all values
    V_m_all = [];
    G_syn_E_E_all = [];
    G_syn_E_I_all = [];
    G_syn_I_E_all = [];
    G_syn_I_I_all = [];
    for n_i = 1:parameters.n
        V_m_all = [V_m_all; [neur_conductances(n_i).V_m]];
        I_ind = find(network.I_indices == n_i, 1);
        if isempty(I_ind)
            G_syn_E_E_all = [G_syn_E_E_all; [neur_conductances(n_i).G_syn_E_E]];
            G_syn_I_E_all = [G_syn_I_E_all; [neur_conductances(n_i).G_syn_I_E]];
        else
            G_syn_E_I_all = [G_syn_E_I_all; [neur_conductances(n_i).G_syn_E_I]];
            G_syn_I_I_all = [G_syn_I_I_all; [neur_conductances(n_i).G_syn_I_I]];
        end    
    end
    clear n_i
    %Plot
    f = figure;
    %Average V_m
    ax1 = subplot(2,3,1);
    plot(-pre_buffer:post_buffer,mean(V_m_all,1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Membrane Potential (V)')
    title('Average Membrane Potential')
    %Average G_E_E
    ax2 = subplot(2,3,2);
    plot(-pre_buffer:post_buffer,mean(G_syn_E_E_all,1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Conductance E-E (S)')
    title('Average E-E Conductance')
    %Average G_E_I
    ax3 = subplot(2,3,3);
    plot(-pre_buffer:post_buffer,mean(G_syn_E_I_all,1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Conductance E-I (S)')
    title('Average E-I Conductance')
    %Average G_I_E
    ax4 = subplot(2,3,4);
    plot(-pre_buffer:post_buffer,mean(G_syn_I_E_all,1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Conductance I-E (S)')
    title('Average I-E Conductance')
    %Average G_I_I
    ax5 = subplot(2,3,5);
    plot(-pre_buffer:post_buffer,mean(G_syn_I_I_all,1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Conductance I-I (S)')
    title('Average I-I Conductance')
    %Average G
    ax6 = subplot(2,3,6);
    plot(-pre_buffer:post_buffer,mean([G_syn_E_E_all;G_syn_E_I_all;G_syn_I_E_all;G_syn_I_I_all],1))
    x_tick_vals = xticks();
    xticklabels(x_tick_vals*parameters.dt)
    xlabel('Time (s)')
    ylabel('Conductance E-I (S)')
    title('Average Conductance')
    %link axes
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
    %save figure
    if parameters.saveFlag
        savefig(f,strcat(parameters.save_path,'/conductance_waveforms.fig'))
        savefig(f,strcat(parameters.save_path,'/conductance_waveforms.svg'))
        close(f)
    end    
end