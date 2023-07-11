function [net_burst_results, test_burst_var] = parallelize_parameter_tests(parameters,...
    parameterSets_vec, ithParamSet, variedParam)
    %_____
    %ABOUT: This function runs through a series of commands to test
    %different parameters in IR_Net, with the goal of finding good
    %parameter ranges for bursting.
    %INPUTS:
    %
    %OUTPUTS:
    %
    %_____
    
    % Set up parameter values for current parameter set
    for i = 1:size(variedParam, 2)
        parameters.(variedParam(i).name) = parameterSets_vec(i,ithParamSet);
    end
    
    % Update any parameters that are dependent on a varied parameter
    parameters = set_dependent_parameters(parameters);
    
    param_save_path = strcat(parameters.save_path,'/',string(ithParamSet));
    if parameters.saveFlag
        if ~isfolder(param_save_path)
            mkdir(param_save_path);
        end
        save(strcat(param_save_path,'/parameters.mat'),'parameters','-v7.3'); 
    end
    
    %Set up storage variables
    net_burst_results = struct;
    test_burst_results = struct;
    
    %Run simulation
    for ithNet = 1:parameters.nNets
        
        %CREATE NETWORK SAVE PATH
        net_save_path = strcat(param_save_path,'/network_',string(ithNet));
        if parameters.saveFlag && ~isfolder(net_save_path)
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
        test_burst_var = struct;
        if parameters.run_light ~= 1
            V_m_var = struct;
            G_var = struct;
            I_var = struct;
            %conns_var = struct;
        end
        for ithTest = 1:parameters.nTrials       
            %Create input conductance variable
            if parameters.inputType == 0
                % Pink noise input
                try %With installed audio package
                    pink_noise = pinknoise(parameters.t_steps+1,parameters.n);
                    pink_noise_scaled = parameters.N_amp*(pink_noise./max(pink_noise));
                catch %For older MATLAB versions or uninstalled audio package
                    pink_noise = dsp.ColoredNoise('Color','pink','NumChannels',parameters.n,'SamplesPerFrame',parameters.t_steps+1);
                    pink_noise_scaled = parameters.N_amp*(pink_noise()./max(pink_noise()));
                    pink_noise_scaled = pink_noise_scaled';
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
                    pink_noise_scaled = pink_noise_scaled';
                end
                theta_wave = (parameters.t_amp/2)*sin(2*pi*parameters.t_freq*(parameters.dt*ones(parameters.n,parameters.t_steps+1).*(1:parameters.t_steps+1))) + parameters.t_amp/2;
                G_in = theta_wave + pink_noise_scaled;
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
                conns(:,:,2) = last_conns; %#ok<*NASGU>
                clear last_conns
            else
                %Create Storage Variables
                V_m = zeros(parameters.n,parameters.t_steps+1); %membrane potential for each neuron at each timestep
                V_m(:,1) = parameters.V_reset + randn([parameters.n,1])*(10^(-3))*sqrt(parameters.dt); %set all neurons to baseline reset membrane potential with added noise
                [V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I, conns] = ir_net_calculator(...
                    parameters, seed, network, V_m); %#ok<*ASGLU>
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
            
            %Calculate bursts
            bursts = burst_calculator(parameters,spikes_bin,net_save_path,ithNet,ithTest);
            test_index = parameters.nTrials*(ithNet-1) + ithTest;
            test_burst_var(test_index).bursts = bursts;
            clear bursts
            
        end %end test loop

        %SAVE NETWORK DATA
        if parameters.saveFlag
            save(strcat(net_save_path,'/burst_var.mat'),'test_burst_var','-v7.3')
            if parameters.run_light ~= 1
                save(strcat(net_save_path,'/V_m_var.mat'),'V_m_var','-v7.3')
                save(strcat(net_save_path,'/G_var.mat'),'G_var','-v7.3')
                save(strcat(net_save_path,'/I_var.mat'),'I_var','-v7.3')
                %save(strcat(net_save_path,'/conns_var.mat'),'conns_var','-v7.3')
            end
        end
        
        %Calculate average burst data for the network and store
        avg_neur_per_burst = 0;
        avg_length_of_burst = 0;
        avg_ibi_of_bursts = 0;
        for ithTest = 1:parameters.nTrials
            num_bursts = length(test_burst_var(ithTest).bursts);
            for b_i = 1:num_bursts
                avg_neur_per_burst = avg_neur_per_burst + (test_burst_var(ithTest).bursts(b_i).num_neur / num_bursts);
                avg_length_of_burst = avg_length_of_burst + (test_burst_var(ithTest).bursts(b_i).length / num_bursts);
                avg_ibi_of_bursts = avg_ibi_of_bursts + (test_burst_var(ithTest).bursts(b_i).ibi / num_bursts);
            end
        end
        
        %Save average burst data for the network
        net_burst_results(ithNet).avg_neur_per_burst = avg_neur_per_burst;
        net_burst_results(ithNet).avg_length_of_burst = avg_length_of_burst;
        net_burst_results(ithNet).avg_ibi_of_bursts = avg_ibi_of_bursts;
        
    end 
    
    if length(net_burst_results) > 1
        avg_avg_neur_per_burst = mean([net_burst_results.avg_neur_per_burst]);
        avg_avg_length_of_burst = mean([net_burst_results.avg_length_of_burst]);
        avg_avg_ibi_of_bursts = mean([net_burst_results.avg_ibi_of_bursts]);
        net_burst_results = struct;
        net_burst_results(1).avg_neur_per_burst = avg_avg_neur_per_burst;
        net_burst_results(1).avg_length_of_burst = avg_avg_length_of_burst;
        net_burst_results(1).avg_ibi_of_bursts = avg_avg_ibi_of_bursts;
    end    
    
    disp(['Parameter set ', num2str(ithParamSet), '/', num2str(size(parameterSets_vec, 2)), ' complete'])

end