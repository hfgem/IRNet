function bursts = burst_calculator(parameters,spike_raster,net_save_path,ithNet,ithTest)
    
    %Calculate where to start looking for bursts
    init_end = ceil(parameters.init_period/parameters.dt);
    
    %Spike raster data
    spike_sum = sum(spike_raster(:,init_end:end),1);
    spikes_exist = find(spike_sum) + init_end - 1;
    spike_diff = spikes_exist(2:end) - spikes_exist(1:end-1);
    large_diff = find(spike_diff*parameters.dt > parameters.burst_t_min);
    
    %Calculate and store burst statistics
    bursts = struct;
    start_ind = 1;
    for d_ind = 1:length(large_diff)
        times = [spikes_exist(start_ind),spikes_exist(large_diff(d_ind))];
        raster = spike_raster(:,times(1):times(2));
        num_neur = sum(sum(raster,2) > 0);
        if num_neur > parameters.burst_n_min*parameters.n
            bursts(end+1).times = times; %#ok<*AGROW>
            bursts(end).raster = raster;
            bursts(end).num_neur = num_neur;
            bursts(end).length = parameters.dt*(bursts(end).times(2) - bursts(end).times(1));
            try
                bursts(end).ibi = (bursts(end).times(1) - bursts(end-1).times(2))*parameters.dt; %Time since last burst in seconds
            catch %If first burst saved
                bursts(end).ibi = bursts(end).times(1)*parameters.dt;
            end
        end
        %Reset the start ind for the next burst
        start_ind = large_diff(d_ind) + 1;
    end    
    bursts(1) = [];
    
    %Plot statistics and individual bursts if plot flag set
    if parameters.plotResults == 1
        if length(fieldnames(bursts)) > 0
            plot_burst_stats(parameters,bursts,ithNet,ithTest,net_save_path)
            plot_bursts(parameters,bursts,ithTest,net_save_path)
        end
    end
    
end