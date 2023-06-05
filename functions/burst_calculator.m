function bursts = burst_calculator(parameters,V_m)
    %Spike raster data
    spike_raster = V_m > parameters.V_th;
    spike_sum = sum(spike_raster,1);
    spikes_exist = find(spike_sum);
    spike_diff = spikes_exist(2:end) - spikes_exist(1:end-1);
    large_diff = find(spike_diff*parameters.dt > parameters.burst_t_min);
    
    %Calculate and store burst data
    bursts = struct;
    start_ind = 1;
    for d_ind = 1:length(large_diff)
        bursts(end+1).times = [spikes_exist(start_ind),spikes_exist(large_diff(d_ind))]; %#ok<*AGROW>
        bursts(end).raster = spike_raster(:,bursts(end).times(1):bursts(end).times(2));
        bursts(end).num_neur = sum(sum(bursts(end).raster,2) > 0);
        bursts(end).length = parameters.dt*(bursts(end).times(2) - bursts(end).times(1));
        try
            bursts(end).ibi = bursts(end).times(1) - bursts(end-1).times(2); %Time since last burst
        catch %If first burst saved
            bursts(end).ibi = bursts(end).times(1);
        end
        if bursts(end).times(1) > parameters.I_E_delay/parameters.dt
            bursts(end).pre_delay = 1;
        else
            bursts(end).pre_delay = 0;
        end
        %Reset the start ind for the next burst
        start_ind = large_diff(d_ind) + 1;
    end    
    bursts(1) = [];

end