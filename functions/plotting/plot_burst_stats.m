function plot_burst_stats(parameters,bursts,ithNet,ithTest,net_save_path)
    %_________
    %ABOUT: This function plots the connection strength changes
    %INPUTS:
    %   - parameters: parameter structure file
    %   - conns: matrix of [n x n x number of saved connection matrices]
    %   - ithNet: number of network structure simulation
    %   - ithTest: number of input current simulation
    %   - net_save_path: save path for network figure results
    %OUTPUTS:
    %   - saved or visualized figure - depending on save flag

    %Calculate average burst length as a function of time
    num_chunks = ceil(length(bursts)/parameters.num_burst_avg);
    burst_avg_len = zeros(1,num_chunks);
    burst_avg_sz = zeros(1,num_chunks);
    burst_avg_ibi = zeros(1,num_chunks);
    for n_i = 1:num_chunks
        burst_ind = max([0,(n_i-1)*parameters.num_burst_avg + 1]);
        burst_end_ind = min([length(bursts),burst_ind + parameters.num_burst_avg]);
        burst_avg_len(n_i) = mean([bursts(burst_ind:burst_end_ind).length]);
        burst_avg_sz(n_i) = mean([bursts(burst_ind:burst_end_ind).num_neur]);
        burst_avg_ibi(n_i) = mean([bursts(burst_ind:burst_end_ind).ibi]);
    end
    
    %Plot results and store if desired
    if parameters.saveFlag
        f = figure('visible','off');
    else
        f = figure;
    end
    subplot(3,2,1)
    plot([bursts.length])
    ylabel('Time (s)')
    title('Burst Lengths')
    subplot(3,2,3)
    plot([bursts.num_neur])
    ylabel('Num Neur')
    title('Burst Size (Number of Neurons)')
    subplot(3,2,5)
    plot([bursts.ibi])
    y_tick_vals = yticks();
    yticklabels(y_tick_vals)
    ylabel('Time (s)')
    title('Time to Next Burst')
    subplot(3,2,2)
    plot(1:num_chunks,burst_avg_len)
    ylabel('Time (s)')
    title(['Average Burst Length ', 'for ' + string(parameters.num_burst_avg) + ' Bursts'])
    subplot(3,2,4)
    plot(1:num_chunks,burst_avg_sz)
    ylabel('Num Neur')
    title(['Average Burst Size ', 'for ' + string(parameters.num_burst_avg) + ' Bursts'])
    subplot(3,2,6)
    plot(1:num_chunks,burst_avg_ibi)
    ylabel('Avg Time (s)')
    title(['Average Burst IBI ','for ' + string(parameters.num_burst_avg) + ' Bursts'])
    sgtitle(strcat('Network ',string(ithNet),', Simulation ',string(ithTest)))
    if parameters.saveFlag
        savefig(f,strcat(net_save_path,'/burst_results_',string(ithTest),'.fig'))
        saveas(f,strcat(net_save_path,'/burst_results_',string(ithTest),'.svg'))
        close(f)
    end
    
end