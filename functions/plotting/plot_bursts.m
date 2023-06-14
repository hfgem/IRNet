function plot_bursts(parameters,bursts,ithTest,net_save_path)
    %_________
    %ABOUT: This function plots the burst rasters
    %INPUTS:
    %   - parameters: parameter structure file
    %   - bursts: structure with burst info
    %   - ithTest: number of input current simulation
    %   - net_save_path: save path for network figure results
    %OUTPUTS:
    %   - saved or visualized figure - depending on save flag

    %Make save directory
    burst_dir = net_save_path + '/burst_rasters';
    mkdir(burst_dir)
    burst_test_dir = burst_dir + '/sim_' + string(ithTest);
    mkdir(burst_test_dir)
    
    %Plot results and store
    if parameters.saveFlag
        plot_ind = randperm(length(bursts),min(parameters.num_rast_to_plot,length(bursts)));
        for i = 1:length(plot_ind)
            b_i = plot_ind(i);
            f = figure('visible','off');
            imagesc(bursts(b_i).raster)
            colormap(flipud(gray))
            colorbar()
            xtick_vals = xticks();
            xticklabels(xtick_vals*parameters.dt + bursts(b_i).times(1)*parameters.dt)
            xlabel('Time (s)')
            ylabel('Neuron Index')
            title(strcat('Burst #',string(i),' Raster Plot'))
            savefig(f,strcat(burst_test_dir,'/burst_',string(b_i),'_raster.fig'))
            saveas(f,strcat(burst_test_dir,'/burst_',string(b_i),'_raster.svg'))
            close(f)
        end
    end
    
end