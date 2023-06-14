function plot_rast(parameters,spikes_bin,ithTest,net_save_path)
    %_________
    %ABOUT: This function plots the full simulation raster
    %INPUTS:
    %   - parameters: parameter structure file
    %   - V_m: matrix with membrane potential of each neuron over time
    %   - ithTest: number of input current simulation
    %   - net_save_path: save path for network figure results
    %OUTPUTS:
    %   - saved or visualized figure - depending on save flag
    
    %Plot results and store
    init_end = ceil(parameters.init_period/parameters.dt);
    if parameters.saveFlag
        f = figure('visible','off');
        imagesc(spikes_bin(:,init_end:end))
        colormap(flipud(gray))
        colorbar()
        xtick_vals = xticks();
        xticklabels(xtick_vals*parameters.dt)
        xlabel('Time (s)')
        ylabel('Neuron Index')
        title('Full Simulation Raster Plot')
        savefig(f,strcat(net_save_path,'/sim_',string(ithTest),'_raster.fig'))
        saveas(f,strcat(net_save_path,'/sim_',string(ithTest),'_raster.svg'))
        close(f)
    end
    
end