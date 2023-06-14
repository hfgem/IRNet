function plot_conns(parameters,conns,ithNet,ithTest,net_save_path)
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

    f = figure('visible','off');
    ax1 = subplot(1,2,1);
    imagesc(squeeze(conns(:,:,1)))
    title('Starting Connectivity')
    colorbar()
    ax2 = subplot(1,2,2);
    imagesc(squeeze(conns(:,:,end)))
    title('Ending Connectivity')
    colorbar()
    linkaxes([ax1,ax2])
    sgtitle(strcat('Network ',string(ithNet),', Simulation ',string(ithTest)))
    if parameters.saveFlag
        savefig(f,strcat(net_save_path,'/conn_changes_',string(ithTest),'.fig'))
        saveas(f,strcat(net_save_path,'/conn_changes_',string(ithTest),'.svg'))
        close(f)
    end
end