function plot_param_test_results(parameters,variedParam,parameterSets_vec,netresults)
     %_________
    %ABOUT: This function plots the results of parameter testing. Pairs of
    %varied parameters are plotted against each other using average values
    %and calculations of good results.
    %   - parameters: parameter structure file
    %   - variedParam: structure of varied parameter names and ranges
    %   - parameterSets_vec: variedParam converted into a matrix of
    %           parameter combinations
    %   - netresults: structure containing average outcomes of network
    %           simulations
    %OUTPUTS:
    %   - saved or visualized figure - depending on save flag

    %Set parameter pairs for visualization
num_params = length(variedParam);
pairs = nchoosek(1:num_params,2);
rescale = size(netresults,2)/(parameters.test_n^length(variedParam));

%Plot results
for p_i = 1:size(pairs,1)
    %Grab indices of 2 varied parameters to compare
    pair_ind_1 = pairs(p_i,1);
    pair_ind_2 = pairs(p_i,2);
    %Grab names of parameters to compare
    param_name_1 = replace(variedParam(pair_ind_1).name,'_',' ');
    param_name_2 = replace(variedParam(pair_ind_2).name,'_',' ');
    %Set up storage of average results
    avg_neur_per_burst_results_mat = zeros(parameters.test_n,parameters.test_n); %param1 x param2
    avg_length_of_burst_results_mat = zeros(parameters.test_n,parameters.test_n); %param1 x param2
    avg_ibi_of_burst_results_mat = zeros(parameters.test_n,parameters.test_n); %param1 x param2
    %Run through all parameter simulation results and save to matrices
    for t_i = 1:size(netresults,2)
        %Get the accurate parameter value index
        param_1_ind = find(variedParam(pair_ind_1).range == parameterSets_vec(pair_ind_1,t_i));
        param_2_ind = find(variedParam(pair_ind_2).range == parameterSets_vec(pair_ind_2,t_i));
        avg_neur_per_burst_results_mat(param_1_ind,param_2_ind) = avg_neur_per_burst_results_mat(param_1_ind,param_2_ind) + (netresults(t_i).avg_neur_per_burst/rescale);
        avg_length_of_burst_results_mat(param_1_ind,param_2_ind) = avg_length_of_burst_results_mat(param_1_ind,param_2_ind) + (netresults(t_i).avg_length_of_burst/rescale);
        avg_ibi_of_burst_results_mat(param_1_ind,param_2_ind) = avg_ibi_of_burst_results_mat(param_1_ind,param_2_ind) + (netresults(t_i).avg_ibi_of_bursts/rescale);
    end
    %Plot
    if parameters.saveFlag
        f = figure('Position', [10 10 800 600],'visible','off');
    else
        f = figure('Position', [10 10 800 600]);
    end
    %Avg Burst Size (# Neurons)
    s1 = subplot(3,2,1);
    imagesc(flipud(avg_neur_per_burst_results_mat))
    c1 = colorbar();
    caxis([0,parameters.n])
    c1.Label.String = 'Avg Number of Neurons';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Avg Num Neur per Burst')
    %Good Params for Neuron Count
    s2 = subplot(3,2,2);
    bin_good_neur_count = avg_neur_per_burst_results_mat >= parameters.min_n_burst;
    imagesc(flipud(bin_good_neur_count))
    c2 = colorbar('Ticks',[0,1],'TickLabels',{'Bad','Good'});
    caxis([0,1])
    c2.Label.String = 'Good Number of Neurons';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Good Parameters for Number of Neurons / Burst')
    %Avg Burst Length
    s3 = subplot(3,2,3);
    imagesc(flipud(avg_length_of_burst_results_mat))
    c3 = colorbar();
    caxis([0,max(avg_length_of_burst_results_mat,[],'all')])
    c3.Label.String = 'Avg Length of Burst (s)';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Avg Burst Length (s)')
    %Good Params for Burst Length
    s4 = subplot(3,2,4);
    bin_good_burst_len = (avg_length_of_burst_results_mat >= parameters.min_burst_len) .* (avg_length_of_burst_results_mat <= parameters.max_burst_len);
    imagesc(flipud(bin_good_burst_len))
    c4 = colorbar('Ticks',[0,1],'TickLabels',{'Bad','Good'});
    caxis([0,1])
    c4.Label.String = 'Good Length of Burst';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Good Parameters for Burst Length')
    %Avg IBI
    s5 = subplot(3,2,5);
    imagesc(flipud(avg_ibi_of_burst_results_mat))
    c5 = colorbar();
    caxis([0,max(avg_ibi_of_burst_results_mat,[],'all')])
    c5.Label.String = 'Avg IBI (s)';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Avg IBI (s)')
    %Good Params for IBI
    s6 = subplot(3,2,6);
    bin_good_ibi = avg_ibi_of_burst_results_mat >= parameters.min_ibi;
    imagesc(flipud(bin_good_ibi))
    c6 = colorbar('Ticks',[0,1],'TickLabels',{'Bad','Good'});
    caxis([0,1])
    c6.Label.String = 'Good Parameters for IBI';
    yticks(1:parameters.test_n)
    yticklabels(fliplr([variedParam(pair_ind_1).range]))
    ylabel(param_name_1)
    xticks(1:parameters.test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xtickangle(45)
    xlabel(param_name_2)
    title('Good Parameters for IBI')
    %Big title
    title_str = [param_name_1,' vs. ',param_name_2];
    sgtitle(title_str)
    %Link axes
    linkaxes([s1,s2,s3,s4,s5,s6])
    %Save
    if parameters.saveFlag
       save_str = [variedParam(pair_ind_1).name,'_vs_',variedParam(pair_ind_2).name,'_avg_net_results'];
       savefig(f, strcat(parameters.save_path,'/',save_str,'.fig'))
       saveas(f, strcat(parameters.save_path,'/',save_str,'.svg'))
    end
end

end