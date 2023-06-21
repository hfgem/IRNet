%Code to test parameters

clear all

%% Get file directory and add functions to path

current_file = matlab.desktop.editor.getActiveFilename;
split_dir = split(current_file,'/');
joined_dir = join(split_dir(1:length(split_dir)-2),'/');
current_dir = strcat(joined_dir{1},'/');

cd(current_dir)

addpath('functions')
addpath('functions/plotting')

clear current_file split_dir joined_dir current_dir

%% Simulation options and save path

parameters.saveFlag = 1; % 1 to save simulation results
parameters.selectPath = 1; % 1 to select save destination, 0 to save in current dir
parameters.plotResults = 1; % 1 to plot basic simulation results

disp("Select Save Path for Results and Parameters")
save_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output stored

disp("Select Path with Base Parameters to Use")
param_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Save Folder'); %Have user input where they'd like the output store

%Upload a parameter file to use as a base (use IR_Net.m to make one if not):
load(strcat(param_path,'/parameters.mat'))
parameters.save_path = save_path;

%% Parameters to modify
%To set up a parameter for testing, use the following format:
%{
variedParam(1).name = 'W_gin'; % 1st parameter to be varied. Must be a field in the parameter structure
variedParam(1).range = linspace(450*10^(-12), 1050*10^(-12), test_n); % set of values to test param1
variedParam(2).name = 'del_G_syn_E_E'; % 2nd parameter to be varied
variedParam(2).range = linspace(450*10^(-12), 1050*10^(-12), test_n); % set of values to test param2
%}

%Number of values to test for each parameter
test_n = 10;
parameters.nNets = 1; %Number of network initializations
parameters.nTrials = 1; %Number of input initializations

%Varied parameter structure
variedParam = struct;
variedParam(1).name = 't_amp'; %Wave amplitude
variedParam(1).range = linspace(4*10^(-9),12*10^(-9),test_n);
variedParam(2).name = 'N_amp'; %Pink noise max amplitude
variedParam(2).range = linspace(0,1*10^(-9),test_n);

%Relevant parameters to above varied ones + simulation
parameters.dt = 1*10^(-3); %timestep (s)
parameters.init_period = 30; %initialization time (s)
parameters.sim_period = 270; %simulation time (s)
parameters.inputType = 2; % 0 = randn(), 1 = poisson, 2 = theta + pink noise
parameters.t_freq = 2; %SWR frequency - Nitzan et al. 2022

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);

% Save cases and updated parameters
if parameters.saveFlag
    save(strcat(save_path,'/parameters.mat'),'parameters')
    save(strcat(save_path,'/variedParam.mat'),'variedParam')
end

disp('Varied Parameters Set')

%% Run Grid Search with Burst Stats Returned

gcp; % starts parallel pool if not already running

% Waitbar code
D = parallel.pool.DataQueue;
h = waitbar(0, 'Starting simulation ...');
num_files = size(parameterSets_vec, 2);
nUpdateWaitbar(num_files, h); % Dummy call to nUpdateWaitbar to initialise
afterEach(D, @nUpdateWaitbar);

tic
netresults = struct;
testresults = struct;
parfor ithParamSet = 1:size(parameterSets_vec, 2) %Run through all parameter combinations
    %For each combination run parallelize_parameter_tests
    [netresults(ithParamSet),testresults(ithParamSet)] = parallelize_parameter_tests(...
                parameters, parameterSets_vec, ithParamSet, variedParam);
    send(D, 1);   
end
runTime = toc;
sprintf('Program Runtime (s) = %.2f',runTime)

%% Analyze Burst Stats

%Set parameter pairs for visualization
num_params = length(variedParam);
pairs = nchoosek(1:num_params,2);
rescale = size(netresults,2)/(test_n^2);
%Below lines need user inputs
[P1,P2] = ind2sub([num_params,test_n],1:size(parameterSets_vec,2)); %update left to equal number of params
param_inds = [P1;P2];

%Plot results
for p_i = 1:size(pairs,1)
    pair_ind_1 = pairs(p_i,1);
    pair_ind_2 = pairs(p_i,2);
    param_name_1 = replace(variedParam(pair_ind_1).name,'_',' ');
    param_name_2 = replace(variedParam(pair_ind_2).name,'_',' ');
    avg_neur_per_burst_results_mat = zeros(test_n,test_n); %param1 x param2
    avg_length_of_burst_results_mat = zeros(test_n,test_n); %param1 x param2
    avg_ibi_of_burst_results_mat = zeros(test_n,test_n); %param1 x param2
    for t_i = 1:size(netresults,2)
        avg_neur_per_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) = avg_neur_per_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) + (netresults(t_i).avg_neur_per_burst/rescale);
        avg_length_of_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) = avg_length_of_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) + (netresults(t_i).avg_length_of_burst/rescale);
        avg_ibi_of_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) = avg_ibi_of_burst_results_mat(param_inds(pair_ind_1,t_i),param_inds(pair_ind_2,t_i)) + (netresults(t_i).avg_ibi_of_bursts/rescale);
    end
    f = figure;
    %Avg Burst Size (# Neurons)
    s1 = subplot(2,2,1);
    imagesc(avg_neur_per_burst_results_mat)
    c1 = colorbar();
    c1.Label.String = 'Avg Number of Neurons';
    yticks(1:test_n)
    yticklabels([variedParam(pair_ind_1).range])
    ylabel(param_name_1)
    xticks(1:test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xlabel(param_name_2)
    title('Avg Num Neur per Burst')
    %Avg Burst Length
    s2 = subplot(2,2,2);
    imagesc(avg_length_of_burst_results_mat)
    c2 = colorbar();
    c2.Label.String = 'Avg Length of Burst (s)';
    yticks(1:test_n)
    yticklabels([variedParam(pair_ind_1).range])
    ylabel(param_name_1)
    xticks(1:test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xlabel(param_name_2)
    title('Avg Burst Length (s)')
    %Avg IBI
    s3 = subplot(2,2,3);
    imagesc(flipud(avg_ibi_of_burst_results_mat))
    c3 = colorbar();
    c3.Label.String = 'Avg IBI (s)';
    yticks(1:test_n)
    yticklabels([variedParam(pair_ind_1).range])
    ylabel(param_name_1)
    xticks(1:test_n)
    xticklabels([variedParam(pair_ind_2).range])
    xlabel(param_name_2)
    title('Avg IBI (s)')
    %Big title
    title_str = [param_name_1,' vs. ',param_name_2];
    sgtitle(title_str)
    %Link axes
    linkaxes([s1,s2,s3])
    %Save
    if parameters.saveFlag
       save_str = [variedParam(pair_ind_1).name,'_vs_',variedParam(pair_ind_2).name,'_avg_net_results'];
       savefig(f, strcat(save_path,'/',save_str,'.fig'))
       saveas(f, strcat(save_path,'/',save_str,'.svg'))
    end
end


%% Functions 
function p = nUpdateWaitbar(data, h)
% https://www.mathworks.com/matlabcentral/answers/660793-help-with-parfor-progress-bar-using-data-queue
    persistent TOTAL COUNT H
    if nargin == 2
        % initialisation mode
        H = h;
        TOTAL = data;
        COUNT = 0;
    else
        % afterEach call, increment COUNT
        COUNT = 1 + COUNT;
        p = COUNT / TOTAL;
        waitbar(p, H, ...
            ['Simulation ', num2str(COUNT), ' of ', num2str(TOTAL), ' complete', newline...
             'Runtime: ', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'), ...
            ', Remaining time: ', datestr(datenum(0,0,0,0,0,toc/p-toc),'HH:MM:SS')]);
    end
end



