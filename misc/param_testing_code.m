%Code to test parameters

clear all
close all

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
parameters.test_n = test_n;
parameters.nNets = 3; %Number of network initializations
parameters.nTrials = 1; %Number of input initializations

%Varied parameter structure
variedParam = struct;
variedParam(1).name = 'N_amp'; %Pink noise max amplitude
variedParam(1).range = round(linspace(0,1*10^(-9),parameters.test_n),3,'significant'); %added rounding
variedParam(2).name = 'tau_I_decay'; %Decay Timescale for Inhibitory Connections
variedParam(2).range = round(linspace(-180,-30,parameters.test_n)./log(0.5),3,'significant'); %added rounding
variedParam(3).name = 'tau_E_decay'; %Decay Timescale for Excitatory Connections
variedParam(3).range = round(linspace(-180,-30,parameters.test_n)./log(0.5),3,'significant'); %added rounding

%Relevant parameters to above varied ones + simulation
parameters.dt = 1*10^(-3); %timestep (s)
parameters.init_period = 30; %initialization time (s)
parameters.sim_period = 270; %simulation time (s)
parameters.inputType = 0; % 0 = randn(), 1 = poisson, 2 = theta + pink noise

%"Good Bursts" Settings
parameters.min_burst_len = 5*10^(-3); %Minimum "good" burst length (s)
parameters.max_burst_len = 500*10^(-3); %Maximum "good" burst length (s)
parameters.min_ibi = 50*10^(-3); %Minimum inter-burst-interval
parameters.min_n_burst = parameters.n*0.05; %Minimum number of neurons in a burst

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
netresults = struct('avg_neur_per_burst',[],'avg_length_of_burst',...
    [],'avg_ibi_of_bursts',[]);
testresults = struct;
%Then parallelize the rest
parfor ithParamSet = 1:size(parameterSets_vec, 2) %Run through all parameter combinations
    %For each combination run parallelize_parameter_tests
    [netresults(ithParamSet), testresults(ithParamSet).results] = parallelize_parameter_tests(...
                parameters, parameterSets_vec, ithParamSet, variedParam);
    send(D, 1);   
end
runTime = toc;
sprintf('Program Runtime (s) = %.2f',runTime)
disp('Now Saving Results')
save(strcat(parameters.save_path,'/netresults.mat'),'netresults','-v7.3'); 
save(strcat(parameters.save_path,'/testresults.mat'),'testresults','-v7.3'); 

%% Determine "Good Parameters"

num_param_sets = size(parameterSets_vec,2);
good_params = [];
for p_i = 1:num_param_sets
    %avg # neurons per burst
    if netresults(p_i).avg_neur_per_burst >= parameters.min_n_burst
        if (parameters.min_burst_len <= netresults(p_i).avg_length_of_burst) && (netresults(p_i).avg_length_of_burst <= parameters.max_burst_len)
            if parameters.min_ibi <= netresults(p_i).avg_ibi_of_bursts
                good_params = [good_params; parameterSets_vec(:,p_i)']; %#ok<SAGROW>
            end    
        end
    end    
    
end    
save(strcat(parameters.save_path,'/good_params.mat'),'good_params','-v7.3'); 

%Display/save which ranges are good
range_file = strcat(parameters.save_path,'/ranges.txt');
varnames = {'Parameter','Minimum','Maximum'};
vartypes = {'string','double','double'};
sz = [size(variedParam,2),3];
range_table = table('Size',sz,'VariableTypes',vartypes,'VariableNames',varnames);
for v_i = 1:size(variedParam,2)
    min_val = min(good_params(:,v_i));
    max_val = max(good_params(:,v_i));
    param_name = replace(variedParam(v_i).name,'_',' ');
    range_table(v_i,:) = {param_name,min_val,max_val};
    string_range = strcat(param_name,' good range: [',...
        string(min_val),',',string(max_val),'].');
    disp(string_range)
end
writetable(range_table,range_file)
%% Plot Burst Stats

plot_param_test_results(parameters,variedParam,parameterSets_vec,netresults)

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



