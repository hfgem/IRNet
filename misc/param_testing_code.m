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
test_n = 2;
parameters.nNets = 1; %Number of network initializations
parameters.nTrials = 1; %Number of input initializations

%Varied parameter structure
variedParam = struct;
variedParam(1).name = 'G_std';
variedParam(1).range = linspace(8*10^(-9),12*10^(-9),test_n);
variedParam(2).name = 'G_mean';
variedParam(2).range = linspace(0,8*10^(-9),test_n);

%Above parameters are for Gaussian input, so:
parameters.inputType = 0; % 0 = randn(), 1 = poisson, 2 = theta + randn()

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
resultsStructLinear = cell(1, size(parameterSets_vec, 2));
parfor ithParamSet = 1:size(parameterSets_vec, 2) %Run through all parameter combinations
    %For each combination run parallelize_parameter_tests
    [resultsStructLinear{ithParamSet}] = parallelize_parameter_tests(...
                parameters, parameterSets_vec, ithParamSet, variedParam);
    send(D, 1);   
end
runTime = toc;
sprintf('Program Runtime (s) = %.2f',runTime)

%% Analyze Burst Stats


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



