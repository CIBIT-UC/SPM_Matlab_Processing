clear, clc, close all

%% Load PhysIO paths
tool_path = 'C:\Users\alexandresayal\Documents\GitHub\tapas\PhysIO\code';

addpath(tool_path);
tapas_physio_init();

%% Load and settings
load('Configs_VP_INHIBITION_SPM.mat')

% ============================= Settings ============================== %
subjectName = 'VPIS20';  % Subject name
fwhm = 5;                      % SS FWHM in mm
transf = 'MNI';                % 'TAL' OR 'MNI'
%=======================================================================%

s = find(not(cellfun('isempty', strfind(datasetConfigs.subjects, subjectName))));

% Be subject specific and remove the anatomical run
datasetConfigs.runs = datasetConfigs.runs{s}(2:end);
datasetConfigs.volumes = datasetConfigs.volumes{s}(2:end);
% datasetConfigs.tr = datasetConfigs.TR{s};
datasetConfigs.tr = 1;
% datasetConfigs.NrOfSlices = datasetConfigs.NrOfSlices{s};
datasetConfigs.NrOfSlices = 42;

%% Provide some information :)
data_configs = struct();
% data_configs.nslices = 42;
% data_configs.tr = 1; % in seconds
data_configs.ndummies = 0;
% data_configs.nscans = 374; % give later, this is cheating ¯\_(?)_/¯
data_configs.onsetslice = 1; % not sure yet what this does

%% Iterate on the runs and generate regressors
nRuns = length(datasetConfigs.runs); % Provide this number manually because reasons

for rr = 1:nRuns
    
    % Provide run-specific data
    data_configs.nslices = datasetConfigs.NrOfSlices;
    data_configs.tr = datasetConfigs.tr;
    data_configs.nscans = datasetConfigs.volumes(rr);
    
    % Create physio regressors
    file_path = fullfile(datasetConfigs.path,subjectName,...
        ['r-' datasetConfigs.runs{rr}],'PHYSIO');
    file_prefix = datasetConfigs.runs{rr};
    
    try
        [Regressors] = generatePhysIORegressors(file_path,file_prefix,data_configs);
        
        % Export txt
        dlmwrite(fullfile(file_path,'physio_regressors.txt'),Regressors,'delimiter','\t','precision','%.10f');
        
        disp('next')
    catch
        fprintf('!!!! FAILED on run %s!!!!!!\n',datasetConfigs.runs{rr})
    end
end
