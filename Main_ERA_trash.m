clear,clc
load('Configs_VP_INHIBITION_SPM.mat')
load(fullfile(datasetConfigs.AnalysisPath,'TCs-MT-5mm.mat'))

runs = datasetConfigs.runs{1,1}(2:end);
delay = 3; %in volumes

baselineIndexes = cell(1,nRuns);

% Load protocol and baseline indexes
for rr = 1:nRuns
   prtmatFile = fullfile(datasetConfigs.path,datasetConfigs.subjects{1},['r-' runs{rr}],'PRT',['RunMRI_D12_R' num2str(rr) '.mat']);
   
   load(prtmatFile)
   
   
    
end
