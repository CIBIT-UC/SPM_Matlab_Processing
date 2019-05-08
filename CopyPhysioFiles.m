clear

load('Configs_VP_INHIBITION_SPM.mat')

rawdatapath = fullfile('F:','RAW_DATA_VP_INHIBITION');

subjectName = 'VPIS20';
rawdatapath = fullfile(rawdatapath,[subjectName '_LOGS']);

s = find(not(cellfun('isempty', strfind(datasetConfigs.subjects, subjectName))));

runs = datasetConfigs.runs{s};
runs = runs(2:end);

for rr = 1:length(runs)
    
    %--- Ask for files
    [file,path] = uigetfile('*.log',...
        sprintf('Select _Info, _RESP, _PULS Files for Run %i',rr), rawdatapath , ...
        'MultiSelect', 'on');
    
    physiofolder = fullfile(datasetConfigs.path,subjectName,['r-' runs{rr}],'PHYSIO');
    
    if ~exist(physiofolder,'dir'); mkdir(physiofolder); end
    
    for ff = 1:length(file)
        
        aux = strsplit(file{ff},'_');
        new_file = [aux{1} '_' runs{rr} '_' aux{end}];
        
        copyfile(fullfile(path,file{ff}),...
            fullfile(physiofolder,new_file));
        
    end
    
end

disp('Done')
