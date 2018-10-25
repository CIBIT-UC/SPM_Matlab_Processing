function [ success , DCMinfo ] = createFolderStructure( datasetConfigs , dataPath , dataTBV , subjectIndex )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

success = false;
DCMinfo = struct();
idx_info = 1;

% -------------------------------------------------------------------------
% Find DICOM files in dataPath
% -------------------------------------------------------------------------

D = dir(fullfile(dataPath,'*.ima'));

% Check if the total number of files is incorrect
if (length(D) < sum(datasetConfigs.volumes)) && (length(D) > 5)
    disp('[createFolderStructure] Fewer files than expected...');
% Maybe the files are *.dcm!
elseif (length(D) < sum(datasetConfigs.volumes)) && (length(D) < 5)
    D = dir(fullfile(dataPath,'*.dcm'));
end

% Extract all file names
files = extractfield(D,'name')';
nFiles = length(files);

% -------------------------------------------------------------------------
% Check if names on files match 
% -------------------------------------------------------------------------

% Retrieve first DICOM file
firstDCMfile = dicominfo(fullfile(D(1).folder,D(1).name));

% Check if the patient given name exists. If it does, the DICOM files were
% not anonimized.
if isfield(firstDCMfile.PatientName,'GivenName')
    % Compare the name given in datasetConfigs with the GivenName+Surname
    % in the DICOM file.
    if ~strcmpi(datasetConfigs.subjects{subjectIndex},...
            [firstDCMfile.PatientName.GivenName firstDCMfile.PatientName.FamilyName])
        
        disp('[createFolderStructure] Check if files correspond to subject!')
        fprintf('Name on DCM files: %s %s\n',...
            firstDCMfile.PatientName.GivenName,firstDCMfile.PatientName.FamilyName);
        fprintf('Name provided: %s\n',datasetConfigs.subjects{subjectIndex});
        
        x = input('[createFolderStructure] Do you wish to proceed anyway (Y/N)?','s');
        switch lower(x)
            case 'y'
                disp('[createFolderStructure] Proceeding...');
            otherwise
                return
        end
    end
else % the GivenName field does not exist
    % When using an anonimization standard, the name in datasetConfigs
    % should match the PatientName.FamilyName field in the DICOM file
    if ~strcmpi(datasetConfigs.subjects{subjectIndex},firstDCMfile.PatientName.FamilyName)
        disp('[createFolderStructure] Check if files correspond to subject!')
        fprintf('Name on DCM files: %s \n',firstDCMfile.PatientName.FamilyName);
        fprintf('Name provided: %s\n',datasetConfigs.subjects{subjectIndex});
        
        x = input('[createFolderStructure] Do you wish to proceed anyway (Y/N)?','s');
        switch lower(x)
            case 'y'
                disp('[createFolderStructure] Proceeding...');
            otherwise
                return
        end
    end   
end

% -------------------------------------------------------------------------
% Extract series
% -------------------------------------------------------------------------

series = zeros(nFiles,1);

% Iterate on the files and search for the series number
% The filenames are formated as <subjectID>.<MR>.<series>.(...)
for ii = 1:nFiles   
    auxnamesplit = strsplit(files{ii},'.');    
    series(ii) = str2double(auxnamesplit{3});  
end

% Find the unique series numbers
[seriesNumbers , seriesIdx] = unique(series);

% Confirm the series numbers with the information in the first DICOM file
% of each series/run
READ_HEADERS = false;

for ii = 1:length(seriesNumbers)
    file_idx = seriesIdx(ii);
    dcmInfo = dicominfo(fullfile(D(file_idx).folder,D(file_idx).name));
    if dcmInfo.SeriesNumber ~= seriesNumbers(ii)
        fprintf('[createFolderStructure] Series numbers do not match between filename and DICOM header: M1=%i M2=%i \n',...
            seriesNumbers(ii),dcmInfo.SeriesNumber);
        
        x = input('[createFolderStructure] Do you wish to read all DICOM headers (a Parallel pool will run this) (Y/N)?','s');
        switch lower(x)
            case 'y'
                READ_HEADERS = true;
                break
            otherwise
                return
        end
    end      
end

% If the user chooses to read all DICOM headers, this block will run
if READ_HEADERS
    series = zeros(nFiles,1);
    parfor ii = 1:nFiles 
        dcmInfo = dicominfo(fullfile(D(ii).folder,D(ii).name));
        series(ii) = dcmInfo.SeriesNumber;
    end
    
    % Find the unique series numbers
    seriesNumbers = unique(series);
end

% Retrieve number of files per series
seriesVolumes = hist(series,length(1:seriesNumbers(end)));
seriesVolumes = seriesVolumes(seriesVolumes~=0);

% -------------------------------------------------------------------------
% Check for incomplete runs or extra runs
% -------------------------------------------------------------------------

nRuns = length(datasetConfigs.runs);

% Number of series larger than expected number of runs
if length(seriesNumbers) > nRuns
    
    % Find series with strange number of volumes
    ignoreS = seriesNumbers(ismember(seriesVolumes,datasetConfigs.volumes) == 0);
    
    % More than one anatomical
    % This is assessed using the number of volumes of the first run
    % (anatomical).
    if sum(seriesVolumes == datasetConfigs.volumes(1)) > 1
        
        boolInput = false;
        disp(['[createFolderStructure] More than one run of anatomical data detected: ' num2str(seriesNumbers(seriesVolumes == datasetConfigs.volumes(1))')])
        while ~boolInput
            x = input('Please input the ones to ignore [<series numbers>]: ','s');
            
            if ~ismember(str2num(x),seriesNumbers(seriesVolumes == datasetConfigs.volumes(1)))
                disp('!---> ERROR: Incorrect series number.');
            else
                ignoreS = [ str2num(x) ignoreS ];
                boolInput = true;
            end
        end
        
    end
    
%     % More than one localiser
%     if sum(seriesVolumes == datasetConfigs.volumes(2)) > 1
%         
%         boolInput = false;
%         disp(['[createFolderStructure] More than one localiser run data detected: ' num2str(seriesNumbers(seriesVolumes == datasetConfigs.volumes(2))')])
%         while ~boolInput
%             x = input('Please input the ones to ignore [<series numbers>]: ','s');
%             
%             if ~ismember(str2double(x),seriesNumbers(seriesVolumes == datasetConfigs.volumes(2)))
%                 disp('!---> ERROR: Incorrect series number.');
%             else
%                 ignoreS = [ str2double(x) ignoreS ];
%                 boolInput = true;
%             end
%         end
%         
%     end
    
    disp(['[createFolderStructure] Ignoring files with series number of ' num2str(ignoreS')]);
    files(ismember(series,ignoreS)) = [];
    idx_to_delete = ismember(seriesNumbers,ignoreS);
    seriesNumbers(idx_to_delete) = [];
    seriesVolumes(idx_to_delete) = [];
    
    % Still more series than expected
    if length(seriesNumbers) > nRuns
        ignoreS = [];
        boolInput = false;
        disp(['[createFolderStructure] ' num2str(length(seriesNumbers) - nRuns) ' extra series remain.']);
        while ~boolInput
            disp(['[createFolderStructure] Current series: ' mat2str(seriesNumbers) '.'])
            x = input('Please input the ones to ignore [<series numbers>]: ','s');
            
            if length(str2num(x)) > length(seriesNumbers) - nRuns
                disp(['!---> ERROR: Too many series to delete. Choose only ' length(seriesNumbers) - nRuns ]);
            elseif ~ismember(str2num(x),seriesNumbers)
                disp('!---> ERROR: Incorrect series number.');
            else
                ignoreS = [ str2num(x) ignoreS ];
                boolInput = true;
            end
        end
        disp(['[createFolderStructure] Ignoring files with series number of ' num2str(ignoreS)]);
        files(ismember(series,ignoreS)) = [];
        idx_to_delete = ismember(seriesNumbers,ignoreS);
        seriesNumbers(idx_to_delete) = [];
        seriesVolumes(idx_to_delete) = [];
    end

% Number of series smaller than expected number of runs 
elseif length(seriesNumbers) < nRuns
    disp('[createFolderStructure] !---> ERROR: Unsufficient data.')
    boolInput = false;
    while ~boolInput
        x = input('[createFolderStructure] Do you wish to proceed anyway (Y/N)?','s');
        switch lower(x)
            case 'y'
                nRuns = length(seriesNumbers);
                boolInput = true;
            otherwise
                return
        end
    end
end

% -------------------------------------------------------------------------
% Check for incorrect number of volumes in all runs
% -------------------------------------------------------------------------

if any(datasetConfigs.volumes ~= seriesVolumes)
    disp('[createFolderStructure] Run volumes do not match the expected:');
    disp(['Expected:  ' num2str(datasetConfigs.volumes)])
    disp(['Input:     ' num2str(seriesVolumes)]);
    
    boolInput = false;
    while ~boolInput
        x = input('[createFolderStructure] Do you wish to proceed anyway (Y/N)?','s');
        switch lower(x)
            case 'y'
                boolInput = true;
            otherwise
                return
        end
    end
end

% -------------------------------------------------------------------------
% Check if Experiment folder exists
% -------------------------------------------------------------------------
if exist(datasetConfigs.path,'dir') ~= 7
    mkdir(datasetConfigs.path,'ANALYSIS');
    disp('[createFolderStructure] Experiment folder created.')
end

% -------------------------------------------------------------------------
% Create Subject Folder
% -------------------------------------------------------------------------
subjectFolder = fullfile(datasetConfigs.path,datasetConfigs.subjects{subjectIndex});

boolInput = false;

if exist(subjectFolder,'dir') == 7
    while ~boolInput
        x = input('[createFolderStructure] Subject Folder already exists. Do you wish to overwrite (Y), stop (N) or proceed (P)?','s');
        switch lower(x)
            case 'y'
                rmdir(subjectFolder,'s')
                mkdir(subjectFolder)
                boolInput = true;
            case 'n'
                return
            case 'p'
                success = true;
                return
            otherwise
                disp('[createFolderStructure] !---> ERROR: Invalid input.')
        end
    end
else
    mkdir(subjectFolder)
end


% -------------------------------------------------------------------------
% Iterate on the runs
% -------------------------------------------------------------------------
for r = 1:nRuns
    
    % Create folder structure
    switch datasetConfigs.runs{r}
        case 'anatomical'
            mkdir(subjectFolder,datasetConfigs.runs{r});
            dataFolder = fullfile(subjectFolder,datasetConfigs.runs{r},'DICOM');
            niifolder = fullfile(subjectFolder,datasetConfigs.runs{r});
        otherwise
            runFolderName = ['r-' datasetConfigs.runs{r}];
            mkdir(subjectFolder,runFolderName);
            mkdir(fullfile(subjectFolder,runFolderName),'DICOM'); % IMA files
            mkdir(fullfile(subjectFolder,runFolderName),'PRT');   % PRT files
            mkdir(fullfile(subjectFolder,runFolderName),'f');     % .nii
            mkdir(fullfile(subjectFolder,runFolderName),'af');    % .nii after STC
%             mkdir(fullfile(subjectFolder,['r-' datasetConfigs.runs{r}]),'raf');   % .nii after realign
            mkdir(fullfile(subjectFolder,runFolderName),'wraf');  % .nii after warp
            mkdir(fullfile(subjectFolder,runFolderName),'swraf'); % .nii after smooth
            
            mkdir(fullfile(subjectFolder,runFolderName),'ANALYSIS'); % GLM stuff

            niifolder = fullfile(subjectFolder,runFolderName,'f');
            dataFolder = fullfile(subjectFolder,runFolderName,'DICOM');
            
            % Copy Protocol (PRT) file
            prtGood = true;
            
            auxDir = dir(fullfile(dataTBV,[datasetConfigs.prtPrefix{r-1} '*.prt']));
            
            if ~isempty(auxDir) % If a .prt file was found
                
                keepDir = auxDir;
                if length(auxDir) > 1 % If more than one was found
                    boolInput = false;
                    disp(['!! More than one ' datasetConfigs.runs{r} ' run prt file detected !!'])
                    disp(extractfield(auxDir,'name')');
                    
                    while ~boolInput
                        x = input(sprintf('Please input the number (1:%i) of the PRT to keep: ',length(auxDir)),'s');
                        
                        if isnan(str2double(x)) || str2double(x) > length(auxDir)
                            fprintf('!---> ERROR: Please insert a number between 1 and %i/n.',length(auxDir));
                        else
                            keepDir = auxDir(str2double(x));
                            boolInput = true;
                        end
                    end
                end
                 
                copyfile(fullfile(keepDir.folder,keepDir.name),...
                    fullfile(subjectFolder,runFolderName,'PRT'));                
                
            else
                fprintf('!! No .prt file found in TBV folder for %s !! \n',datasetConfigs.runs{r});
                prtGood = false;
            end

            % Convert PRT to mat file
            if prtGood
                prtFile = xff(fullfile(subjectFolder,runFolderName,'PRT',keepDir.name));
                
                names = prtFile.ConditionNames';
                onsets = cell(1,length(names));
                durations = cell(1,length(names));
                onoff = prtFile.OnOffsets;
                
                for c = 1:length(names)                   
                    onsets{c} = onoff(onoff(:,1)==c,2)';
                    durations{c} = (onoff(onoff(:,1)==c,3)-onoff(onoff(:,1)==c,2)+1)';                   
                end
                
                save(fullfile(subjectFolder,runFolderName,'PRT',[keepDir.name(1:end-4) '.mat']),'names','onsets','durations')
            end
            
    end
        
    % Copy DICOM files of the series/run
    fprintf('Copying %s files...\n',datasetConfigs.runs{r});
    search_name = [auxnamesplit{1} '.' auxnamesplit{2} '.' num2str(seriesNumbers(r),'%.4i') '*'];
    copyfile( fullfile(dataPath,search_name) , dataFolder );
    
    % Convert to .nii
    auxdir = dir(fullfile(dataFolder,search_name));    
    fileList = cellfun(@(x) fullfile(auxdir(1).folder,x),{auxdir.name}','UniformOutput',false);
    convertDCMtoNII( fileList , niifolder );
    
    % Rename .nii files beautifully
%     oldFileNames = dir(fullfile(niifolder,'*.nii'));
%     niiname = sprintf('%s_%s_%04i',datasetConfigs.subjects{subjectIndex},datasetConfigs.runs{r},seriesNumbers(r));
%     
%     switch datasetConfigs.runs{r}
%         case 'anatomical'
%            movefile(fullfile(niifolder,oldFileNames.name),fullfile(niifolder,[niiname '.nii'])); 
%         otherwise
%            for jj = 1:numel(oldFileNames)
%                movefile(fullfile(niifolder,oldFileNames(jj).name),fullfile(niifolder,sprintf('%s_%04i.nii',niiname,jj)));
%            end
%     end

    if ~ strcmp(datasetConfigs.runs{r},'anatomical')
        dcmHeader = dicominfo(fullfile(auxdir(1).folder,auxdir(1).name));
        DCMinfo(idx_info).sliceTimes = dcmHeader.Private_0019_1029;
        DCMinfo(idx_info).sliceNumber = length(DCMinfo(idx_info).sliceTimes);
        [~,DCMinfo(idx_info).sliceVector] = sort(DCMinfo(idx_info).sliceTimes);
        DCMinfo(idx_info).TR = dcmHeader.RepetitionTime / 1000;
        DCMinfo(idx_info).TA = DCMinfo(idx_info).TR-(DCMinfo(idx_info).TR/DCMinfo(idx_info).sliceNumber);
        DCMinfo(idx_info).RefSlice = DCMinfo(idx_info).sliceVector(1);
        DCMinfo(idx_info).EchoTime = dcmHeader.EchoTime / 1000;
        
        idx_info = idx_info + 1;
    end

end

save(fullfile(subjectFolder,'DCMinfo.mat'),'DCMinfo');

mkdir(fullfile(subjectFolder,'batch'));
mkdir(fullfile(subjectFolder,'mr'));

success = true;
disp('[createFolderStructure] Folder structure creation completed.')

end

