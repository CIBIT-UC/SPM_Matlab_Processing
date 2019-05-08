
marsbar('on')

load('Configs_VP_INHIBITION_SPM.mat')

runs = datasetConfigs.runs{1,1}(2:end);
nRuns = length(runs);
nSubjects = length(datasetConfigs.subjects);
nRois = 2 + 1;
delay = 3; % haemodynamic delay in volumes

roiSizeTag = '5mm';

tcArray = cell(nSubjects,nRuns,nRois); % TCs from all subjects, all runs, all rois

pscArray = cell(nSubjects,nRuns,nRois); % PSC from all subjects, all runs, all rois

% Iterate on the subjects
for subjectIndex = 1:nSubjects
    
    subjectName = datasetConfigs.subjects{subjectIndex};
    
    dataRootSubject = fullfile(datasetConfigs.path,subjectName);
    
    roiAnalysisPath = fullfile(dataRootSubject,['r-' runs{1}],'ANALYSIS');
    
    clear roiArray
    roiArray{1} = maroi(fullfile(roiAnalysisPath,[subjectName '-leftMT-' roiSizeTag '_roi.mat' ]));
    roiArray{2} = maroi(fullfile(roiAnalysisPath,[subjectName '-rightMT-' roiSizeTag '_roi.mat' ]));
    roiArray{3} = label( roiArray{1} | roiArray{2} , ['bilateralMT-' roiSizeTag]); % bilateral MT
    
    % Iterate on the runs
    for runIndex = 1:nRuns
        runAnalysisPath = fullfile(dataRootSubject,['r-' runs{runIndex}],'ANALYSIS');
        
        % Protocol - baseline indexes
        prtmatFile = fullfile(dataRootSubject,['r-' runs{runIndex}],'PRT',['RunMRI_D12_R' num2str(runIndex) '.mat']);
        load(prtmatFile)
        offsets = (onsets{1,1}+durations{1,1}-1);
        onsets = onsets{1,1};
        [t1, t2] = ndgrid(onsets, 0:(offsets(1)-onsets(1)));
        baselineVols = t1 + t2 + delay;
        
        % SPM design file
        D = mardo(fullfile(runAnalysisPath,'SPM.mat'));
        
        for roiIndex = 1:nRois
            
            % Extract data
            Y = get_marsy(roiArray{roiIndex}, D, 'mean'); % OR 'eigen1'
            
            auxTC = summary_data(Y);
            tcArray{subjectIndex,runIndex,roiIndex} = auxTC;  % get summary time course(s)
            
            baselineMean = mean(mean(auxTC(baselineVols)));
            pscArray{subjectIndex,runIndex,roiIndex} = (auxTC - baselineMean) / baselineMean; % calculate PSC
            
        end
        
    end
    
end

save(fullfile(datasetConfigs.AnalysisPath,'TCs-MT-5mm.mat'),'tcArray','pscArray','nSubjects','nRois','nRuns','delay')
