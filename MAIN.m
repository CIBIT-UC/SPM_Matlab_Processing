tic
clear, clc;
spm('Defaults','fMRI');
spm_jobman('initcfg');

spmInstallPath = 'C:\Users\alexandresayal\Documents\MATLAB\spm12';

% Force folder
cd('C:\Users\alexandresayal\Documents\GitHub\SPM_Matlab_Processing');

%% ============================= Settings ============================== %%

load('Configs_VP_MBEPI_SPM.mat')

dataPath = 'F:\RAW_DATA_MB-EPI_TESTS\TESTEMB_03OUT2018';
dataTBV = 'F:\RAW_DATA_MB-EPI_TESTS\PRT'; 

%% ============================= Settings ============================== %%
subjectName = 'Pilot03_2';  % Subject name
%=========================================================================%

subjectIndex = find(not(cellfun('isempty', strfind(datasetConfigs.subjects, subjectName))));

if isempty(subjectIndex)
    error('Invalid Subject Name')
end

% --- Select Run Sequence
datasetConfigs.runs = datasetConfigs.runs{subjectIndex};
datasetConfigs.volumes = datasetConfigs.volumes{subjectIndex};
datasetConfigs.prtPrefix = datasetConfigs.prtPrefix{subjectIndex};

dataRootSubject = fullfile(datasetConfigs.path, subjectName);

%% -- Create Folder Structure
[ success , DCMinfo ] = createFolderStructure( datasetConfigs , dataPath , dataTBV , subjectIndex );

assert(success,'Folder creation aborted.');
clear dataPath dataTBV;

%% Retrieve functional runs info
functionalRuns = ( dir ( fullfile( dataRootSubject, 'r-*' ) ) );
numFunctionalRuns = length( functionalRuns );
alignRunIdx = 1;
alignRun = functionalRuns(alignRunIdx).name; % for coregister

%% Retrieve functional files cell
scansList = cell(1,numFunctionalRuns);
prefix = 'f';
add = ',1';

for r = 1:numFunctionalRuns
    auxDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'f',[prefix '*.nii']));
    scansList{1,r} = cellfun(@(x) fullfile(auxDir(1).folder,[x add]),{auxDir.name}','UniformOutput',false);
end

%% Retrieve structural file cell
scansListAnatomical = cell(1,1);

auxDir =  dir(fullfile(dataRootSubject,'anatomical','s*.nii'));
scansListAnatomical{1,1} = fullfile(auxDir(1).folder,[auxDir(1).name ',1']);

%% Change directory
% Not sure yet if it is useful or not
% If functions are needed they will need to be added to the matlab path
cd(fullfile(dataRootSubject,'batch'));
diary(sprintf('diary_%s_%s.txt',subjectName,datestr(now,'ddmmmmyyyy_HHMM')))

%% -- Run step 1
% Run preprocess_1 batch script
% Performs slice time correction and motion correction

% Create preprocess_1.mat
matlabbatch = cell(1,2);
matlabbatch{1}.spm.temporal.st.scans = scansList;
matlabbatch{1}.spm.temporal.st.nslices = DCMinfo.sliceNumber;
matlabbatch{1}.spm.temporal.st.tr = DCMinfo.TR;
matlabbatch{1}.spm.temporal.st.ta = DCMinfo.TA;
matlabbatch{1}.spm.temporal.st.so = DCMinfo.sliceVector;
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';
for r = 1:numFunctionalRuns
    matlabbatch{2}.spm.spatial.realign.estimate.data{r}(1) = cfg_dep(['Slice Timing: Slice Timing Corr. Images (Sess ' num2str(r) ')'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{r}, '.','files'));
end
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.rtm = 0;
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estimate.eoptions.weight = '';

save(fullfile(dataRootSubject,'batch','preprocess_1.mat'),'matlabbatch');

% Load and run
load(fullfile(dataRootSubject,'batch','preprocess_1.mat'));
spm_jobman('run', matlabbatch);

%% -- Organise files
for r = 1:numFunctionalRuns
    try
        movefile(...
            fullfile(dataRootSubject,functionalRuns(r).name,'f','af*.nii'),...
            fullfile(dataRootSubject,functionalRuns(r).name,'af'))
    catch
        fprintf('No files found - %s \n',functionalRuns(r).name);
    end
end

%% Retrieve functional files cell 2
scansList2 = cell(1,numFunctionalRuns);
scansListFull = cell(0,1);
prefix = 'af';
add = ',1';

for r = 1:numFunctionalRuns
    auxDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'af',[prefix '*.nii']));
    scansList2{1,r} = cellfun(@(x) fullfile(auxDir(1).folder,[x add]),{auxDir.name}','UniformOutput',false);
    scansListFull = [scansListFull ; scansList2{1,r}];
end

%% -- Run step 2
% Run preprocess_2 batch script
% Performs coregistration, segmentation, normalisation of anatomical and functional images

clear matlabbatch

matlabbatch{1}.spm.spatial.coreg.estimate.ref = scansListAnatomical;
matlabbatch{1}.spm.spatial.coreg.estimate.source = scansList2{1,alignRunIdx}(1,1);
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

matlabbatch{2}.spm.spatial.preproc.channel.vols = scansListAnatomical;
matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{2}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {[fullfile(spmInstallPath,'tpm','TPM.nii') ',1']};
matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {[fullfile(spmInstallPath,'tpm','TPM.nii') ',2']};
matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm =  {[fullfile(spmInstallPath,'tpm','TPM.nii') ',3']};
matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm =  {[fullfile(spmInstallPath,'tpm','TPM.nii') ',4']};
matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm =  {[fullfile(spmInstallPath,'tpm','TPM.nii') ',5']};
matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm =  {[fullfile(spmInstallPath,'tpm','TPM.nii') ',6']};
matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{2}.spm.spatial.preproc.warp.write = [1 1];

matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol = scansListAnatomical;
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = scansListAnatomical;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spmInstallPath,'tpm','TPM.nii')};
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                                                       78 76 85];
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.interp = 7;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = scansList2{1,alignRunIdx}(1,1);
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample = scansListFull;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spmInstallPath,'tpm','TPM.nii')};
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                                                       78 76 85];
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.vox = [3 3 3];
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.prefix = 'wr';

save(fullfile(dataRootSubject,'batch','preprocess_2.mat'),'matlabbatch');

% Load and run
load(fullfile(dataRootSubject,'batch','preprocess_2.mat'));
spm_jobman('run', matlabbatch);

%% Organise and Retrieve files
scansList3 = cell(1,numFunctionalRuns);
scansListFull2 = cell(0,1);
prefix = 'wraf';
add = ',1';

for r = 1:numFunctionalRuns
        try
        movefile(...
            fullfile(dataRootSubject,functionalRuns(r).name,'af','wraf*.nii'),...
            fullfile(dataRootSubject,functionalRuns(r).name,'wraf'))
    catch
        fprintf('No files found - %s \n',functionalRuns(r).name);
        end
    
    auxDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'wraf',[prefix '*.nii']));
    scansList3{1,r} = cellfun(@(x) fullfile(auxDir(1).folder,[x add]),{auxDir.name}','UniformOutput',false);
    scansListFull2 = [scansListFull2 ; scansList3{1,r}];
end

%% Run step 3 - Smooting
% Run preprocess_3 batch script
% Performs smoothing of functional images

clear matlabbatch

matlabbatch{1}.spm.spatial.smooth.data = scansListFull2;
matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

save(fullfile(dataRootSubject,'batch','preprocess_3.mat'),'matlabbatch');

% Load and run
load(fullfile(dataRootSubject,'batch','preprocess_3.mat'));
spm_jobman('run', matlabbatch);

%% -- Organise files
for r = 1:numFunctionalRuns
    try
        movefile(...
            fullfile(dataRootSubject,functionalRuns(r).name,'wraf','swraf*.nii'),...
            fullfile(dataRootSubject,functionalRuns(r).name,'swraf'))
    catch
        fprintf('No files found - %s \n',functionalRuns(r).name);
    end
end

%% -- Retrieve files - after all processing steps
scansList_F = cell(1,numFunctionalRuns);
scansListFull_F = cell(0,1);
prefix = 'swraf';
add = ',1';

for r = 1:numFunctionalRuns    
    
    auxDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'swraf',[prefix '*.nii']));
    scansList_F{1,r} = cellfun(@(x) fullfile(auxDir(1).folder,[x add]),{auxDir.name}','UniformOutput',false);
    scansListFull_F = [scansListFull_F ; scansList_F{1,r}];
    
end

%% -- Single run analysis
% - Model specification and estimation

for r = 1:numFunctionalRuns    
    clear matlabbatch
    prtmatDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'PRT','*.mat'));
    motionRegDir = dir(fullfile(dataRootSubject,functionalRuns(r).name,'f','rp*.txt'));
    fname = strsplit(functionalRuns(r).name,'-');
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(dataRootSubject,functionalRuns(r).name,'ANALYSIS')};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = DCMinfo.sliceNumber;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scansList_F{1,r};
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(prtmatDir.folder,prtmatDir.name)};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(motionRegDir.folder,motionRegDir.name)};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    save(fullfile(dataRootSubject,'batch',['analysis_' fname{2} '.mat']),'matlabbatch');

    % Load and run
    load(fullfile(dataRootSubject,'batch',['analysis_' fname{2} '.mat']));
    spm_jobman('run', matlabbatch);

end

%% -- End main stuff
elapsedTime = toc;
clear matlabbatch
disp('Done');
beep;
save(sprintf('workspace_%s_%s',subjectName,datestr(now,'ddmmmmyyyy_HHMM')));
diary off

%% -- CUSTOM STUFF -- %%

%% -- Create localiser contrast
clear matlabbatch

matlabbatch{1}.spm.stats.con.spmmat = {fullfile(dataRootSubject,'r-loc','ANALYSIS','SPM.mat')};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Motion > Static';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 -1 1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch);

%% -- Create runsB1234 contrasts

for r = 2:numFunctionalRuns
    clear matlabbatch
    
    load(fullfile(dataRootSubject,functionalRuns(r).name,'ANALYSIS','SPM.mat'));
    
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(dataRootSubject,functionalRuns(r).name,'ANALYSIS','SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Real Motion > Static';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = define_contrasts(SPM,{'Adapt_','NonAdapt'},{'Static'});
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'NonAdapt > Adapt';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = define_contrasts(SPM,{'NonAdapt'},{'Adapt_'});
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    spm_jobman('run', matlabbatch);

end

%% -- Create runsB1234 FFX design
clear matlabbatch

if ~exist(fullfile(dataRootSubject,'mr','ffx_RunB1234'),'dir')
    mkdir(fullfile(dataRootSubject,'mr','ffx_RunB1234'));
end

matlabbatch{1}.spm.stats.mfx.ffx.dir = {fullfile(dataRootSubject,'mr','ffx_RunB1234')};
matlabbatch{1}.spm.stats.mfx.ffx.spmmat = {
                                           fullfile(dataRootSubject,'r-runB1','ANALYSIS','SPM.mat')
                                           fullfile(dataRootSubject,'r-runB2','ANALYSIS','SPM.mat')
                                           fullfile(dataRootSubject,'r-runB3','ANALYSIS','SPM.mat')
                                           fullfile(dataRootSubject,'r-runB4','ANALYSIS','SPM.mat')
                                           };
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('FFX Specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch);

%% -- Create contrasts for runsB1234 FFX
clear matlabbatch
    
load(fullfile(dataRootSubject,'mr','ffx_RunB1234','SPM.mat'));

matlabbatch{1}.spm.stats.con.spmmat = {fullfile(dataRootSubject,'mr','ffx_RunB1234','SPM.mat')};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Real Motion > Static';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = define_contrasts(SPM,{'Adapt_','NonAdapt'},{'Static'});
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'NonAdapt > Adapt';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = define_contrasts(SPM,{'NonAdapt'},{'Adapt_'});
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch);
