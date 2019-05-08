clear ; clc, close all;

%% Load Stuff and Set Stuff

% --- Configuration data
load('Configs_VP_INHIBITION_SPM.mat')

% --- Subject Name
subjectName = 'VPIS07';
subjectIndex = find(not(cellfun('isempty', strfind(datasetConfigs.subjects, subjectName))));

% --- Select Run Sequence
datasetConfigs.runs = datasetConfigs.runs{subjectIndex};
% datasetConfigs.volumes = datasetConfigs.volumes{subjectIndex};
% datasetConfigs.prtPrefix = datasetConfigs.prtPrefix{subjectIndex};

% --- I/O Folders
if ismac
%     ioFolder = '/Users/alexandresayal/Google Drive/GitHub_DATA/ICNAS_VisualPerception/Inhibition/io_folder';
%     outputFolder = '/Users/alexandresayal/Google Drive/GitHub_DATA/ICNAS_VisualPerception/Inhibition/output_images';
%     datasetConfigs.path = '/Volumes/E/DATA_VP_Inhibition';
else
    ioFolder = fullfile('F:\Google Drive\GitHub_DATA\ICNAS_VisualPerception\Inhibition_SPM','io_folder');
    outputFolder = fullfile('F:\Google Drive\GitHub_DATA\ICNAS_VisualPerception\Inhibition_SPM','images');
end

% --- VOI File
% voiFilePath = fullfile(datasetConfigs.path,subjectName,'MultiRun');
% voiFileName = 'VOIs_Loc_edit.voi';
% voiFile = xff(fullfile(voiFilePath,voiFileName));

% --- Runs to process
selectedRuns = datasetConfigs.runs(2:end);
nTrialsPerRun = 2;
nTrials = length(selectedRuns) * nTrialsPerRun;

% --- Structures to keep data
TCP = struct(); % TimeCourse and Protocol data
% TCP.MultiRun.Baseline = zeros(1,voiFile.NrOfVOIs);
% TCP.MultiRun.Average = zeros(1,voiFile.NrOfVOIs);
ERA = struct(); % ERA data

% --- Haemodynamic Delay
delay = 3; %in volumes

%% PRTs and VOI TimeCourses

for r = selectedRuns
    
    % --- Read Protocol
    prtFilePath = fullfile(datasetConfigs.path,subjectName,['r-' r{:}],'PRT');
    prtFileName = ['RunMRI_D12_R' r{:}(end) '.prt'];    
    [ cond_names , TCP.(r{:}).intervalsPRT , TCP.(r{:}).intervals ] = readProtocol( prtFilePath , prtFileName , datasetConfigs.TR );
    
end

%     
%     % --- Read VTC
%     vtcFileName = [subjectName '_' r{:} '_SCCTBL_3DMCTS_LTR_THPGLMF2c_MNI_SD3DVSS5.00mm.vtc'];
%     vtcFile = xff(fullfile(prtFilePath,vtcFileName));
%     
%     % --- Read VOIs TC
%     TCP.(r{:}).VOITC = vtcFile.VOITimeCourseOrig(voiFile);
%     
%     % --- Calculate Baseline Value (average of static periods)
%     baseVols = TCP.(r{:}).intervals==find(ismember(cond_names, 'Static'));
%     baseVols = [ zeros(1,delay,'logical') baseVols(1:end-delay) ];
%     TCP.(r{:}).Baseline = mean(TCP.(r{:}).VOITC(baseVols,:),1);
%     
%     TCP.MultiRun.Baseline = TCP.MultiRun.Baseline + TCP.(r{:}).Baseline;
%     
%     % --- Calculate average of TC
%     TCP.(r{:}).Average = mean(TCP.(r{:}).VOITC(delay+1:end,:),1);
%     
%     TCP.MultiRun.Average = TCP.MultiRun.Average + TCP.(r{:}).Average;
%     
%     % --- Memory is valuable
%     vtcFile.ClearObject;
%     
% end
% 
% TCP.MultiRun.Baseline = TCP.MultiRun.Baseline / length(selectedRuns);
% TCP.MultiRun.Average = TCP.MultiRun.Average / length(selectedRuns);

%% Load TC data and protocols
load(fullfile(datasetConfigs.AnalysisPath,'TCs-MT-5mm.mat'))

voiNames = {'leftMT','rightMT','bilateralMT'};

%% Create ERA

% --- Trial Types
% trialTps = cond_names(2:4);
trialTestTps = cond_names(5:10);

% ---
extraB = 3;

% --- Iterate
for v = 1:nRois
        
    for t = 1:length(trialTestTps)
        
        idx = 1;

        % --- Initialise matrices.
        nTrialVols = 30; % 3 (static) + 12 (adaptation) + 6 (test) + 3 (MAE) + 3 (report) + 3 (static) = 30 volumes
        ERA.(voiNames{v}).(trialTestTps{t}).data = zeros(nTrials,nTrialVols);
        % For the mean (nTrials+1), sem (nTrials+2), median (nTrials+3) and semedian (nTrials+4) of the above trials.
        ERA.(voiNames{v}).(trialTestTps{t}).stats = zeros(4,nTrialVols);
        
        r_idx = 1;
        for r = selectedRuns
            
            % --- Volumes of interest
            % I know there are two trials per run, so I spare a for loop :)
            int_aux = [];
            int_aux(1,:) = TCP.(r{:}).intervalsPRT.(trialTestTps{t})(1,1) - 12 - extraB : TCP.(r{:}).intervalsPRT.(trialTestTps{t})(1,2) + 3 + 3 + extraB ;
            int_aux(2,:) = TCP.(r{:}).intervalsPRT.(trialTestTps{t})(2,1)- 12 - extraB : TCP.(r{:}).intervalsPRT.(trialTestTps{t})(2,2) + 3 + 3 + extraB ;

            % --- Compensate Haemodynamic Delay
            int_aux = int_aux + delay;
                        
            % --- Calculate BOLD Signal Var (Baseline per run)
            ERA.(voiNames{v}).(trialTestTps{t}).data(idx+0,:) = pscArray{subjectIndex,r_idx,v}(int_aux(1,:)) * 100;
            
            ERA.(voiNames{v}).(trialTestTps{t}).data(idx+1,:) = pscArray{subjectIndex,r_idx,v}(int_aux(2,:)) * 100;
                       
            idx = idx + nTrialsPerRun;
            r_idx = r_idx + 1;
            
        end
        
        % --- Mean of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).stats(1,:) = mean(ERA.(voiNames{v}).(trialTestTps{t}).data(1:idx-1,:));
        
        % --- Std Error of the Mean of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).stats(2,:) = std(ERA.(voiNames{v}).(trialTestTps{t}).data(1:idx-1,:)) / sqrt(nTrials); 
        
        % --- Median of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).stats(3,:) = median(ERA.(voiNames{v}).(trialTestTps{t}).data(1:idx-1,:));
        
        % --- Std Error of the Median of Trials
        ERA.(voiNames{v}).(trialTestTps{t}).stats(4,:) = 1.253*std(ERA.(voiNames{v}).(trialTestTps{t}).data(1:idx-1,:)) / sqrt(nTrials);
        
    end
    
end

%% Plot Averages for Bilateral MT - MEAN - Baseline per Run

fig1 = figure('Name','ERA','Position',[100 100 1300 1000]);
xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh aCoh';'Coh aInCoh';'Coh aNA';'InCoh aCoh';'InCoh aInCoh';'InCoh aNA'};
% ------------------------------------------------------------------------%
subplot(2,1,1)
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).stats(2,:));
    end
    
    hold on;    
end
hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal variation (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'Coherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([1 2 3]))
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
subplot(2,1,2)
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA.bilateralMT.(trialTestTps{tt}).stats(1,:),...
        ERA.bilateralMT.(trialTestTps{tt}).stats(2,:),...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA.bilateralMT.(trialTestTps{tt}).stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).stats(2,:)) > absMax
        absMax = max(ERA.bilateralMT.(trialTestTps{tt}).stats(1,:)+ERA.bilateralMT.(trialTestTps{tt}).stats(2,:));
    end
    
    hold on;    
end
hold off
xx = [xvector(1)-1 xvector(end)+1];
yy = [-1 ceil(absMax+0.5)];
xlim(xx); ylim(yy);
xlabel('Time (volumes)')
ylabel({'BOLD signal variation (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.2,'Motion','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.2,'InCoherent','HorizontalAlignment','center','FontSize',12)
text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([4 5 6]))
% ------------------------------------------------------------------------%

suptitle(sprintf('ERA of Bilateral MT - Subject %s',subjectName))

%% Save above figure
% export_fig(fullfile(outputFolder,[ subjectName '_FullTC_BilateralMT_D' num2str(delay)]),'-eps','-pdf','-transparent');
print(fig1,fullfile(outputFolder,[ subjectName '_FullTC_BilateralMT_D' num2str(delay)]),'-dpng')

%% Calculate AuC of Test blocks
x_auc = 15:22; % test block + 1 before + 1 after
Y_AUC = struct();
absMin13 = 99;
absMin46 = 99;

for tt = 1:3
    
    % --- Retrive TC values and SEM values
    Y_AUC.(trialTestTps{tt}).psc = ERA.bilateralMT.(trialTestTps{tt}).stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc_sem = ERA.bilateralMT.(trialTestTps{tt}).stats(2,x_auc);
    
    % --- Retrieve TC values normalising with the average of the three last
    % points before test of each trial type
    Y_AUC.(trialTestTps{tt}).psc_norm = Y_AUC.(trialTestTps{tt}).psc - mean(ERA.bilateralMT.(trialTestTps{tt}).stats(1,x_auc(1)-2:x_auc(1)));
    
    % --- Find minimum
    if min(Y_AUC.(trialTestTps{tt}).psc_norm) < absMin13
        absMin13 = min(Y_AUC.(trialTestTps{tt}).psc_norm);
    end
    
end

for tt = 4:6
    
    % --- Retrive TC values and SEM values
    Y_AUC.(trialTestTps{tt}).psc = ERA.bilateralMT.(trialTestTps{tt}).stats(1,x_auc);
    Y_AUC.(trialTestTps{tt}).psc_sem = ERA.bilateralMT.(trialTestTps{tt}).stats(2,x_auc);
    
    % --- Retrieve TC values normalising with the average of the three last
    % points before test of each trial type
    Y_AUC.(trialTestTps{tt}).psc_norm = Y_AUC.(trialTestTps{tt}).psc - mean(ERA.bilateralMT.(trialTestTps{tt}).stats(1,x_auc(1)-2:x_auc(1)));
    
    % --- Find minimum
    if min(Y_AUC.(trialTestTps{tt}).psc_norm) < absMin46
        absMin46 = min(Y_AUC.(trialTestTps{tt}).psc_norm);
    end
    
end

%% Remove minimum value of all trials
% This will make the minimum point of the trials = 0 --> push down the
% curves to avoid unnecessary area :)

for tt = 1:3
    
    Y_AUC.(trialTestTps{tt}).psc_norm0 = Y_AUC.(trialTestTps{tt}).psc_norm - absMin13;  
    
end

for tt = 4:6
    
    Y_AUC.(trialTestTps{tt}).psc_norm0 = Y_AUC.(trialTestTps{tt}).psc_norm - absMin46;  
    
end

%% Prepare The ultimate plot
AUC = struct();
trialLabels = {'Coh aCoh';'Coh aInCoh';'Coh aNA';'InCoh aCoh';'InCoh aInCoh';'InCoh aNA'};
clrMap = lines;
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-0.2 1.2];
lwidth = 1.5;

% import beta values
% load(fullfile(ioFolder,'BETAS_test_persubject.mat'))
% betas = BETAS(subjectIndex,:);

%% Plot The ultimate plot
fig_ultimate = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.5 0.8]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(2,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.Coh_aCoh.psc_norm0,Y_AUC.Coh_aCoh.psc_sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.Coh_aInCoh.psc_norm0,Y_AUC.Coh_aInCoh.psc_sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.Coh_aNA.psc_norm0,Y_AUC.Coh_aNA.psc_sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal variation (%)'},'FontSize',12)
title('ERA')
s1.FontSize = 12;
% ------------------------------------------------------------------------%
% -- AUC Coherent --------------------------------------------------------%
s2 = subplot(2,2,2);
yy_auc = 0;
for ii = 1:3
    AUC.(trialTestTps{ii}) = trapz(x_auc,Y_AUC.(trialTestTps{ii}).psc_norm0);
    
    b = bar(ii,AUC.(trialTestTps{ii}),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
    text(ii,AUC.(trialTestTps{ii})+0.25,num2str(AUC.(trialTestTps{ii})),...
        'HorizontalAlignment','center','FontSize',12)
    hold on
    
    if AUC.(trialTestTps{ii}) > yy_auc
        yy_auc = AUC.(trialTestTps{ii});
    end
    
end
hold off

xticks(1:3)
xticklabels(trialLabels(1:3))
ylim([0 ceil(yy_auc)+1])
ylabel('Area','FontSize',12)
title('Area under Curve','FontSize',12)
s2.FontSize = 12;
% ------------------------------------------------------------------------%
% -- Betas Coherent ------------------------------------------------------%
% s3 = subplot(2,3,3);
% yy_auc = 0;
% for ii = 1:3
%     
%     b = bar(ii,betas(ii+3),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
%     text(ii,betas(ii+3)+0.25,num2str(betas(ii+3)),'HorizontalAlignment','center',...
%         'FontSize',12)
%     hold on
%     
%     if betas(ii+3) > yy_auc
%         yy_auc = betas(ii+3);
%     end
%     
% end
% hold off
% 
% xticks(1:3)
% xticklabels(trialLabels(1:3))
% ylim([0 ceil(yy_auc)+1])
% ylabel('beta value')
% title('Betas','FontSize',12)
% s3.FontSize = 12;
% ------------------------------------------------------------------------%
% -- ERA InCoherent ------------------------------------------------------%
s4 = subplot(2,2,3);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC.InCoh_aCoh.psc_norm0,Y_AUC.InCoh_aCoh.psc_sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC.InCoh_aInCoh.psc_norm0,Y_AUC.InCoh_aInCoh.psc_sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC.InCoh_aNA.psc_norm0,Y_AUC.InCoh_aNA.psc_sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal variation (%)'},'FontSize',12)
s4.FontSize = 12;
% ------------------------------------------------------------------------%
% -- AUC InCoherent ------------------------------------------------------%
s5 = subplot(2,2,4);
yy_auc = 0;
for ii = 4:6
    AUC.(trialTestTps{ii}) = trapz(x_auc,Y_AUC.(trialTestTps{ii}).psc_norm0);
    
    b = bar(ii,AUC.(trialTestTps{ii}),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
    text(ii,AUC.(trialTestTps{ii})+0.25,num2str(AUC.(trialTestTps{ii})),...
        'HorizontalAlignment','center','FontSize',12)
    hold on
    
    if AUC.(trialTestTps{ii}) > yy_auc
        yy_auc = AUC.(trialTestTps{ii});
    end
    
end
hold off

xticks(4:6)
xticklabels(trialLabels(4:6))
ylim([0 ceil(yy_auc)+1])
ylabel('Area under curve')
s5.FontSize = 12;
% ------------------------------------------------------------------------%
% -- Betas InCoherent ----------------------------------------------------%
% s6 = subplot(2,3,6);
% yy_auc = 0;
% for ii = 4:6
%     
%     b = bar(ii,betas(ii+3),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
%     text(ii,betas(ii+3)+0.25,num2str(betas(ii+3)),'HorizontalAlignment','center',...
%         'FontSize',12)
%     hold on
%     
%     if betas(ii+3) > yy_auc
%         yy_auc = betas(ii+3);
%     end
%     
% end
% hold off
% 
% xticks(4:6)
% xticklabels(trialLabels(4:6))
% ylim([0 ceil(yy_auc)+1])
% ylabel('beta value')
% s6.FontSize = 12;
% ------------------------------------------------------------------------%
suptitle(sprintf('ERA and AUC - Subject %s',subjectName))

%% Export The ultimate Figure
% export_fig(fullfile(outputFolder,[ subjectName '_FullTC_BilateralMT_D' num2str(delay)]),'-eps','-pdf','-transparent');
print(fig_ultimate,fullfile(outputFolder,[ subjectName '_UltimateFigure_BilateralMT_D' num2str(delay)]),'-dpng')

%% Save dataset
save(fullfile(ioFolder,[subjectName '_bold.mat']),'ERA','Y_AUC','trialLabels','AUC','x_auc','trialLabels','clrMap','trialTestTps','nTrialVols','subjectIndex','subjectName')
    