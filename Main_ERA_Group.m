clear,clc,close all

%%
N = 20; % number of subjects

%% Initialize stuff

% --- I/O Folders
if ismac
%     ioFolder = '/Users/alexandresayal/Google Drive/GitHub_DATA/ICNAS_VisualPerception/Inhibition/io_folder';
%     outputFolder = '/Users/alexandresayal/Google Drive/GitHub_DATA/ICNAS_VisualPerception/Inhibition/output_images';
%     datasetConfigs.path = '/Volumes/E/DATA_VP_Inhibition';
else
    ioFolder = fullfile('F:\Google Drive\GitHub_DATA\ICNAS_VisualPerception\Inhibition_SPM','io_folder');
    outputFolder = fullfile('F:\Google Drive\GitHub_DATA\ICNAS_VisualPerception\Inhibition_SPM','images_group');
end

% --- Import beta values
% load(fullfile(ioFolder,'BETAS_test_persubject.mat'))

delay = 3;

ERA_G = struct();

AUC_G = struct();

Y_AUC_G = struct();

for ss = 1:N
    
    load(fullfile(ioFolder,sprintf('VPIS%.2i_bold.mat',ss)))
    
    for tt = 1:length(trialTestTps)
    
        ERA_G.BilateralMT.(trialTestTps{tt}).mean(ss,:) = ERA.bilateralMT.(trialTestTps{tt}).stats(1,:);
        
        Y_AUC_G.(trialTestTps{tt}).data(ss,:) = Y_AUC.(trialTestTps{tt}).psc_norm0;
        
        AUC_G.(trialTestTps{tt}).data(ss,:) = AUC.(trialTestTps{tt});
        
    end
       
end

%% Calculate stats
for tt = 1:length(trialTestTps)
    
    ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean = mean(ERA_G.BilateralMT.(trialTestTps{tt}).mean);
    ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem = std(ERA_G.BilateralMT.(trialTestTps{tt}).mean) / sqrt(N);
    
    Y_AUC_G.(trialTestTps{tt}).stats.mean = mean(Y_AUC_G.(trialTestTps{tt}).data);
    Y_AUC_G.(trialTestTps{tt}).stats.sem = std(Y_AUC_G.(trialTestTps{tt}).data) / sqrt(N);
    
    AUC_G.(trialTestTps{tt}).stats.mean = mean(AUC_G.(trialTestTps{tt}).data);
    AUC_G.(trialTestTps{tt}).stats.sem = std(AUC_G.(trialTestTps{tt}).data) / sqrt(N);
    
end

%% Betas and differences

% BETAS_G.stats.mean = mean(BETAS);
% BETAS_G.stats.sem = std(BETAS) / sqrt(N);
% 
% % Calculate difference between betas
% BETAS_G.diff.labels = {'Coh aCoh - Coh','Coh aInCoh - InCoh','Coh aNa - Na','InCoh aCoh - Coh','InCoh aInCoh - InCoh','InCoh aNa - Na'};
% BETAS_G.diff.data = [BETAS(:,4)-BETAS(:,1) , BETAS(:,5)-BETAS(:,2) , BETAS(:,6)-BETAS(:,3) , ...
%                      BETAS(:,7)-BETAS(:,1) , BETAS(:,8)-BETAS(:,2) , BETAS(:,9)-BETAS(:,3) ];
%                  
% BETAS_G.diff.stats.mean = mean(BETAS_G.diff.data); 
% BETAS_G.diff.stats.sem = std(BETAS_G.diff.data) / sqrt(N); 

%% Ssafety measures
clear ERA BETAS AUC Y_AUC

%% Plot Averages for Bilateral MT - MEAN - Baseline per Run

fig1 = figure('Name','Group ERA','Position',[100 100 1300 1000]);
xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
% ------------------------------------------------------------------------%
s1=subplot(2,1,1);
absMax = 0;
for tt = [1 2 3]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem);
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
s1.FontSize = 16;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt,:));
    
    if max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem) > absMax
        absMax = max(ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean+ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem);
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
s2.FontSize = 16;
% ------------------------------------------------------------------------%

suptitle(sprintf('Group ERA of Bilateral MT - N = %i',N))

%% Save above figure
% export_fig(fullfile(outputFolder,sprintf('Group_N%i_FullTC_BilateralMT_D%i',N,delay))),'-eps','-pdf','-transparent');
print(fig1,fullfile(outputFolder,sprintf('Group_N%i_FullTC_BilateralMT_D%i',N,delay)),'-dpng')

%% Prepare The ultimate plot
trialLabels = {'Coh \rightarrow Coh';'InCoh \rightarrow Coh';'NA \rightarrow Coh';'Coh \rightarrow InCoh';'InCoh \rightarrow InCoh';'NA \rightarrow InCoh'};
clrMap = lines;
x_auc = 15:22; % test block + 1 before + 1 after
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-0.2 1.2];
lwidth = 1.5;

%% Plot The ultimate plot
fig_ultimate = figure('Name','The Ultimate Figure','units','normalized','outerposition',[0 0 0.8 0.9]);
movegui('center')

% ------------------------------------------------------------------------%
% -- ERA Coherent --------------------------------------------------------%
s1 = subplot(2,3,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.Coh_aCoh.stats.mean,Y_AUC_G.Coh_aCoh.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.Coh_aInCoh.stats.mean,Y_AUC_G.Coh_aInCoh.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.Coh_aNA.stats.mean,Y_AUC_G.Coh_aNA.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal variation (%)'},'FontSize',12)
title('ERA')
s1.FontSize = 14;
% ------------------------------------------------------------------------%
% -- AUC Coherent --------------------------------------------------------%
s2 = subplot(2,3,2);
yy_auc = 0;
for ii = 1:3
    
    errorbar(ii,AUC_G.(trialTestTps{ii}).stats.mean,AUC_G.(trialTestTps{ii}).stats.sem,...
        'Marker','none','Color',clrMap(ii,:),'LineWidth',1.5)
    
    hold on
    
    b = bar(ii,AUC_G.(trialTestTps{ii}).stats.mean,'FaceColor',clrMap(ii,:),'BarWidth',0.4);
    
    text(ii,AUC_G.(trialTestTps{ii}).stats.mean+AUC_G.(trialTestTps{ii}).stats.sem+0.25,num2str(AUC_G.(trialTestTps{ii}).stats.mean),...
        'HorizontalAlignment','center','FontSize',12)
    
    hold on
    
    if AUC_G.(trialTestTps{ii}).stats.mean > yy_auc
        yy_auc = AUC_G.(trialTestTps{ii}).stats.mean;
    end
    
end
hold off

xticks(1:3)
xticklabels(trialLabels(1:3))
ylim([0 ceil(yy_auc)+1])
ylabel('Area','FontSize',12)
title('Area under Curve','FontSize',12)
s2.FontSize = 14;
% ------------------------------------------------------------------------%
% -- Betas Coherent ------------------------------------------------------%
% s3 = subplot(2,3,3);
% yy_beta = [0 0];
% text_disp = 0.25;
% for ii = 1:3
%     
%     errorbar(ii,BETAS_G.diff.stats.mean(ii),BETAS_G.diff.stats.sem(ii),...
%         'Marker','none','Color',clrMap(ii,:),'LineWidth',1.5)
%     
%     hold on
%     
%     b = bar(ii,BETAS_G.diff.stats.mean(ii),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
%     
%     if BETAS_G.diff.stats.mean(ii) < 0 ; text_disp = text_disp * -1; end
%     
%     text(ii,BETAS_G.diff.stats.mean(ii)+text_disp,num2str(BETAS_G.diff.stats.mean(ii)),'HorizontalAlignment','center',...
%         'FontSize',12)
%     
%     hold on
%     
%     if BETAS_G.diff.stats.mean(ii) > yy_beta(2)
%         yy_beta(2) = BETAS_G.diff.stats.mean(ii);
%     end
%     
%     if BETAS_G.diff.stats.mean(ii) < yy_beta(1)
%         yy_beta(1) = BETAS_G.diff.stats.mean(ii);
%     end
%     
% end
% hold off
% 
% xticks(1:3)
% xticklabels(BETAS_G.diff.labels(1:3))
% ylim([floor(yy_beta(1))-1 ceil(yy_beta(2))+1])
% ylabel('beta value')
% title('Betas','FontSize',12)
% s3.FontSize = 14;
% ------------------------------------------------------------------------%
% -- ERA InCoherent ------------------------------------------------------%
s4 = subplot(2,3,4);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.InCoh_aCoh.stats.mean,Y_AUC_G.InCoh_aCoh.stats.sem,'Color',clrMap(4,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.InCoh_aInCoh.stats.mean,Y_AUC_G.InCoh_aInCoh.stats.sem,'Color',clrMap(5,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.InCoh_aNA.stats.mean,Y_AUC_G.InCoh_aNA.stats.sem,'Color',clrMap(6,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(4:6),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal variation (%)'},'FontSize',12)
s4.FontSize = 14;
% ------------------------------------------------------------------------%
% -- AUC InCoherent ------------------------------------------------------%
s5 = subplot(2,3,5);
yy_auc = 0;
for ii = 4:6
    
    errorbar(ii,AUC_G.(trialTestTps{ii}).stats.mean,AUC_G.(trialTestTps{ii}).stats.sem,...
        'Marker','none','Color',clrMap(ii,:),'LineWidth',1.5)
    
    hold on
    
    b = bar(ii,AUC_G.(trialTestTps{ii}).stats.mean,'FaceColor',clrMap(ii,:),'BarWidth',0.4);
    
    text(ii,AUC_G.(trialTestTps{ii}).stats.mean+AUC_G.(trialTestTps{ii}).stats.sem+0.25,num2str(AUC_G.(trialTestTps{ii}).stats.mean),...
        'HorizontalAlignment','center','FontSize',12)
    
    hold on
    
    if AUC_G.(trialTestTps{ii}).stats.mean > yy_auc
        yy_auc = AUC_G.(trialTestTps{ii}).stats.mean;
    end
    
end
hold off

xticks(4:6)
xticklabels(trialLabels(4:6))
ylim([0 ceil(yy_auc)+1])
ylabel('Area under curve')
s5.FontSize = 14;
% ------------------------------------------------------------------------%
% -- Betas InCoherent ----------------------------------------------------%
% s6 = subplot(2,3,6);
% yy_beta = [0 0];
% text_disp = 0.25;
% for ii = 4:6
%     
%     errorbar(ii,BETAS_G.diff.stats.mean(ii),BETAS_G.diff.stats.sem(ii),...
%         'Marker','none','Color',clrMap(ii,:),'LineWidth',1.5)
%     
%     hold on
%     
%     b = bar(ii,BETAS_G.diff.stats.mean(ii),'FaceColor',clrMap(ii,:),'BarWidth',0.4);
%     
%     if BETAS_G.diff.stats.mean(ii) < 0; text_disp = text_disp * -1; end
%     
%     text(ii,BETAS_G.diff.stats.mean(ii)+text_disp,num2str(BETAS_G.diff.stats.mean(ii)),'HorizontalAlignment','center',...
%         'FontSize',12)
%     
%     hold on
%     
%     if BETAS_G.diff.stats.mean(ii) > yy_beta(2)
%         yy_beta(2) = BETAS_G.diff.stats.mean(ii);
%     end
%     
%     if BETAS_G.diff.stats.mean(ii) < yy_beta(1)
%         yy_beta(1) = BETAS_G.diff.stats.mean(ii);
%     end
%     
% end
% hold off
% 
% xticks(4:6)
% xticklabels(BETAS_G.diff.labels(4:6))
% ylim([floor(yy_beta(1))-1 ceil(yy_beta(2))+1])
% ylabel('beta value')
% title('Betas','FontSize',12)
% s6.FontSize = 14;
% ------------------------------------------------------------------------%

%% Export The ultimate Figure
% export_fig(fullfile(outputFolder,sprintf('Group_N%i_UltimateFigure_BilateralMT_D',N,delay)),'-eps','-pdf','-transparent');
print(fig_ultimate,fullfile(outputFolder,sprintf('Group_N%i_UltimateFigure_BilateralMT_D%i',N,delay)),'-dpng')

%% Export workspace
save(fullfile(ioFolder,sprintf('Group_N%i_BilateralMT.mat',N)),'AUC_G','ERA_G','Y_AUC_G','N')


