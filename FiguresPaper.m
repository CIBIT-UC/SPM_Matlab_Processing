% clear,clc,close all

% RUN AFTER GROUP ANALYSIS

%% Figure 1
fig1 = figure('Name','Group ERA','Position',[100 100 1100 1000]);
xvector = -2:nTrialVols-3;
clrMap = lines;
trialLabels = {'Coherent \rightarrow Coherent';'InCoherent \rightarrow Coherent';'NonAdapting \rightarrow Coherent';'Coherent \rightarrow InCoherent';'InCoherent \rightarrow InCoherent';'NonAdapting \rightarrow InCoherent'};
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
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.4,{'Test with','Coherent motion'},'HorizontalAlignment','center','FontSize',12)
% text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
% text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([1 2 3]))
s1.FontSize = 12;
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
s2=subplot(2,1,2);
absMax = 0;
for tt = [4 5 6]
    
    errorbar(xvector,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.mean,...
        ERA_G.BilateralMT.(trialTestTps{tt}).stats.sem,...
        '.-','LineWidth',1.5,'MarkerSize',15,'Color',clrMap(tt-3,:));
    
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
ylabel({'BOLD signal change (%)'})

xticks([1 6 12 13 18 19 21 22 24])

line(xx,zeros(length(xx),1),'LineStyle',':','Color','k')

line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([12.5 12.5],yy,'LineStyle','--','Color','k')
line([18.5 18.5],yy,'LineStyle','--','Color','k')
line([21.5 21.5],yy,'LineStyle','--','Color','k')
line([24.5 24.5],yy,'LineStyle','--','Color','k')

text(6,yy(2)-0.4,'First motion period','HorizontalAlignment','center','FontSize',12)
text(15.5,yy(2)-0.4,{'Test with','InCoherent motion'},'HorizontalAlignment','center','FontSize',12)
% text(20,yy(2)-0.2,'MAE','HorizontalAlignment','center','FontSize',12)
% text(23,yy(2)-0.2,'Report','HorizontalAlignment','center','FontSize',12)

legend(trialLabels([4 5 6]))
s2.FontSize = 12;
% ------------------------------------------------------------------------%

%% Calculate point to point t-tests
% Bonferroni correction with 3 comparions

TTestRes = zeros(4,8);

for jj = 1:8
    
   [~,TTestRes(2,jj)] = ttest(Y_AUC_G.Coh_aCoh.data(:,jj),Y_AUC_G.Coh_aInCoh.data(:,jj));
   
   [~,TTestRes(4,jj)] = ttest(Y_AUC_G.InCoh_aCoh.data(:,jj),Y_AUC_G.InCoh_aInCoh.data(:,jj));
   
   TTestRes(2,jj) = TTestRes(2,jj) * 3;
   TTestRes(1,jj) = TTestRes(2,jj) < 0.05;
   TTestRes(4,jj) = TTestRes(4,jj) * 3;
   TTestRes(3,jj) = TTestRes(4,jj) < 0.05;
      
end

%% Figure 2
figure('position',[150 150 1200 500])
clrMap = lines;
x_auc = 15:22; % test block + 1 before + 1 after
xvector = 0:length(x_auc)-1;
xx = [-1 length(x_auc)];
yy = [-0.2 1.2];
lwidth = 1.5;

s1 = subplot(1,2,1);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.Coh_aCoh.stats.mean,Y_AUC_G.Coh_aCoh.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.Coh_aInCoh.stats.mean,Y_AUC_G.Coh_aInCoh.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.Coh_aNA.stats.mean,Y_AUC_G.Coh_aNA.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold on

plot(find(TTestRes(1,:))-1,0.88*ones(1,length(find(TTestRes(1,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(1:3),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s1.FontSize = 12;

%-------------------------------------------------------
s4 = subplot(1,2,2);
line([0.5 0.5],yy,'LineStyle','--','Color','k')
line([6.5 6.5],yy,'LineStyle','--','Color','k')
line(xx,[0 0],'LineStyle',':','Color','k')

hold on

e1 = errorbar(xvector,Y_AUC_G.InCoh_aCoh.stats.mean,Y_AUC_G.InCoh_aCoh.stats.sem,'Color',clrMap(1,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e2 = errorbar(xvector,Y_AUC_G.InCoh_aInCoh.stats.mean,Y_AUC_G.InCoh_aInCoh.stats.sem,'Color',clrMap(2,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);
e3 = errorbar(xvector,Y_AUC_G.InCoh_aNA.stats.mean,Y_AUC_G.InCoh_aNA.stats.sem,'Color',clrMap(3,:),'LineWidth',lwidth,'Marker','.','MarkerSize',25);

hold on

plot(find(TTestRes(3,:))-1,0.88*ones(1,length(find(TTestRes(3,:)))),'*','Color','k','LineWidth',1,'MarkerSize',8)

hold off

xlim(xx), ylim(yy)
xticks(xvector)
legend([e1 e2 e3],trialLabels(4:6),'FontSize',12)
box on
xlabel('Time (volumes)','FontSize',12)
ylabel({'BOLD signal change (%)'},'FontSize',12)
s4.FontSize = 12;
