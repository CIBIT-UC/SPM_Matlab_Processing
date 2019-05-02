subjectName = 'VPIS08';

% Read SPM data
load([subjectName '-MT-5mm-SPM.mat'])

% Read BV
vtcFile = xff(fullfile('E:','DATA_VP_Inhibition',subjectName,'run-run1-data','PROJECT','ANALYSIS',...
    'VPIS08_run1_SCCTBL_3DMCTS_LTR_THP0.01Hz_TAL_SD3DVSS5.00mm.vtc'));
voiFile = xff(fullfile('E:','DATA_VP_Inhibition',subjectName,'run-run1-data','PROJECT','ANALYSIS',...
    'VOIs-MT-5mm.voi'));

TCs = vtcFile.VOITimeCourseOrig(voiFile,Inf);

%% Plots 
figure,
subplot(2,2,1)
plot(mean(TCs{1},2) / mean(mean(TCs{1})) * 100 )
title('BV - leftMT')
subplot(2,2,2)
plot(mean(TCs{2},2) / mean(mean(TCs{2})) * 100 )
title('BV - rightMT')
subplot(2,2,3)
plot(mean(leftMT_5mm_tc,2))
title('SPM - leftMT')
subplot(2,2,4)
plot(mean(rightMT_5mm_tc,2))
title('SPM - rightMT')

%%
figure('position',[150 150 1200 600])
yyaxis left
plot(mean(TCs{1},2), 'linewidth',1.5)
ylabel('BV')
yyaxis right
plot(mean(leftMT_5mm_tc,2), 'linewidth',1.5)
ylabel('SPM')

title('LeftMT')

%%
figure('position',[150 150 1200 600])
yyaxis left
plot(mean(TCs{2},2), 'linewidth',1.5)
ylabel('BV')
yyaxis right
plot(mean(rightMT_5mm_tc,2), 'linewidth',1.5)
ylabel('SPM')

title('RightMT')