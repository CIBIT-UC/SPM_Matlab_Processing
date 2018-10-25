prtFile = xff('C:\Users\alexandresayal\DATA_AdaptationCheck_SPM\Calibtest\SC_RunC1.prt');

names = prtFile.ConditionNames';
onsets = cell(1,length(names));
durations = cell(1,length(names));
onoff = prtFile.OnOffsets;

for c = 1:length(names)
    onsets{c} = onoff(onoff(:,1)==c,2)';
    durations{c} = (onoff(onoff(:,1)==c,3)-onoff(onoff(:,1)==c,2)+1)';
end

save('C:\Users\alexandresayal\DATA_AdaptationCheck_SPM\Calibtest\SC_RunC1.mat','names','onsets','durations')