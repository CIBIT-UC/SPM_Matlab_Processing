function plotMotionCurves(mcFolder, runName)

   auxDir = dir(fullfile(mcFolder,'rp*.txt'));

   loadmot = load(fullfile(auxDir(1).folder,auxDir(1).name));
   figure;
   subplot(2,1,1);
   plot(loadmot(:,1:3));
   grid on;
%    ylim([-3 3]);  % enable to always scale between fixed values
   title(['Motion parameters: shifts (top, in mm) and rotations (bottom, in dg)'], 'interpreter', 'none');
   
   subplot(2,1,2);
   plot(loadmot(:,4:6)*180/pi);
   grid on;
%     ylim([-3 3]);
   title(['Data from ' runName], 'interpreter', 'none');
%    mydate = date;
   
%    motname = [mcFldr filesep 'motion_sub' sprintf('%02.0f', i) '_' mydate '.png'];
%    print(printfig, '-dpng', '-noui', '-r100', motname);

end

