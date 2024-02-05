function  PlottingVSDI(AvrgROIsWaves,BloqueWave,nombre
%PLOTTINGVSDI Summary of this function goes here
%   Detailed explanation goes here

 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
 for i=1:size(AvrgROIsWavesI,2)
     Rois=squeeze(AvrgROIsWavesI(:,i,:))
for ii=3:size(AvrgROIsWavesI,3)
    y=Rois(ii)*[0:0.001:1];
    grid on
    plot(y,'k-','LineWidth',2)
    hold on
end 
% lgd{ii}=
% lgdn=string(lgd)
% legend(lgdn,'FontSize',25)
% legend show
% legend('location','northeast')
% legend('Orientation','vertical')

% title(nombre,'FontSize',20);
xlabel('Time','FontSize',20);
ylabel('%F','FontSize',20);
 end 
% fullFileName = fullfile(RutaGuardadoGraficas, strcat(nombre, '.jpg'));
% saveas(gcf, fullFileName);
end
