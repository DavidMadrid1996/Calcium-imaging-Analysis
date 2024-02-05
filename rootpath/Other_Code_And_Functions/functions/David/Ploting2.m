function  Plotting2(AverageSlopeR,BloqueWave,nombre,RutaGuardadoGraficas)
%PLOTTING Summary of this function goes here
%   Detailed explanation goes here
figure
 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for ii=1:length(AverageSlopeR)
    y=AverageSlopeR(ii)*[0:0.001:1];
    grid on
    plot(y,'r-','LineWidth',2)
    hold on
end 
lgdn=string(lgd)
legend(lgdn,'FontSize',25)
legend show
legend('location','northeast')
legend('Orientation','vertical')

title(nombre,'FontSize',20);
xlabel('Time','FontSize',20);
ylabel('%F','FontSize',20);
fullFileName = fullfile(RutaGuardadoGraficas, strcat(nombre, '.jpg'));
saveas(gcf, fullFileName);
end
