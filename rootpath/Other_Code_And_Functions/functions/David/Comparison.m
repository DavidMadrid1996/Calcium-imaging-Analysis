function [outputArg1,outputArg2] = Comparison(BloqueWaveM,BloqueWave,titulo,RutaGuardadoGraficas)
%COMPARISON Summary of this function goes here
%   Detailed explanation goes here
[rows columns]=size(BloqueWave);
figure
for ii=2:columns
    grid on
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    plot(BloqueWaveM(:,1),BloqueWaveM(:,ii),'LineWidth',2)
    lgd{ii-1}=BloqueWave.Properties.VariableNames{ii}
    hold on
    

    % xlabel(Parametro,'FontSize',20);
    % ylabel('Amplitud','FontSize',20);
end
xline(300,'LineWidth',3)
lgdn=string(lgd)
legend(lgdn,'FontSize',25)
legend show
legend('location','northeast')
legend('Orientation','vertical')
title(titulo,'FontSize',20);
xlabel('Time','FontSize',20);
ylabel('%F','FontSize',20);
fullFileName = fullfile(RutaGuardadoGraficas, strcat(titulo, '.jpg'));
saveas(gcf, fullFileName);

end

