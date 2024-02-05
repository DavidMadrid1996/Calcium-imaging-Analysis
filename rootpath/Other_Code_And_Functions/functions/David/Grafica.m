function [outputArg1,outputArg2] = Grafica(inputArg1,inputArg2,TimeVectorM,nombre,onset,maxi,mini,RutaGuardadoGraficas)
%PLOTING Summary of this function goes here
%   Detailed explanation goes here
figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual');
hold on
grid on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
plot(TimeVectorM,inputArg1,'LineWidth',3);
hold on
plot(TimeVectorM,inputArg2,'LineWidth',3);
title(nombre,'FontSize',20);
xlabel('TimeVectorM','FontSize',20);
ylabel('%F','FontSize',20);
legend('Fenil','Sacarosa');
legend('location','southoutside');
legend('Orientation','horizontal');
xline(onset,'LineWidth',3,'DisplayName','Onset');
axis([0 max(TimeVectorM) mini maxi]);
fullFileName = fullfile(RutaGuardadoGraficas, nombre);

saveas(gcf, fullFileName);
end

