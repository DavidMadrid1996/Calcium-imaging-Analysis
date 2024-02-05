function  DotPlotting(AverageSlopeR,titulo,VSDI,parametro)
%DOTPLOTTING Summary of this function goes here
%   Detailed explanation goes here
c={'r','b','k','g','y','c', 'm','r','b','k','g','y','c'}
P=[1 2 3 4 5 6 7 8 9 10 11]
[row colum]=size(AverageSlopeR);
for ii=1:colum
 plot(P(ii),AverageSlopeR(:,ii),strcat(c{ii},'.'),'MarkerSize',30)
hold on
end 
title(titulo)
ylabel(parametro,'FontSize',20)
names =VSDI.roi.labels%{'dm4L_L'; 'dm4M_L'; 'dm4L_R'; 'dm4M_R'; 'dm3_L';'dm3_R';'dm2_L';'dm2_R';'dld_L'...
%     ;'dm1_L';'dm1_r'};
set(gca,'xtick',[1:11],'xticklabel',names,'FontSize',20)

end

