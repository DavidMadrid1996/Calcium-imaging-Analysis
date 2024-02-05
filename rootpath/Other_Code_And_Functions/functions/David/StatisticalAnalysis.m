function [sTable] = StatisticalAnalysis(B,VSDI)
%STATISTICALANALYSIS Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(B,2)
   for j=1:size(B,2)
        [h,p] = ttest2(B(:,i),B(:,j),'Alpha',0.00001);
        Decisiondm4M_L=h;
        ValorSdm4M_L=p;
        
        SignificantMatrix(i,j,1)=Decisiondm4M_L;
        SignificantMatrix(i,j,2)=ValorSdm4M_L;
   end 
end 
rowNames = VSDI.roi.labels;
colNames = VSDI.roi.labels;
sTable = array2table(SignificantMatrix(:,:,1),'RowNames',rowNames,'VariableNames',colNames)
figure
imagesc(SignificantMatrix(:,:,1))
names =VSDI.roi.labels
set(gca,'xtick',[1:11],'xticklabel',names,'FontSize',20)
set(gca,'xaxisLocation','top')
set(gca,'ytick',[1:11],'yticklabel',names,'FontSize',20)
% colormap(jet)
% sgtitle('Significancia entre ROIs')
end

