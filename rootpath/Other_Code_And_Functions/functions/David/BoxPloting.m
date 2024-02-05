function  [B,media]= BoxPloting(AvrgROIsWavesI,titulo,VSDI,parametro)
%BOXPLOTING Summary of this function goes here
%   Detailed explanation goes here
row=size( AvrgROIsWavesI,1);
colum=size( AvrgROIsWavesI,2);
threeD=size( AvrgROIsWavesI,3);
P=[1 2 3 4 5 6]
names =VSDI.roi.labels

for i=1:colum
   A =  AvrgROIsWavesI(1,i,:);
   A=squeeze(A);
   B(:,i)=A
end 
boxplot(B,names)
hold on 
media=mean(B)
plot(media,'k-','LineWidth',3)
title(titulo)
ylabel(parametro,'FontSize',20)
set(gca,'xtick',[1:11],'xticklabel',names,'FontSize',20)
end
% ∆F/ms