function [MaxiEachColumn,GeneralMaxIndex,GeneralMaxi] = MaxiMini(BloqueWavesAM)
%MAXIMINI Summary of this function goes here
%   Detailed explanation goes here
[row column] =size (BloqueWavesAM)
for i=2:column
MaxiEachColumn(i-1)=max(max(BloqueWavesAM(:,i)));
GeneralMaxIndex=find(MaxiEachColumn==max(MaxiEachColumn));
GeneralMaxi=find(BloqueWavesAM(:,GeneralMaxIndex+1)==max(MaxiEachColumn))

end 
end 

