function [Slope] = Derivative(BloqueWaves)
%DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
[row colum]=size(BloqueWaves);
for ii=2:colum
    for iii=1:row-1
    Slope(iii,ii-1)=(BloqueWaves(iii+1,ii)-BloqueWaves(iii,ii))/(BloqueWaves(2,1)-BloqueWaves(1,1));
    end 
end 
end

