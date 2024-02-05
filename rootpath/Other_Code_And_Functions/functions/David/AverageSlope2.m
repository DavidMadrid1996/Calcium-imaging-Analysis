function [AverageSlopeR] = AverageSlope2(SlopeR)
%AVERAGESLOPE Summary of this function goes here
%   Detailed explanation goes here
[row colum]=size(SlopeR);
for ii=1:row
    AverageSlopeR(ii)=mean(SlopeR(ii,:));
end 
end
