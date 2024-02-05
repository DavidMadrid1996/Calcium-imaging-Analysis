function [SlopeR] = RecorteSlope2(Slope,R)
%RECORTESLOPE Summary of this function goes here
%   R, VA EN FRAMES
[row colum]=size(Slope);
for ii=1:row
   SlopeR(ii,:) = Slope(ii,R);
end 
end