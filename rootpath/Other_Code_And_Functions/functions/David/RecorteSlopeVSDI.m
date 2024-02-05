function  [SlopeR] = RecorteSlopeVSDI(Slope,R)
%RECORTESLOPEVSDI Summary of this function goes here
%   Detailed explanation goes here
[row colum]=size(Slope);
for ii=1:colum
   SlopeR(1:length(R),ii) = Slope(R,ii);
end 
end

