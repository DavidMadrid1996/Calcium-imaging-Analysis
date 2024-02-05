function [Slope] = Derivative2(Matrix)
%DERIVATIVE2 Summary of this function goes here
%   Detailed explanation goes here
[row colum]=size(Matrix);
for ii=1:row
    for iii=1:colum-1
    Slope(ii,iii)=(Matrix(ii,iii+1)-Matrix(ii,iii))/5;
    end 
end 
end

