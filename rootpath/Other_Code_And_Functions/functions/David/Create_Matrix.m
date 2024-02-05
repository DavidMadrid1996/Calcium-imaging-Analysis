function [BloqueWave] = Create_Matrix(Wave)
%CREATE_MATRIX Summary of this function goes here
%   Detailed explanation goes here
S1=readtable(Wave);

BloqueWave=S1(4:end,:)

BloqueWave.Properties.VariableNames={'Time' 'DM4' 'DM2' 'DLD' };

end

