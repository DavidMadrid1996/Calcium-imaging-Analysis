function [TimeVector] = Time_Vector(Sabor1)
%TIME_VECTOR Summary of this function goes here
%   Detailed explanation goes here
S1=readtable(Sabor1);
TimeVector=S1(4:end,1);
TimeVector.Properties.VariableNames= {'Tiempos'};
end

