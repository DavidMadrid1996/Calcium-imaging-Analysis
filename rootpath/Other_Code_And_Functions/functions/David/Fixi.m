function [BloqueWavesAMF] = Fixi(BloqueWavesAM)
%FIXI Summary of this function goes here
%   Detailed explanation goes here
PorMil=BloqueWavesAM.*10000;
Fijado=fix(PorMil);
BloqueWavesAMF=Fijado./10000;
end

