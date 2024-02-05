function [MaxiMini] = MAXIMINVSDI(Waves)
%MAXIMINVSDI Summary of this function goes here
%   Detailed explanation goes here
[row colum]=size(Waves);
for i=1:colum
    MaxiMini(1,i)=max(Waves(:,i));
    MaxiMini(2,i)=min(Waves(:,i));
end 
end

