%{

The lentgh of the behaviour video with ligth must be the same than the duration of the confocal video. 
If the duration are not the same there was some problem with the aqusition
and can comprimise details analysis as frecuency of the swimming and
syncronization of the activity. 

--> If we consider the duration of the confocal video the correct one and
the duration of the behaviour video is shorter then we can enlarge the
behaviour video  to match the duration of the confocal video

%}

S_1_Loading_Animal

% First we calculated already the delay btw the confocal video and the
% behavioural video and is store in the variable CALCIUM.delay. Now we know
% the to the last value of the behavioural video must be added CALCIUM.delay 
% to get the last value of the confocal video. For each time of the CALCIUMroiTS.deeplabcut.lasertimesync
% we must add an acumulation of the quantity sr_rescaled defined below so
% at the last value this acumulation reach the value CALCIUM.delay 

sr_rescaled = -CALCIUM.delay/length((CALCIUMroiTS.deeplabcut.lasertimesync));
time_scaled = NaN(1,length(CALCIUMroiTS.deeplabcut.lasertimesync(1))); 
time_scaled(1) = sr_rescaled + CALCIUMroiTS.deeplabcut.lasertimesync(1); 
suma =  sr_rescaled; 


for i = 2:length((CALCIUMroiTS.deeplabcut.lasertimesync))
    suma = suma + sr_rescaled; 
   time_scaled(i) = suma + CALCIUMroiTS.deeplabcut.lasertimesync(i) ; 
end 


             % Correction problem with the trasposition
        if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
            CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
%            CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
        else
        end


        
figure(1)
plot(CALCIUMroiTS.diff_perc03.times,CALCIUMroiTS.diff_perc03.data(1,:)),hold on
%   wave=wave(1:end-1)
plot(time_scaled,wave)

CALCIUMroiTS.deeplabcut.time_scaled = time_scaled'; 
CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);

% Finally we recalculate the frecuency of each swimming episode

for k=1:length(fieldnames(CALCIUM.EPISODE_X)) % we will do it across swimming episodes
    temp1 = fieldnames(CALCIUM.EPISODE_X);
    myfield1 = char(temp1(k));
    temp =  CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1));
    for kk = 1:length(CALCIUM.EPISODE_X.(myfield1))-1
        CALCIUM.FREQ_SWIM_Y_RESCALED.(myfield1)(kk) = (1/(temp(kk+1)-temp(kk)))/2;
    end

end
CALCIUMimg('save',CALCIUM,[],list,nfish);
figure(2)
plot(CALCIUM.FREQ_SWIM_Y_RESCALED.swm_1), hold on
 plot(CALCIUM.FREQ_SWIM_Y.swm_1)





