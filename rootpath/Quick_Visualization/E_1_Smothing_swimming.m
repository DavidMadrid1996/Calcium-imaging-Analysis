%% Smooth Swim Trace
figure
wave_smooth = medfilt1(wave,10);
plot(CALCIUMroiTS.deeplabcut.time_scaled,wave_smooth,'LineWidth',1)