function plot_ERPs(waves,timebase, trials2plot,  roi2plot, labels, titletxt, stim_dur)
% plot_ERPs(waves,timebase, trials2plot,  roi2plot, labels, titletxt, stim_dur)

% INPUT
% waves = waves;
% labels = VSDroiTS.roi.labels;
% trials2plot = VSDI.nonanidx;
% roi2plot = [1,5];
% timebase = VSDI.timebase;
% titletxt = 'titulo';
%'stim_dur' : in ms

nroi = length(roi2plot);

for roi = 1:nroi

figure 
sgtitle(strcat(titletxt,'.',labels{roi})) ; 

select2plot = squeeze(waves(:,roi,trials2plot));

% plot ALL TRIALS
plot(timebase, select2plot); hold on
%plot mean
plot(timebase, mean(select2plot,2, 'omitnan'), 'Linewidth', 3, 'Color', '#292924')

% plot settings
% ylim([ymin ymax])
xline(0, '--'); xline(stim_dur, '--'); 
xlabel('Time')
ylabel('% \Delta F') % note that matlab interprets "\mu" as the Greek character for micro

end

end