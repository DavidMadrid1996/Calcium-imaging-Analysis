function Subplot_ERPs(waves,timebase, trials2plot,  roi2plot, labels, titletxt, stim_dur,general_title)
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
figure 
for roi = 1:nroi


sgtitle(general_title) ; 

select2plot = squeeze(waves(:,roi,trials2plot));

% plot ALL TRIALS
subplot(2,3,roi)
plot(timebase, select2plot); hold on

% surf(select2plot)

title(strcat(titletxt,'.',labels{roi})) ; 
%plot mean
 plot(timebase, mean(select2plot,2, 'omitnan'), 'Linewidth', 3, 'Color', '#292924')

% plot settings
 ylim([min(min(min(waves))) max(max(max(waves)))])
xline(0, '--'); xline(stim_dur, '--'); 
xlabel('Time')
ylabel('% \Delta F') % note that matlab interprets "\mu" as the Greek character for micro


end

end
