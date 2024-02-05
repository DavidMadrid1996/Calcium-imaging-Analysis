
%% Quick loading of the animal. Copy the path to the files

load("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\VALID_RECONSTRUCTION_Final_Script_to_score_cells\GCaMP_221221_2_BAD_VISIBILITY_TAIL\Calcium_ourToolbox\Calcium_data\dataCALCIUM\221221_2_10\CALCIUM_221221_2_10.mat");
load("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\VALID_RECONSTRUCTION_Final_Script_to_score_cells\GCaMP_221221_2_BAD_VISIBILITY_TAIL\Calcium_ourToolbox\Calcium_data\datawaves\221221_2_10\CALCIUMRoiTS_221221_2_10.mat");
[CALCIUMroiTS] =VSDroiTS; 
sr=1/CALCIUMroiTS.deeplabcut.sr;
wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)

%% Quick check of the swimming
close 
figure(1)
yyaxis left
plot( CALCIUMroiTS.deeplabcut.time_scaled,wave),hold on
yline(0)
for ii = 1:length(fieldnames(CALCIUM.EPISODE_X))
    temp1 = fieldnames(CALCIUM.EPISODE_X);
    myfield1 = char(temp1(ii)); % extracting the fieldnames of each cell of the loop
    plot(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),'LineWidth',2)
    scatter(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),30,'filled')
         text(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),string(CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.( myfield1)))
end

          % Correction problem with the trasposition
        if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
            CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
%             CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
        else
        end
        yyaxis right
plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(15,:),'g','LineWidth',2)

        % Examples of traces that you want to visualize

        close all
        temp = -CALCIUMroiTS.diff_perc03.data(CALCIUM.ShortDelay.swm_1,:); color='b'; 
         clear tempfilter

        temp = - CALCIUMroiTS.diff_perc03.data(CALCIUM.LongDelay.swm_1,:);color='g'; 
         clear tempfilter

        temp = - CALCIUMroiTS.diff_perc03.data(CALCIUM.NonActive.swm_1([1 5 2 ]),:);color='k'; 
        clear tempfilter


        prompt = 'Do you want to normalize it?';
        str = input(prompt,'s');
        for i = 1:size(temp,1)
            if str=='Y'
                if max(temp(i,:))<0

                    tempfilter(i,:) = medfilt1(temp(i,:)./min(temp(i,:)),8);
                else
                    tempfilter(i,:) = medfilt1(temp(i,:)./max(temp(i,:)),8);
                end
            elseif str=='N'
                tempfilter(i,:) = medfilt1(temp(i,:))
            end
            plot(CALCIUM.timebase, tempfilter(i,:),color),hold on
        end
  plot(CALCIUM.timebase,mean(tempfilter,1),color,'LineWidth',3)



      

        writematrix(right,"C:\Users\Usuario\OneDrive\Escritorio\Calcium imaging graph.xlsx",'Sheet','Calcium Traces','Range','A2')

   %% Loading movie 
clear movies4D
data = bfopen(char(strcat(CALCIUM.CZI_directory_loop,'.czi')));
% czifinfo(char(strcat(CALCIUM.CZI_directory_loop,'.czi')))

for j = 1:length(data{1,1})
  movies4D (:,:,j) = data{1,1}{j,1};
end 

 movies4D = im2double(movies4D);

 th = 0.01
 figure(1)
  imagesc(mean(movies4D(:,:,[200:230]),3),[0 th]); 
colormap('turbo')
colorbar
   figure(2)
    imagesc(mean(movies4D(:,:,[239:259]),3),[0 th]);
colormap('turbo')
colorbar

