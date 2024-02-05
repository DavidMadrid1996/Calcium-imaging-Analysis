
%% Importing brain data
close
brain = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Wholebrain');
chx10 = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Chx10');
mauthnerR = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Mauthner Right');
mauthnerL = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Mauthner Left');
backfill = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Backfill');
% double = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx");

colors = [0.3010 0.7450 0.9330]; 
% Bilateral --> [0.3010 0.7450 0.9330]
% left --> [0.8500 0.3250 0.0980]

%% Selecting type of plot

prompt = 'Do you want to plot in 2D or 3D ->  Select 2 or 3 : ';
str = str2double(input(prompt,'s'));

if str==2 % This will be 2D ploting

    figure(1)

    % Ploting scafold
    subplot(6,6,[4 5 10 11 16  17 22 23]) %subplot(5,5,[4  9  14  19  24 25])

    % Ploting brain
    % scatter(brain(:,1),brain(:,2),'o','MarkerEdgeAlpha','0.05')

    % Holding
    hold on

    % Ploting chx10 population
    scatter(chx10(:,1),chx10(:,2),15,[168./255, 168./255, 168./255],'filled')

    % Ploting mauthner cells
%     scatter(mauthnerR(:,1),mauthnerR(:,2),'b')
%     scatter(mauthnerL(:,1),mauthnerL(:,2),'b')

    % Rescaling
%     ylim ([-0.5 3.5])
   ylim ([-0.5 2.5])




S_11_Group_analysis_files_2


% binsize = 5; % for the histogram of probability
binwidth = 0.1; 

        %% Ploting all cells
% 
%         for i = 1:length(TOTAL_COORDINATES)
%             figure(1)
%             subplot(6,6,[4 5 10 11 16  17 22 23])
%             scatter(TOTAL_COORDINATES((i),1),TOTAL_COORDINATES((i),2),40,	[0.4940, 0.1840, 0.5560],'filled'), hold on % Coordinates of the lateral refere
% 
% %             subplot(6,6,[1 2 3])
% %             plot(-CALCIUMroiTS.diff_perc03.data((i),:),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2), hold on
% %             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
%         end
% %         subplot(6,6,[1 2 3])
% %         title('Al cells')
%         if isempty(i)
%         else
%             % Histogram for all the neurons vertical Axis
%             % For All neurons
%             subplot(6,6,[6 12 18 24 ])
%             Hist(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),2),-0.5,2.5, 20,[0.4940, 0.1840, 0.5560],'V'),hold on
%             figure(2),hold on
% %             h1 =histogram(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),2), 'Normalization','probability', 'DisplayStyle','bar'); 
% %             h1.FaceColor = [0.4940, 0.1840, 0.5560]; 
% %             h1.BinWidth = binwidth ;  
%             xlim([-0.5,2.5]);
% 
%             [f,xi,bw] = ksdensity(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),2));
%             hold on
%             plot(xi,f./10,'color',[0.4940, 0.1840, 0.5560])
% %             xlim([0 100])
% 
%             % Histogram for all the neurons horizontal Axis
%             figure(1)
%             subplot(6,6,[ 28 29 34 35])
%             Hist(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),1),-1.5,1.5, 20,[0.4940, 0.1840, 0.5560],'H'),hold on
%             figure(3),hold on
% %             h2 = histogram(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),1),  'Normalization','probability', 'DisplayStyle','bar'); 
% %             h2.FaceColor = h1.FaceColor; 
% %             h2.BinWidth = binwidth ;  
%             xlim([-1.5,1.5]);
% 
%             [f,xi,bw] = ksdensity(TOTAL_COORDINATES(1:length(TOTAL_COORDINATES),1));
%             hold on
%             plot(xi,f./10,'color',[0.4940, 0.1840, 0.5560])
% %             ylim([0 100])
% 
%         end






  

%         %% Ploting reduction of activity
% 
%         for i = 1:length(ReductionActivity_COORDINATES)
%             figure(1)
%             subplot(6,6,[4 5 10 11 16  17 22 23 ])
%              scatter(ReductionActivity_COORDINATES((i),1),ReductionActivity_COORDINATES((i),2),  40,'r','filled'), hold on % Coordinates of the lateral refere
% 
% %             subplot(6,6,[19 20 21])
% %             plot(-CALCIUMroiTS.diff_perc03.data(ReductionActivity.(myfield1)(i),:),'r','LineWidth',2),hold on
% %             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
%         end
% %         subplot(6,6,[19 20 21])
% %         title('Reduction Activity')
% 
%         if isempty(i)
%         else
%             subplot(6,6,[6 12 18 24 ])
%             Hist(ReductionActivity_COORDINATES(:,2),0,3.5, 20,'r','V')
% 
%             subplot(6,6,[ 28 29 34 35])
%             Hist(ReductionActivity_COORDINATES(:,1),-2.5,2.5, 20,'r','H')
%         end


%         %% Ploting pre activity
% 
%         for i = 1:length(PreActivity_COORDINATES)
%             figure(1)
% 
%             subplot(6,6,[4 5 10 11 16  17 22 23  ])
% 
%              scatter(PreActivity_COORDINATES((i),1),PreActivity_COORDINATES((i),2),40  ,'m','filled'), hold on % Coordinates of the lateral refere
% 
% %             subplot(6,6,[25 26 27])
% %             plot(-CALCIUMroiTS.diff_perc03.data(PreActivity.(myfield1)(i),:),'m','LineWidth',2),hold on
% %             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
%         end
% %         subplot(6,6,[25 26 27])
% %         title('Pre Activity')
% 
% 
%         if isempty(i)
%         else
%             subplot(6,6,[6 12 18 24 ])
%             Hist(PreActivity_COORDINATES(:,2),0,3.5, 20,'m','V')
% 
%             subplot(6,6,[ 28 29 34 35])
%             Hist(PreActivity_COORDINATES(:,1),-2.5,2.5, 20,'m','H')
%         end

        %% Ploting non active

        for i = 1:size(NonActive_COORDINATES,1)
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23  ])
             scatter(NonActive_COORDINATES((i),1),NonActive_COORDINATES((i),2),40  ,'k','filled'), hold on % Coordinates of the lateral refere

%             subplot(6,6,[31 32 33])
%             plot(-CALCIUMroiTS.diff_perc03.data(NonActive.(myfield1)(i),:),'k','LineWidth',2),hold on
%             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
%         subplot(6,6,[31 32 33])
%         title('Non Active')


        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(NonActive_COORDINATES(:,2),-0.5,2.5, 20,'k','V')
            figure(2)
%             h3 = histogram(NonActive_COORDINATES(:,2),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h3.FaceColor = 'k'; 
%            h3.BinWidth = binwidth ; 
           [f,xi,bw] = ksdensity(NonActive_COORDINATES(:,2));
            hold on
            plot(xi,f./10,'color','k')



            figure(1)
            subplot(6,6,[ 28 29 34 35])
            Hist(NonActive_COORDINATES(:,1),-1.5,1.5, 20,'k','H')
            figure(3)
%             h4 = histogram(NonActive_COORDINATES(:,1),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h4.FaceColor = h3.FaceColor; 
%             h4.BinWidth = binwidth ; 
              [f,xi,bw] = ksdensity(NonActive_COORDINATES(:,1));
            hold on
            plot(xi,f./10,'color','k')
        end


        %%         %% Ploting long delay


        for i = 1:size(LongDelay_COORDINATES,1)
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23 ])
             scatter(LongDelay_COORDINATES((i),1),LongDelay_COORDINATES((i),2),40,[0.9290, 0.6940, 0.1250],'filled'), hold on % Coordinates of the lateral refere
%             subplot(6,6,[7 8 9])
pause(1)
%             plot(-CALCIUMroiTS.diff_perc03.data(LongDelay.(myfield1)(i),:),'color',[0.9290, 0.6940, 0.1250],'LineWidth',2), hold on
%             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
%         subplot(6,6,[7 8 9])
%         title('Long delay')

        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(LongDelay_COORDINATES(:,2),-0.5,2.5, 20,[0.9290, 0.6940, 0.1250],'V')
            figure(2)
%             h5 = histogram(LongDelay_COORDINATES(:,2),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h5.FaceColor = [0.9290, 0.6940, 0.1250]; 
%             h5.BinWidth = binwidth ; 
          [f,xi,bw] = ksdensity(LongDelay_COORDINATES(:,2));
            hold on
            plot(xi,f./10,'color',[0.9290, 0.6940, 0.1250])



            figure(1)
            subplot(6,6,[ 28 29 34 35])
            Hist(LongDelay_COORDINATES(:,1),-1.5,1.5, 20,[0.9290, 0.6940, 0.1250],'H')
            figure(3)
%             h6 = histogram(LongDelay_COORDINATES(:,1),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h6.FaceColor = h5.FaceColor;
%             h6.BinWidth = binwidth ; 
                [f,xi,bw] = ksdensity(LongDelay_COORDINATES(:,1));
            hold on
            plot(xi,f./10,'color',[0.9290, 0.6940, 0.1250])

        end


              %% Ploting no/short delay

        for i = 1:size(ShortDelay_COORDINATES,1)
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23 ])
            scatter(ShortDelay_COORDINATES((i),1),ShortDelay_COORDINATES((i),2),  40,colors,'filled'), hold on % Coordinates of the lateral refere

%             subplot(6,6,[13 14 15])
%             plot(-CALCIUMroiTS.diff_perc03.data(ShortDelay.(myfield1)(i),:),'g','LineWidth',2),hold on
%             xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
%         subplot(6,6,[13 14 15])
%         title('No / Short delay')

        if isempty(i)
        else

            subplot(6,6,[6 12 18 24 ])
            Hist(ShortDelay_COORDINATES(:,2),-0.5,2.5, 20,colors,'V')
            figure(2)
%             h7 = histogram(ShortDelay_COORDINATES(:,2),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h7.FaceColor = colors; 
%             h7.BinWidth = binwidth ; 
           [f,xi,bw] = ksdensity(ShortDelay_COORDINATES(:,2));
            hold on
            plot(xi,f./10,'color',colors)


            figure(1)
            subplot(6,6,[ 28 29 34 35])
            Hist(ShortDelay_COORDINATES(:,1),-1.5,1.5, 20,colors,'H')
            figure(3)
%             h8 = histogram(ShortDelay_COORDINATES(:,1),  'Normalization','probability', 'DisplayStyle','bar'); 
%             h8.FaceColor = h7.FaceColor;
%             h8.BinWidth = binwidth ; 
             [f,xi,bw] = ksdensity(ShortDelay_COORDINATES(:,1));
            hold on
            plot(xi,f./10,'color',colors)
        end


        %% Creating folder to saved the plot
%         if exist('CATEGORY', 'dir')
%         else
%             mkdir (fullfile(path.data,'dataCALCIUM',CALCIUM.ref,'CATEGORY')); % Creating the directory for dataCALCIUM
%         end
%         saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),'CATEGORY',strcat('CATEGORIZED_',char(CALCIUM.list(nfish,1)),'.fig')))
%         saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),'CATEGORY',strcat('CATEGORIZED_',char(CALCIUM.list(nfish,1)),'.eps')))
%         uiwait
end 