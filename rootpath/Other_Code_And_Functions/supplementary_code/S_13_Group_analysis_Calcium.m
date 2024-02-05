% %% MAXIMUM CALCIUM RESPONSE IN ALL THE RECORDINGS (FOR ALL THE RECORDINGS OF THE SAME ANIMAL)
% user_settings
% for iiii=1:length(list)
%     nfish =iiii; % Number of the fish to load (this is the trial for the same fish)
%
%     [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
%     CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
%     sr=1/CALCIUMroiTS.deeplabcut.sr;
%
%
%     MAXCALCIUM_RESPONSE_TOTAL(iiii) = CALCIUM.MAXCALCIUM_RESPONSE;
% end
%
%
%
% for iiii=1:length(list)
%     nfish =iiii; % Number of the fish to load (this is the trial for the same fish)
%
%     [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
%     CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
%     sr=1/CALCIUMroiTS.deeplabcut.sr;
%
%     CALCIUM.MAXCALCIUM_RESPONSE_TOTAL = max(MAXCALCIUM_RESPONSE_TOTAL);
%
%
%     CALCIUMimg('save',CALCIUM,[],list,nfish);
% end
%
% clear wave

%% Data extraction. RUN ALWAYS THIS FIRST
S_11_Group_analysis_files

close
brain = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Wholebrain');
chx10 = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Chx10');
mauthnerR = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Mauthner Right');
mauthnerL = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Mauthner Left');
backfill = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx",'Backfill');
% double = xlsread("E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\ImarisMatlab_script_to_import_cells_and_surface_from_IMARIS\Norm brain leander\Brain Chx10 + backfill coords.xlsx");

%% Selecting type of plot

prompt = 'Do you want to plot in 2D or 3D ->  Select 2 or 3 : ';
str = str2double(input(prompt,'s'));




if str==2 % This will be 2D ploting

    figure(1)

    % Ploting scafold
    subplot(6,6,[4 5 10 11 16  17 22 23]) %subplot(5,5,[4  9  14  19  24 25])

    % Ploting brain
    %     scatter(brain(:,1),brain(:,2),'o','MarkerEdgeAlpha','0.05')

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




    for j = 1:size(ValidPathdataCalcium,1)

        load(char(fullfile(ValidPathdataCalcium(j,:),strcat('CALCIUM_',strcat(Animal_List_unique(j,1),...
            Animal_List_unique(j,2),'.mat')))));
        [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(j,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(j,1),...
            Animal_List_unique(j,2),'.mat')))));
        [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;

        LongDelay =  CALCIUM.LongDelay;
        ShortDelay =  CALCIUM.ShortDelay;
        ReductionActivity=  CALCIUM.ReductionActivity;
        PreActivity=  CALCIUM.PreActivity;
        ActiveNeurons = CALCIUM.ActiveNeurons;
        TotalNeurons = CALCIUM.TotalNeurons;
        NonActive = CALCIUM.NonActive;

        COORDINATES =CALCIUM.ZAbsoluteCoordinates;

        myfield1 = char(Animal_List_unique(j,3));

        % Correction problem with the trasposition
        if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
            CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
%             CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
        else
        end






        %% Ploting all cells

        for i = 1:length(COORDINATES)
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23])
            scatter(COORDINATES((i),1),COORDINATES((i),2),40,	[0.4940, 0.1840, 0.5560],'filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[1 2 3])
            plot(-CALCIUMroiTS.diff_perc03.data((i),:),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2), hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[1 2 3])
        title('Al cells')
        if isempty(i)
        else

            % Histogram for all the neurons vertical Axis
            % For All neurons
            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES(1:length(COORDINATES),2),-0.5,3.5, 20,[0.4940, 0.1840, 0.5560],'V'),hold on

            % Histogram for all the neurons horizontal Axis
            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES(1:length(COORDINATES),1),-2.5,2.5, 20,[0.4940, 0.1840, 0.5560],'H'),hold on

        end


        %% Ploting long delay


        for i = 1:length(LongDelay.(myfield1))
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23 ])
            scatter(COORDINATES(LongDelay.(myfield1)(i),1),COORDINATES(LongDelay.(myfield1)(i),2),40,	[0.9290, 0.6940, 0.1250],'filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[7 8 9])
            plot(-CALCIUMroiTS.diff_perc03.data(LongDelay.(myfield1)(i),:),'color',[0.9290, 0.6940, 0.1250],'LineWidth',2), hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[7 8 9])
        title('Long delay')

        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES(LongDelay.(myfield1),2),0,3.5, 20,[0.9290, 0.6940, 0.1250],'V')

            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES(LongDelay.(myfield1),1),-2.5,2.5, 20,[0.9290, 0.6940, 0.1250],'H')
        end



        %% Ploting no/short delay

        for i = 1:length(ShortDelay.(myfield1))
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23 ])
            scatter(COORDINATES(ShortDelay.(myfield1)(i),1),COORDINATES(ShortDelay.(myfield1)(i),2),40,'g','filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[13 14 15])
            plot(-CALCIUMroiTS.diff_perc03.data(ShortDelay.(myfield1)(i),:),'g','LineWidth',2),hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[13 14 15])
        title('No / Short delay')

        if isempty(i)
        else

            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES(ShortDelay.(myfield1),2),0,3.5, 20,'g','V')

            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES(ShortDelay.(myfield1),1),-2.5,2.5, 20,'g','H')
        end


        %% Ploting reduction of activity

        for i = 1:length(ReductionActivity.(myfield1))
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23 ])
            scatter(COORDINATES(ReductionActivity.(myfield1)(i),1),COORDINATES(ReductionActivity.(myfield1)(i),2),40,'r','filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[19 20 21])
            plot(-CALCIUMroiTS.diff_perc03.data(ReductionActivity.(myfield1)(i),:),'r','LineWidth',2),hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[19 20 21])
        title('Reduction Activity')

        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES(ReductionActivity.(myfield1),2),0,3.5, 20,'r','V')

            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES(ReductionActivity.(myfield1),1),-2.5,2.5, 20,'r','H')
        end


        %% Ploting pre activity

        for i = 1:length(PreActivity.(myfield1))
            figure(1)

            subplot(6,6,[4 5 10 11 16  17 22 23  ])

            scatter(COORDINATES(PreActivity.(myfield1)(i),1),COORDINATES(PreActivity.(myfield1)(i),2),40,'m','filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[25 26 27])
            plot(-CALCIUMroiTS.diff_perc03.data(PreActivity.(myfield1)(i),:),'m','LineWidth',2),hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[25 26 27])
        title('Pre Activity')


        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES(PreActivity.(myfield1),2),0,3.5, 20,'m','V')

            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES(PreActivity.(myfield1),1),-2.5,2.5, 20,'m','H')
        end

        %% Ploting non active

        for i = 1:length(NonActive.(myfield1))
            figure(1)
            subplot(6,6,[4 5 10 11 16  17 22 23  ])
            scatter(COORDINATES(NonActive.(myfield1)(i),1),COORDINATES(NonActive.(myfield1)(i),2),40,'k','filled'), hold on % Coordinates of the lateral refere

            subplot(6,6,[31 32 33])
            plot(-CALCIUMroiTS.diff_perc03.data(NonActive.(myfield1)(i),:),'k','LineWidth',2),hold on
            xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
        end
        subplot(6,6,[31 32 33])
        title('Non Active')


        if isempty(i)
        else
            subplot(6,6,[6 12 18 24 ])
            Hist(COORDINATES([NonActive.(myfield1)],2),0,3.5, 20,'k','V')

            subplot(6,6,[ 28 29 34 35])
            Hist(COORDINATES([NonActive.(myfield1)],1),-2.5,2.5, 20,'k','H')
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
elseif  str==3 % This will be the 3D plot

    figure(1)
    scatter3(brain(:,1),brain(:,2),brain(:,3),'o','MarkerEdgeAlpha','0.05'), hold on
    ylim ([-0.5 3.5])
    scatter3(mauthnerR(:,1),mauthnerR(:,2),mauthnerR(:,3),'b')
    scatter3(mauthnerL(:,1),mauthnerL(:,2),mauthnerL(:,3),'b')
    scatter3(chx10(:,1),chx10(:,2),chx10(:,3),15,[168./255, 168./255, 168./255],'filled')

    %% Ploting all cells
    for i = 1:length(COORDINATES)
        scatter3(COORDINATES((i),1),COORDINATES((i),2),...
            COORDINATES((i),3),40,[0.4940, 0.1840, 0.5560],'filled'), hold on
    end

    %% Ploting long delay
    for i = 1:length(LongDelay)
        scatter3(COORDINATES(LongDelay.(myfield1)(i),1),COORDINATES(LongDelay.(myfield1)(i),2),COORDINATES(LongDelay(i),3),...
            40,	[0.9290, 0.6940, 0.1250],'filled'), hold on
    end

    %% Ploting no/short delay
    for i = 1:length(ShortDelay)
        scatter3(COORDINATES(ShortDelay.(myfield1)(i),1),COORDINATES(ShortDelay.(myfield1)(i),2),COORDINATES(ShortDelay(i),3)...
            ,40,'g','filled'), hold on
    end

    %% Ploting reduction of activity
    for i = 1:length(ReductionActivity)
        scatter3(COORDINATES(ReductionActivity.(myfield1)(i),1),COORDINATES(ReductionActivity.(myfield1)(i),2),COORDINATES(ReductionActivity(i),3),...
            40,'r','filled'), hold on
    end

    %% Ploting pre activity
    for i = 1:length(PreActivity)
        scatter3(COORDINATES(PreActivity.(myfield1)(i),1),COORDINATES(PreActivity.(myfield1)(i),2),COORDINATES(PreActivity(i),3),...
            40,'m','filled'), hold on
    end

    %% Ploting non active
    for i = 1:length(NonActive)
        scatter3(COORDINATES(NonActive.(myfield1)(i),1),COORDINATES(NonActive.(myfield1)(i),2),COORDINATES(NonActive(i),3),...
            40,'k','filled'), hold on
    end

end
end