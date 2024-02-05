
%% Here you have the possibility of selecting which coordinates you want to plot.
% --> If you select CALCIUM.ZAbsoluteCoordinates the neurons will be plot it on
% top of the imaris data
% --> If yo select CalculatedCoordinates the neurons of the imaris data
% will be color coded based on their function and no new neurons will be
% ploted
COORDINATES =CALCIUM.ZAbsoluteCoordinates; 

%%
if isfield(CALCIUM, 'LongDelay') % if 1 --> Then means that it was saved,
    % thus there is data about the neurons active, non actived, etc.
    % Then just reload is needed

    LongDelay =  CALCIUM.LongDelay;
    ShortDelay =  CALCIUM.ShortDelay;
    ReductionActivity=  CALCIUM.ReductionActivity;
    PreActivity=  CALCIUM.PreActivity;
    ActiveNeurons = CALCIUM.ActiveNeurons;
    TotalNeurons = CALCIUM.TotalNeurons;
    NonActive = CALCIUM.NonActive;

else % Means there are no previous data about the neurons,
    % it must be created. For now is based on manual identification
    % of cells DURING SWIMMING. After neurons will be saved

    myfield = 'swm_3'; 

    error('You must introduce the neurons')

    LongDelay.(myfield)=  [];
    ShortDelay.(myfield) = [1 2 3 4 6 7 8 9 10 11 35 48 46 47 45 42];
    ReductionActivity.(myfield)= [] ;
    PreActivity.(myfield)=  [] ;

    ActiveNeurons.(myfield) = [LongDelay.(myfield) ShortDelay.(myfield) ReductionActivity.(myfield) PreActivity.(myfield)]; % Neurons which fall in one of the categories and then are active
    TotalNeurons.(myfield) = 1:length(COORDINATES); % All the neurons recorded
    NonActive.(myfield) =  find(~(ismember(TotalNeurons.(myfield),ActiveNeurons.(myfield))));  % % Neurons which do not fall in any category and then are not active

    CALCIUM.LongDelay = LongDelay ;
    CALCIUM.ShortDelay = ShortDelay ;
    CALCIUM.ReductionActivity = ReductionActivity ;
    CALCIUM.PreActivity = PreActivity ;
    CALCIUM.ActiveNeurons = ActiveNeurons ;
    CALCIUM.TotalNeurons = TotalNeurons;
    CALCIUM.NonActive = NonActive;

    CALCIUMimg('save',CALCIUM,[],list,nfish);
end


%% Loading raw data from imaris

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

prompt = strcat('You have','_',num2str(length(fieldnames(CALCIUM.EPISODE_X))),'_',' which one do you want to plot? ex: swm_1, swm_3,...');
str1 = input(prompt,'s');


if str==2 % This will be 2D ploting

    figure(1)

    % Ploting scafold
    subplot(6,6,[4 5 10 11 16  17 22 23]) %subplot(5,5,[4  9  14  19  24 25])

    % Ploting brain
    scatter(brain(:,1),brain(:,2),'o','MarkerEdgeAlpha','0.05')

    % Holding
    hold on

    % Ploting chx10 population
    scatter(chx10(:,1),chx10(:,2),15,[168./255, 168./255, 168./255],'filled')

    % Ploting mauthner cells
    scatter(mauthnerR(:,1),mauthnerR(:,2),'b')
    scatter(mauthnerL(:,1),mauthnerL(:,2),'b')

    % Rescaling
    ylim ([-0.5 3.5])


    % Correction problem with the trasposition
    if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
        CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
        CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
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


    for i = 1:length(LongDelay.(str1))
        figure(1)
        subplot(6,6,[4 5 10 11 16  17 22 23 ])
        scatter(COORDINATES(LongDelay.(str1)(i),1),COORDINATES(LongDelay.(str1)(i),2),40,	[0.9290, 0.6940, 0.1250],'filled'), hold on % Coordinates of the lateral refere

        subplot(6,6,[7 8 9])
        plot(-CALCIUMroiTS.diff_perc03.data(LongDelay.(str1)(i),:),'color',[0.9290, 0.6940, 0.1250],'LineWidth',2), hold on
        xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    end
    subplot(6,6,[7 8 9])
    title('Long delay')

    if isempty(i)
    else
        subplot(6,6,[6 12 18 24 ])
        Hist(COORDINATES(LongDelay.(str1),2),0,3.5, 20,[0.9290, 0.6940, 0.1250],'V')

        subplot(6,6,[ 28 29 34 35])
        Hist(COORDINATES(LongDelay.(str1),1),-2.5,2.5, 20,[0.9290, 0.6940, 0.1250],'H')
    end



    %% Ploting no/short delay

    for i = 1:length(ShortDelay.(str1))
        figure(1)
        subplot(6,6,[4 5 10 11 16  17 22 23 ])
        scatter(COORDINATES(ShortDelay.(str1)(i),1),COORDINATES(ShortDelay.(str1)(i),2),40,'g','filled'), hold on % Coordinates of the lateral refere

        subplot(6,6,[13 14 15])
        plot(-CALCIUMroiTS.diff_perc03.data(ShortDelay.(str1)(i),:),'g','LineWidth',2),hold on
        xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    end
    subplot(6,6,[13 14 15])
    title('No / Short delay')

    if isempty(i)
    else

        subplot(6,6,[6 12 18 24 ])
        Hist(COORDINATES(ShortDelay.(str1),2),0,3.5, 20,'g','V')

        subplot(6,6,[ 28 29 34 35])
        Hist(COORDINATES(ShortDelay.(str1),1),-2.5,2.5, 20,'g','H')
    end


    %% Ploting reduction of activity

    for i = 1:length(ReductionActivity.(str1))
        figure(1)
        subplot(6,6,[4 5 10 11 16  17 22 23 ])
        scatter(COORDINATES(ReductionActivity.(str1)(i),1),COORDINATES(ReductionActivity.(str1)(i),2),40,'r','filled'), hold on % Coordinates of the lateral refere

        subplot(6,6,[19 20 21])
        plot(-CALCIUMroiTS.diff_perc03.data(ReductionActivity.(str1)(i),:),'r','LineWidth',2),hold on
        xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    end
    subplot(6,6,[19 20 21])
    title('Reduction Activity')

    if isempty(i)
    else
        subplot(6,6,[6 12 18 24 ])
        Hist(COORDINATES(ReductionActivity.(str1),2),0,3.5, 20,'r','V')

        subplot(6,6,[ 28 29 34 35])
        Hist(COORDINATES(ReductionActivity.(str1),1),-2.5,2.5, 20,'r','H')
    end


    %% Ploting pre activity

    for i = 1:length(PreActivity.(str1))
        figure(1)

        subplot(6,6,[4 5 10 11 16  17 22 23  ])

        scatter(COORDINATES(PreActivity.(str1)(i),1),COORDINATES(PreActivity.(str1)(i),2),40,'m','filled'), hold on % Coordinates of the lateral refere

        subplot(6,6,[25 26 27])
        plot(-CALCIUMroiTS.diff_perc03.data(PreActivity.(str1)(i),:),'m','LineWidth',2),hold on
        xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    end
    subplot(6,6,[25 26 27])
    title('Pre Activity')


    if isempty(i)
    else
        subplot(6,6,[6 12 18 24 ])
        Hist(COORDINATES(PreActivity.(str1),2),0,3.5, 20,'m','V')

        subplot(6,6,[ 28 29 34 35])
        Hist(COORDINATES(PreActivity.(str1),1),-2.5,2.5, 20,'m','H')
    end

    %% Ploting non active

    for i = 1:length(NonActive.(str1))
        figure(1)
        subplot(6,6,[4 5 10 11 16  17 22 23  ])
        scatter(COORDINATES(NonActive.(str1)(i),1),COORDINATES(NonActive.(str1)(i),2),40,'k','filled'), hold on % Coordinates of the lateral refere

        subplot(6,6,[31 32 33])
        plot(-CALCIUMroiTS.diff_perc03.data(NonActive.(str1)(i),:),'k','LineWidth',2),hold on
        xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    end
    subplot(6,6,[31 32 33])
    title('Non Active')


    if isempty(i)
    else
        subplot(6,6,[6 12 18 24 ])
        Hist(COORDINATES([NonActive.(str1)],2),0,3.5, 20,'k','V')

        subplot(6,6,[ 28 29 34 35])
        Hist(COORDINATES([NonActive.(str1)],1),-2.5,2.5, 20,'k','H')
    end


    %% Creating folder to saved the plot
    if exist('CATEGORY', 'dir')
    else
        mkdir (fullfile(path.data,'dataCALCIUM',CALCIUM.ref,'CATEGORY')); % Creating the directory for dataCALCIUM
    end
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),'CATEGORY',strcat('CATEGORIZED_',char(CALCIUM.list(nfish,1)),'.fig')))
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),'CATEGORY',strcat('CATEGORIZED_',char(CALCIUM.list(nfish,1)),'.eps')))
    uiwait

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
        scatter3(COORDINATES(LongDelay.(str1)(i),1),COORDINATES(LongDelay.(str1)(i),2),COORDINATES(LongDelay(i),3),...
            40,	[0.9290, 0.6940, 0.1250],'filled'), hold on 
    end

    %% Ploting no/short delay
    for i = 1:length(ShortDelay)
        scatter3(COORDINATES(ShortDelay.(str1)(i),1),COORDINATES(ShortDelay.(str1)(i),2),COORDINATES(ShortDelay(i),3)...
            ,40,'g','filled'), hold on 
    end

 %% Ploting reduction of activity
    for i = 1:length(ReductionActivity)
        scatter3(COORDINATES(ReductionActivity.(str1)(i),1),COORDINATES(ReductionActivity.(str1)(i),2),COORDINATES(ReductionActivity(i),3),...
            40,'r','filled'), hold on 
    end

     %% Ploting pre activity
    for i = 1:length(PreActivity)
        scatter3(COORDINATES(PreActivity.(str1)(i),1),COORDINATES(PreActivity.(str1)(i),2),COORDINATES(PreActivity(i),3),...
            40,'m','filled'), hold on 
    end

     %% Ploting non active
    for i = 1:length(NonActive)
        scatter3(COORDINATES(NonActive.(str1)(i),1),COORDINATES(NonActive.(str1)(i),2),COORDINATES(NonActive(i),3),...
            40,'k','filled'), hold on 
    end

end


   
    
  
   
