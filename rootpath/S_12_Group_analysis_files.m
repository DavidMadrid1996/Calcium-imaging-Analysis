%Last revision: 2023/11/29

%% Data extraction.
clear
close all
S_11_Group_analysis_files % @Check it

%% Prelocating variables
suma_coordinates = 0; 
suma_LongDelay = 0; 
suma_ShortDelay = 0;
suma_ReductionActivity = 0;
suma_PreActivity =0 ;
suma_ActiveNeurons = 0;
suma_TotalNeurons = 0;
suma_NonActive = 0;

% Will create a matrix which will allow us to aline all
% the neurons that we want. All neurons will be aline with the swimming
% starting at the middle column
Data_ShortDelay_Blank = zeros(500,5000); 
Data_LongDelay_Blank = zeros(500,5000); 
Data_COORDINATES_Blank = zeros(500,5000); 
m = 2500; 

%%


    for j = 1:size(ValidPathdataCalcium,1)

        % loading animal
        load(char(fullfile(ValidPathdataCalcium(j,:),strcat('CALCIUM_',strcat(Animal_List_unique(j,1),...
            Animal_List_unique(j,2),'.mat')))));
        [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(j,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(j,1),...
            Animal_List_unique(j,2),'.mat')))));
        [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
        myfield1 = char(Animal_List_unique(j,3));
        COORDINATES =CALCIUM.ZAbsoluteCoordinates;


         % Extracting delay vector. 
         clear temp_Delay
        for i = 1:length(fieldnames(CALCIUM.DELAY))
            temp = fieldnames(CALCIUM.DELAY); 
            temp2 = char(temp (i));
            temp3 = nonzeros(CALCIUM.DELAY.(temp2).(myfield1)); 

            if isempty(temp3) |  isnan(temp3)
                temp3 = NaN; 
                temp_Delay(i) = temp3; 

            elseif temp3>10
                temp3 = NaN; 
                temp_Delay(i) = temp3; 

                elseif temp3<0
                temp3 = NaN; 
                temp_Delay(i) = temp3; 

            else 
                temp4 =  temp3>0; 
                temp5 = temp3(temp4); 
                temp_Delay(i) =temp5(1) ;
            end
        end

                % Correction problem with the trasposition
        if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
            CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
%             CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
        else
        end


        % Finding the time point in the calcium imaging when the swim
        % starts
        time_for_raster = find(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)(1))>CALCIUM.timebase);
        how_much_frames = 50; 

        %

        TOTAL_COORDINATES(1+suma_coordinates:length(COORDINATES)+suma_coordinates,:) =  COORDINATES(:,:);

        Data_COORDINATES(1+suma_coordinates:length(COORDINATES)+suma_coordinates,:) = ...
            CALCIUMroiTS.diff_perc03.data(:,(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

        Data_COORDINATES_Blank(1+suma_coordinates:length(COORDINATES)+suma_coordinates, m-time_for_raster(end):m+size(CALCIUMroiTS.diff_perc03.data,2)-time_for_raster(end)-1) =...
            CALCIUMroiTS.diff_perc03.data ; 
          


        suma_coordinates= length(COORDINATES) + suma_coordinates ; 


        %
        LongDelay_COORDINATES(1+suma_LongDelay:length(COORDINATES(CALCIUM.LongDelay.(myfield1)))+suma_LongDelay,:) =  COORDINATES(CALCIUM.LongDelay.(myfield1),:);

        Data_LongDelay(1+suma_LongDelay:length(COORDINATES(CALCIUM.LongDelay.(myfield1)))+suma_LongDelay,:) = ...
            CALCIUMroiTS.diff_perc03.data(CALCIUM.LongDelay.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

        
        
         Data_LongDelay_Blank(1+suma_LongDelay:length(CALCIUM.LongDelay.(myfield1))+suma_LongDelay, m-time_for_raster(end):m+size(CALCIUMroiTS.diff_perc03.data,2)-time_for_raster(end)-1) = ...
            CALCIUMroiTS.diff_perc03.data(CALCIUM.LongDelay.(myfield1),:); % This is the most general

        for h = 1:length(COORDINATES(CALCIUM.LongDelay.(myfield1)))
            Name__LongDelay{h+suma_LongDelay,:} = strcat('CALCIUM_',strcat(Animal_List_unique(j,1),...
                Animal_List_unique(j,2),'.mat'),'_',myfield1,'_LongDelay','_Cell','_', num2str(CALCIUM.LongDelay.(myfield1)(h))); 
        end


       

        Time_LongDelay(1+suma_LongDelay:length(COORDINATES(CALCIUM.LongDelay.(myfield1)))+suma_LongDelay,:) = temp_Delay(CALCIUM.LongDelay.(myfield1))'; 
        suma_LongDelay = length(COORDINATES(CALCIUM.LongDelay.(myfield1))) + suma_LongDelay ;



        %
        ShortDelay_COORDINATES(1+suma_ShortDelay:length(COORDINATES(CALCIUM.ShortDelay.(myfield1)))+suma_ShortDelay,:) =  COORDINATES(CALCIUM.ShortDelay.(myfield1),:);
        
        
        Data_ShortDelay(1+suma_ShortDelay:length(COORDINATES(CALCIUM.ShortDelay.(myfield1)))+suma_ShortDelay,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.ShortDelay.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);


         Data_ShortDelay_Blank(1+suma_ShortDelay:length(CALCIUM.ShortDelay.(myfield1))+suma_ShortDelay, m-time_for_raster(end):m+size(CALCIUMroiTS.diff_perc03.data,2)-time_for_raster(end)-1) = ...
            CALCIUMroiTS.diff_perc03.data(CALCIUM.ShortDelay.(myfield1),:); % This is the most general

         for h = 1:length(COORDINATES(CALCIUM.ShortDelay.(myfield1)))
             Name__ShortDelay{h+suma_ShortDelay,:} = strcat('CALCIUM_',strcat(Animal_List_unique(j,1),...
                 Animal_List_unique(j,2),'.mat'),'_',myfield1,'_ShortDelay','_Cell','_', num2str(CALCIUM.ShortDelay.(myfield1)(h))); 
         end

        
        Time_ShortDelay(1+suma_ShortDelay:length(COORDINATES(CALCIUM.ShortDelay.(myfield1)))+suma_ShortDelay,:) = temp_Delay(CALCIUM.ShortDelay.(myfield1))'; 
        suma_ShortDelay= length(COORDINATES(CALCIUM.ShortDelay.(myfield1))) + suma_ShortDelay ;


        
        %
        ReductionActivity_COORDINATES(1+suma_ReductionActivity:length(COORDINATES(CALCIUM.ReductionActivity.(myfield1)))+suma_ReductionActivity,:) =  COORDINATES(CALCIUM.ReductionActivity.(myfield1),:);
        
         Data_ReductionActivity(1+suma_ReductionActivity:length(COORDINATES(CALCIUM.ReductionActivity.(myfield1)))+suma_ReductionActivity,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.ReductionActivity.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

        suma_ReductionActivity = length(COORDINATES(CALCIUM.ReductionActivity.(myfield1))) + suma_ReductionActivity ;


         %
        PreActivity_COORDINATES(1+suma_PreActivity:length(COORDINATES(CALCIUM.PreActivity.(myfield1)))+suma_PreActivity,:) =  COORDINATES(CALCIUM.PreActivity.(myfield1),:);
        
         Data_PreActivity(1+suma_PreActivity:length(COORDINATES(CALCIUM.PreActivity.(myfield1)))+suma_PreActivity,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.PreActivity.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

        suma_PreActivity = length(COORDINATES(CALCIUM.PreActivity.(myfield1))) + suma_PreActivity ; 


        %
        ActiveNeurons_COORDINATES(1+suma_ActiveNeurons:length(COORDINATES(CALCIUM.ActiveNeurons.(myfield1)))+suma_ActiveNeurons,:) =  COORDINATES(CALCIUM.ActiveNeurons.(myfield1),:);
        
        Data_ActiveNeurons(1+suma_ActiveNeurons:length(COORDINATES(CALCIUM.ActiveNeurons.(myfield1)))+suma_ActiveNeurons,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.ActiveNeurons.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

        suma_ActiveNeurons = length(COORDINATES(CALCIUM.ActiveNeurons.(myfield1))) + suma_ActiveNeurons ; 


         %
        
         TotalNeurons_COORDINATES(1+suma_TotalNeurons:length(COORDINATES(CALCIUM.TotalNeurons.(myfield1)))+suma_TotalNeurons,:) =  COORDINATES(CALCIUM.TotalNeurons.(myfield1),:);
        
          Data_TotalNeurons(1+suma_TotalNeurons:length(COORDINATES(CALCIUM.TotalNeurons.(myfield1)))+suma_TotalNeurons,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.TotalNeurons.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);

         suma_TotalNeurons = length(COORDINATES(CALCIUM.TotalNeurons.(myfield1))) + suma_TotalNeurons ; 


        %
        NonActive_COORDINATES(1+suma_NonActive:length(COORDINATES(CALCIUM.NonActive.(myfield1)))+suma_NonActive,:) =  COORDINATES(CALCIUM.NonActive.(myfield1),:);
        
         Data_NonActive(1+suma_NonActive:length(COORDINATES(CALCIUM.NonActive.(myfield1)))+suma_NonActive,:) = ...
        CALCIUMroiTS.diff_perc03.data(CALCIUM.NonActive.(myfield1),(time_for_raster(end)-1):time_for_raster(end)+how_much_frames);


         for h = 1:length(COORDINATES(CALCIUM.NonActive.(myfield1)))
             Name__NonActive{h+suma_NonActive,:} = strcat('CALCIUM_',strcat(Animal_List_unique(j,1),...
                 Animal_List_unique(j,2),'.mat'),'_',myfield1,'_NonActiveNonActive','_Cell','_', num2str(CALCIUM.NonActive.(myfield1)(h))); 
         end

        suma_NonActive = length(COORDINATES(CALCIUM.NonActive.(myfield1))) + suma_NonActive ; 





    end 