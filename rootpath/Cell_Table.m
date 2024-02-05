%% David Madrid . Last Rev 26/06/2022
clear
close all
user_settings;
list

%%                              LOADING FISH IN A LOOP

count = 0; 
for iiii=1:length(list)
    

    nfish =iiii; % Number of the fish to load (this is the trial for the same fish)
    Name = strcat('A',num2str(iiii)); 

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    %% For each episode will be created a matrix. The first columns of the matrix (coordinates)
    %  of the cell will be the same, the only it will change is the rest of the
    %  columns that will be specific for each episode of s

    for i=1:length(fieldnames(CALCIUM.EPISODE_X))
        temp = fieldnames(CALCIUM.EPISODE_X );
        myfield = char(temp (i));
        TEMP = table;
        TEMP(:,[1 2 3]) =  array2table(CALCIUM.Coordinates); % Coordinates

        for k = 1:size(CALCIUMroiTS.diff_perc03.data,2)
            my_field1 = strcat('CELL',num2str(k));
            TEMP(k,4) = array2table(CALCIUM.RESPONSE.(my_field1)(i)); % If there is response or not
            if  table2array(TEMP(k,4)) == 0
                TEMP(k,5) = array2table(0);
            else
                TEMP(k,5) =  array2table(CALCIUM.BASELINE_TO_PEAK.(my_field1)(table2array(TEMP(k,4)))); % The amplitud of the response
            end

        end



        TEMP(:,6) = deal(CALCIUM.LABELLINGBH(i)); % Type of behaviour


        TEMP(:,7) = array2table(mean(CALCIUM.FREQ_SWIM_Y.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))'); % Mean Freq
        TEMP(:,8) = array2table(max(CALCIUM.FREQ_SWIM_Y.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))'); % Max Freq
        TEMP(:,9) = array2table(mean(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))'); % Mean amplitud
        TEMP(:,10) = array2table(max(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))'); % Max amplitud
        TEMP(:,11) = array2table((mean(CALCIUM.FREQ_SWIM_Y.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))')./CALCIUM.NMAX_FREQ); % Mean Freq
        TEMP(:,12) = array2table((max(CALCIUM.FREQ_SWIM_Y.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))')./CALCIUM.NMAX_FREQ); % Max Freq
        TEMP(:,13) = array2table((mean(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))')./CALCIUM.NMAX_SWIMMING); % Mean amplitud
        TEMP(:,14) = array2table((max(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield)).*ones(1,size(CALCIUMroiTS.diff_perc03.data,2))')./CALCIUM.NMAX_SWIMMING); % Max amplitud



        TEMP.Properties.VariableNames = {'X Coordinates','Y Coordinates','Z Coordinates','Response','Ca2+ Transient','Type behaviour','Mean Freq','Max Freq'...
            ,'Mean Ampl','Max Ampl','Mean Freq Norm','Max Freq Norm'...
            ,'Mean Ampl Norm','Max Ampl Norm'};

        CALCIUM.DATA.(myfield)= TEMP;
        CALCIUMimg('save',CALCIUM,[],list,nfish);

        ANIMAL(count+1:size(TEMP,1)+count,1:size(TEMP,2)) = TEMP;
        count = count +size(TEMP,1); 


clear TEMP
    end

end
