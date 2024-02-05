%% Correction problem with the trasposition
if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
    CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
    %CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
else
end

%% Filtering Calcium imaging movie (median filter)
close all
% Median value by default
median = 10;

%----------------------------------------
% Control: is the data already filtered ?
% If 'data_No_filtered' exits as a field means that the data was filtered
% already and there is a copy of the original data.

if isfield(CALCIUMroiTS.diff_perc03,'data_No_filtered') % Checking if the field exits
    prompt = {['You have already filtered your data, do you want to get the original data back? if yes type Y, if no type N']};
    dlgtitle = 'Input';
    dims = [1 35];
    definput ={'N'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if isequal(answer,{'Y'}) % If yes, then asking
        CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data_No_filtered; %Getting back the non filtered data
         error('You must erase manually CALCIUMroiTS.diff_perc03.data_No_filtered')
    else
    end

else  % If is not a field

    CALCIUMroiTS.diff_perc03.data_No_filtered = CALCIUMroiTS.diff_perc03.data; %This will
    % save the data that is not filtered

    % 1. Visualization of original data

    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(211)
    plot(-CALCIUMroiTS.diff_perc03.data')
    %plot(-CALCIUMroiTS.diff_perc03.data_No_filtered')
    title('Non filtered data')

    % User interface
    prompt = {['Do you wanna filter the data?: if yes type Y, if no type N']};
    dlgtitle = 'Input';
    dims = [1 35];
    definput ={'Y'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    if isequal(answer,{'Y'})

        % Applying the median filter and visualization
        for i = 1:size(CALCIUMroiTS.diff_perc03.data,1)
            temporal_wave(i,:) = medfilt1(CALCIUMroiTS.diff_perc03.data(i,:),median);
        end
        subplot(212)
        plot(-temporal_wave')
        title('filtered data')

        prompt = 'If you want to keep filtering type "Y", if you are satisfy type "N" : ';
        str = input(prompt,'s');

        if isempty(str)
            str = 'N'; % you are done
            close
        else
            while isequal(str,'Y')

                prompt = {'Change median threshold?'};
                dlgtitle = 'Input';
                dims = [1 35];

                definput ={sprintf('%.0f',median)};
                answer = inputdlg(prompt,dlgtitle,dims,definput);
                median = str2double(answer(1));

                for i = 1:size(CALCIUMroiTS.diff_perc03.data,1)
                    temporal_wave(i,:) = medfilt1(CALCIUMroiTS.diff_perc03.data(i,:),median);
                end

                subplot(212)
                plot(-temporal_wave)


                prompt = 'If you want to keep filtering type "Y", if you are satisfy type "N" : ';
                str = input(prompt,'s');
            end

        end
        CALCIUMroiTS.diff_perc03.data = temporal_wave;
        CALCIUMroiTS.diff_perc03.Filter_Applied = median;
        CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish); %Saving the data filtered

    elseif isequal(answer,{'N'})
    end
end

close all
clear temporal_wave

%% Creating Linear regime
% The linear regime allow us to identify in a clearer way the general
% shape of the transient and also to detect easily the onset of it

close all
clear TF_table
Type = 'Linear'; 
Th = 2000; 
segline = NaN(size(CALCIUMroiTS.diff_perc03.data));
for i = 1:size(CALCIUMroiTS.diff_perc03.data,1)
    A = -CALCIUMroiTS.diff_perc03.data(i,:);
    [TF,S1,S2] = ischange(A,Type,'Threshold',Th);
    segline(i,:) = S1.*(1:length(A)) + S2;
    plot(1:length(A),A,1:length(A),segline(i,:))
    legend('Data','Linear Regime')
    TF_table(:,i) = array2table(TF'); % This indicates when the sudden transitions occur
   uiwait
end
CALCIUMroiTS.diff_perc03.Linear_Regime = segline; 
CALCIUMroiTS.diff_perc03.Linear_Regime_Type = Type;
CALCIUMroiTS.diff_perc03.Linear_Regime_Th = Th; 
CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish); %Saving the data filtered
%imagesc(segline)


%% Amplitudes of the peaks and onset

close all
CALCIUMpeak = 60; % @SET value
for i = 1:size(CALCIUMroiTS.diff_perc03.data,1) % Loop across all the cells

    % First we indentify the peaks of the CALCIUM WAVE
    figure('units','normalized','outerposition',[0 0 1 1])
    CALCIUMpeak = 60; % @SET value
    findpeaks(-CALCIUMroiTS.diff_perc03.data(i,:),'MinPeakProminence',CALCIUMpeak); hold on

    %Maybe is posible that you want to change the value of the threshold
    %for the peak per each cell so we will open a prompt for introducing
    %the new threshold value

    prompt = {'Enter Threshold:','Are you satisfy with the result?:'};
    dlgtitle = 'Input';
    dims = [1 35];

    definput ={sprintf('%.0f',CALCIUMpeak),'N'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    CALCIUMpeak = str2double(answer(1));
    cla

    % Conditional to repit again the process if the answer
    % is no. In this part, you can change dynamicly the value of the
    % threshold in the GUI

    while isequal(answer(2),{'N'})
        findpeaks(-CALCIUMroiTS.diff_perc03.data(i,:),'MinPeakProminence',CALCIUMpeak); hold on
        prompt = {'Enter Threshold:','Are you satisfy with the result?:'};
        dlgtitle = 'Input';
        dims = [1 35];

        definput ={sprintf('%.0f',CALCIUMpeak),'N'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        CALCIUMpeak = str2double(answer(1));
        cla
    end

    findpeaks(-CALCIUMroiTS.diff_perc03.data(i,:),'MinPeakProminence',CALCIUMpeak); hold on


    % Once the threshold is set, we can remove some peaks. First we create
    % the variable peaks and locs
    [pksCAL,locsCAL] = findpeaks(-CALCIUMroiTS.diff_perc03.data(i,:),'MinPeakProminence',CALCIUMpeak);

    %Creating a dialog box


    prompt = {'You want to remove some peaks:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput ={'Y'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    while isequal(answer(1),{'Y'})
        [x_time,~] = getpts;

        % Finding closest match
        removing = NaN(1,length(x_time)); % Prelocate
        for j=1:length(x_time)
            removing(j) = find(abs(locsCAL-x_time(j))==min(abs(locsCAL-x_time(j))));
        end

        % Removing this index values from the peaks and locs vectors
        pksCAL(removing) = [];
        locsCAL(removing) = [];
        cla
        findpeaks(-CALCIUMroiTS.diff_perc03.data(i,:),'MinPeakProminence',CALCIUMpeak);
        scatter(locsCAL,pksCAL,30,'filled','g')
        pause(0.5)
        prompt = {'You want to remove some peaks:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput ={'Y'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
    end

    % Testing

    scatter(locsCAL,pksCAL,30,'filled','g')

    clear x_time y_time removing
    clear x_time

    % Because each cell can have diferent number of peaks, we will store
    % the data in a structure


    my_field = strcat('CELL',num2str(i)); % Creating dinamicly the

    AMPLITUDE_Y.(my_field) = pksCAL; % we put
    AMPLITUDE_X.(my_field) = CALCIUMroiTS.diff_perc03.times(locsCAL); % we put



%     ALL THE COMMENTED CODE BELLOW IS FOR THE ONSET OF THE CALCIUM
%     ACTIVITY, IT IS NOT READY YET
% 
%     % Finding onset of the calcium wave
% 
%     for kk = 1:length(pksCAL) % Loop across waves (peaks)
% 
%         % Visualization
%         plot(segline(i,:)) % Ploting linear regime
%         plot(diff(segline(i,:))) % Ploting slope of linear regime
% 
%         Increase_Slope = pksCAL(kk)*0.01; % This is the value that the slope of the segline
%         % must take to be detected
% 
%         temp = find(diff(segline(i,:))>Increase_Slope); 
%         V1 = find(diff(temp)==1); 
%         temp(V1(1))
%         scatter( temp(V1(1)),segline(i,temp(V1(1))),80,'filled')
% 
% % 
% %         slope = find(diff(segline(i,:))>10) % The median filter
% %         if size(slope,2)>1
% %             slope2 = find( medfilt1(diff(segline(i,:)))>1)
% %             scatter(slope2(1),segline(i,slope(1)),70,'filled')
% %         else
% %             % will remove sudden but briefs changes in the slope
% %             scatter(slope(1),segline(i,slope(1)),70,'filled')
% %         end
% %         CALCIUM.ONSET_X.(my_field)(kk) = ;
% % 
% %         scatter(,'filled','r')
% 
% 
%     end


saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
    char(CALCIUM.list(nfish,1)),strcat('Peaks_and_onset',my_field, char(CALCIUM.list(nfish,1)))))




    % Finally we clear the pksCAL and locsCAL
    clear locsCAL pksCAL sumation
    w = waitforbuttonpress; % Once the scatter is plotted just click any key to go to next celL
    close all


end

%Saving
CALCIUM.AMPLITUDE_Y = AMPLITUDE_Y;
CALCIUM.AMPLITUDE_X = AMPLITUDE_X;
CALCIUMimg('save',CALCIUM,[],list,nfish);

%%                   EXTRACTING THE BASLINE FOR THE CALCIUM IMAGING

% Extracting the baseline for each cell in order to calculate the correct
% amplitude to the peak.

%Prelocation of the variables
baseline_frames = NaN(2,size(CALCIUMroiTS.diff_perc03.data,2));
xi_baseline_cutted = NaN(2,size(CALCIUMroiTS.diff_perc03.data,2));
TEMP = CALCIUM.AMPLITUDE_X;

% Loop across cells
for i = 1:size(CALCIUMroiTS.diff_perc03.data,1)
    temp2 = fieldnames(CALCIUM.AMPLITUDE_X);
    myfield2 = char(temp2(i)); % extracting the name of each episode
    tempor = 1:length(CALCIUM.AMPLITUDE_X.(myfield2));

    if size(tempor)==[1 0] % if not transiente
        jjj=1
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(i,:))
        disp('Click once at the begining of the baseline, and twice at the end')
        [x_time,~] = getpts; %@

        % Control just to be sure we select 2 points
        if length(x_time)>2
            error('You click more the two points')
        elseif length(x_time)<2
            error('You click just the first point')
        end


        baseline_frames (:,i) = x_time; %Saving value in the variable baseline_frames
        %
        %     round(baseline_frames) % to have the value of the frame

        for j=1:size(baseline_frames,1) %Finding the index in the time vector that is
            % closer to the values we just pick

            substraction = abs(CALCIUMroiTS.diff_perc03.times-baseline_frames (j,i));
            xi_baseline_cutted (j,i) = find(substraction==min(substraction));
            clear substraction % Erase temporal variable
        end
        temp = mean(-CALCIUMroiTS.diff_perc03.data(i,[xi_baseline_cutted(1,i):xi_baseline_cutted(2,i)]));
        baseline_CALCIUM.(myfield2)(jjj) = temp;
        %                 clear temp

        close

    else
        for jj = 1:length(CALCIUM.AMPLITUDE_X.(myfield2))


            figure('units','normalized','outerposition',[0 0 1 1])
            plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(i,:)), hold on
            scatter(CALCIUM.AMPLITUDE_X.(myfield2)(jj),CALCIUM.AMPLITUDE_Y.(myfield2)(jj),100)
            title(strcat('baseline for the transient',num2str(jj)))
            disp('Click once at the begining of the baseline, and twice at the end')
            [x_time,~] = getpts; %@

            % Control just to be sure we select 2 points
            if length(x_time)>2
                error('You click more the two points')
            elseif length(x_time)<2
                error('You click just the first point')
            end


            baseline_frames (:,i) = x_time; %Saving value in the variable baseline_frames
            %
            %     round(baseline_frames) % to have the value of the frame

            for j=1:size(baseline_frames,1) %Finding the index in the time vector that is
                % closer to the values we just pick

                substraction = abs(CALCIUMroiTS.diff_perc03.times-baseline_frames (j,i));
                xi_baseline_cutted (j,i) = find(substraction==min(substraction));
                clear substraction % Erase temporal variable
            end

            %Computing value of the baseline

            %temp will be a temporary vector that include the values of the calcium
            %wave btw the specific index that we selected. we will calculate the
            %average of this to have the avrg baseline for EACH cell.

            temp = mean(-CALCIUMroiTS.diff_perc03.data(i,[xi_baseline_cutted(1,i):xi_baseline_cutted(2,i)]));
            baseline_CALCIUM.(myfield2)(jj) = temp;
            clear temp

            close
        end
    end
end

user_settings
CALCIUM.baseline_amplitude = baseline_CALCIUM;
CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie


%%                    BASELINE TO PEAK FOR EACH CELL CALCIUM IMAGING

%Having the value of the baseline for each cell and having also all the
%peaks in amplitude for each cell, we can calculate now the
%baseline-to-peak distance for each cell. the number of distance per cell
%will be equal to the number of amplitude peaks. Because each cell can have
%different peaks we can not save all the cells in a matrix, we have to use
%a structure

% Loop across cells
for k = 1:size(CALCIUMroiTS.diff_perc03.data,1)

    my_field = strcat('CELL',num2str(k)); % Creating dinamicly the name of the
    % cell

    % Performing the baseline to peak per each cell, for that we use the
    % variable that we created perviously called AMPLITUDE_Y.(my_field)
    % where my field is CELL1,CELL2,... So we are performing the
    % substraction for the baseline of each cell with all the peaks of the
    % same cell to have the baseline to peak for that cell, and we are
    % looping across cell

    temp = CALCIUM.AMPLITUDE_Y.(my_field)'-CALCIUM.baseline_amplitude.(my_field)';
    if isempty(temp)

        BASELINE_TO_PEAK.(my_field) = 0; 
    else
        BASELINE_TO_PEAK.(my_field) = temp;
    end


end

% Saving

CALCIUM.BASELINE_TO_PEAK = BASELINE_TO_PEAK;
CALCIUMimg('save',CALCIUM,[],list,nfish);





% %% SET UP A THRESHOLD FOR THE CALCIUM (IN THIS RECORDING)
% % When you are detecting the peaks of each cell the y axis of each cell
% % is diferent and then you dont know if the peak that you are selecting
% % is big or small. Set a threshold that will be the 20% of the maximum
% % response any cell in the recording
% 
% 
% for i = 1:length(fieldnames(CALCIUM.BASELINE_TO_PEAK))
%     temp5 = fieldnames(CALCIUM.BASELINE_TO_PEAK);
%     temp6 = char(temp5(i));
%     if isempty(CALCIUM.BASELINE_TO_PEAK.(temp6))
%     else
%         maxcalcium_response(i) = max(CALCIUM.BASELINE_TO_PEAK.(temp6));
%     end
% end
% 
% MAXCALCIUM_RESPONSE = max(maxcalcium_response);
% CALCIUM.MAXCALCIUM_RESPONSE = MAXCALCIUM_RESPONSE;
% CALCIUMimg('save',CALCIUM,[],list,nfish);




