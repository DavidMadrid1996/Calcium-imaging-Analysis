%Last revision: 2023/10/21

%% Integration calcium imaging and deeplabcut signal

% clearing variables that it will be used
close all
clear WINDOW_X WINDOW_X_SWIM WINDOW_Y WINDOW_Y_SWIM DELAY


% Setting time thresholds:

% 1. We have a period of x seconds after the finish of the swimming episode for the calcium signal to appear
cte = 3; % [s]

% 2. Having the CALCIUM peak we will take a retrospectic window period of x secs
restric = 3; % [s]

% Stabilising temp variable
TEMP = CALCIUM.AMPLITUDE_X;


% transposition correction
if size(CALCIUMroiTS.diff_perc03.data,1) > size(CALCIUMroiTS.diff_perc03.data,2)
    CALCIUMroiTS.diff_perc03.data=CALCIUMroiTS.diff_perc03.data';
else
end


% First wi will loop across the cells
for i = 1:length(fieldnames(CALCIUM.AMPLITUDE_X))

    % Extracting cell names
    temp1 = fieldnames(CALCIUM.AMPLITUDE_X);

    % Extracting the fieldnames of eathe current cell
    myfield1 = char(temp1(i));

    % The CALCIUM.AngHeadTail3D will be used as swim trace
    wave = squeeze(CALCIUM.AngHeadTail3D);


    % If there is no peak was detected in the trace for this cell means
    % that the cell was not responding. In this case the variable DELAY
    % (this is the delay between the peak of the transiente and the swimming)
    % will be NaN for the non responding cell
    if isempty(CALCIUM.AMPLITUDE_X.(myfield1)) % so if is true DELAY for all episodes
        % of swimming of this cell will be NaN
        for h= 1:length(fieldnames(CALCIUM.EPISODE_X))  % Loop across the number of episodes
            temp2 = fieldnames(CALCIUM.EPISODE_X);
            myfield2 = char(temp2(h));
            DELAY.(myfield1).(myfield2)(1) = NaN; % Setting DELAY per each episode to NaN
            RESPONSE.(myfield1).(myfield2)(1) = NaN; 
        end


    else % In the case that is not empty means that the cell responded at least once

        % We will loop then across all the response to see which episode of
        % swimming is responsable of the this response and what is the time
        % delay between the response and the swimming episode

        yyaxis right
        plot((CALCIUMroiTS.diff_perc03.times),-CALCIUMroiTS.diff_perc03.data(i,:),'b'), hold on
        yyaxis left
        plot(CALCIUMroiTS.deeplabcut.time_scaled, wave-CALCIUM.baseline_swimming,'g','LineWidth',1);

        for j=1:length(CALCIUM.AMPLITUDE_X.(myfield1)) % number of responses of the current cell

        % Scatter current response    
        yyaxis right
        scatter((CALCIUM.AMPLITUDE_X.(myfield1)(j)),CALCIUM.AMPLITUDE_Y.(myfield1)(j) ,100,'magenta','filled')


            % Restriction to the left: from the current calcium response to
            % x seconds into the past
            leftrstr = (TEMP.(myfield1)(j))-restric ;
            window = [leftrstr (TEMP.(myfield1)(j))]; % Creating the retrospectic window

            % Swimming episodes loop
            for h = 1:length(fieldnames(CALCIUM.EPISODE_X))
                % Extracting the name of each episode
                temp2 = fieldnames(CALCIUM.EPISODE_X);
                % Selecting the current episode of swimming
                myfield2 = char(temp2(h)); 

                % Now we create a variable that contain the time points in
                % seconds of the current swimming episode
                unico = (CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield2))); 
                % we add x constant number of seconds defined at the beginning
                ctevalues = unico(end) + cte; 

                % Finally we create a propspectiv time window that will
                % start with the beginning of the swimming episode and it
                % will move forward the length of the episode plus the
                % constant x number of seconds
                safetime = [unico(1) ctevalues];  % Creating the prospective window

                % Next we defines a variable that will be safe as the
                % moment in the calcium recording when the movement occurs
                CA_movement(i,h) = unico(1);

                % We define the variable delay for the current response of
                % the cell with respect to the current swimming episode:
                % Prelocating
                DELAY.(myfield1).(myfield2)(j) = zeros(1,1);

                % So this variable will respresent for the current cell,
                % for the current episode of swimming, if the peak of
                % transient "colocalized" temporally with the swimming
                % episode. because you are more recently inside a Swimming episodes loop
                % means that in the second iteration of the loop (it will be for example swm_3)
                % the variable will produce the delay between the current
                % peak and this next episode of swimming. So how the DELAY
                % variable must be understand is for example:

                % DELAY.CELL1.swm_1= [1.194140849420853,8.097140849420853,11.256140849420852]
                % DELAY.CELL1.swm_3= [-8.909057606177596,-2.006057606177595,1.152942393822404]

                % means that the first peak occurs 1.1941 second after swm_1 
                % and 8.9090 seconds (--8.9090) before of swm_3. So the columns
                % indicate the delay of each peak with respect to that
                % swimming episode
                DELAY.(myfield1).(myfield2)(j) = (TEMP.(myfield1)(j)) -CA_movement(i,h);


                % This is for the response structure


                if (TEMP.(myfield1)(j))>=safetime(1)...
                        && (TEMP.(myfield1)(j))<=safetime(2)
                    % If is 1 means that the code founded the location of one
                    % peak of the cell in the period created, so this means the
                    % cell is active at the same time (+ cts seconds) that
                    % swimming is ocurring


                    RESPONSE.(myfield1).(myfield2)(j) = 1;
                    %  we create this structure that will tell us
                    % how many times a swimming precedes a peak. For
                    % example if two swimming episodes preceded the
                    % third peak it will be [3,3]

                  

                    % THIS is for the WINDOW
                    if window(1)>safetime(end) || window(end)<safetime(1) % when there is not overlap btw
                        % the vectors means that there is not timepoints overlaping
                        % so the previous to the calcium activity there is not
                        % swimming response
                    else
                         % elements in the swimming episode bigger then the window 1
                        w1 = (CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield2)))>window(1);

                         % elements in the swimming episode smaller then the window 2
                        w2 = (CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield2)))<window(2); 
                        
                        % this are only the timepoints of the swimming included in the window
                        % this are only the VALUES of the swimming included in the window
                        WINDOW_X.(myfield1).(myfield2) = (CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield2))).*w1.*w2;  
                        WINDOW_Y.(myfield1).(myfield2) = (CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2)).*w1.*w2;
                        clear w1 w2

                        %%%%Ploting for visualization
                        yyaxis left
                        hold on
                        safetimevector = linspace(safetime(1),safetime(2),100);
                        A = CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2)(1).*ones(1,length(safetimevector));
                        scatter(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield2)),CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield2),120,"red",'filled')
                        scatter(safetimevector,A,'r')
                        if WINDOW_X.(myfield1).(myfield2)==0
                        else
                            xline(nonzeros(WINDOW_X.(myfield1).(myfield2)),'r','LineWidth',4)
                        end
                        yyaxis right
                        hold on
                        window_vector = linspace(window(1),window(2),100);

                        C = CALCIUM.AMPLITUDE_Y.(myfield1)(j).*ones(1,length(window_vector));
                        scatter(window_vector,C,'m')
                        xline(CALCIUM.AMPLITUDE_X.(myfield1)(j),'M','LineWidth',4)


                    end



                else % Means that there is not ovelaping between the peak and the swimming
                    RESPONSE.(myfield1).(myfield2)(j) = 0; 

                end % overlaping
                 Onset
            end % current cswim episode
        end % Current peak
    end % Current cell
    %
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
        char(CALCIUM.list(nfish,1)),strcat('Windows',myfield1, char(CALCIUM.list(nfish,1)),'.fig')))
%     uiwait
    close all
end

% Finally saving all variables inside CALCIUM strucuture
CALCIUM.DELAY = DELAY;
CALCIUM.CaMovementSeconds = CA_movement;
CALCIUM.RESPONSE = RESPONSE;
CALCIUM.WINDOW_X = WINDOW_X;
CALCIUM.WINDOW_Y = WINDOW_Y;
CALCIUMimg('save',CALCIUM,[],list,nfish);


% TO SUM UP: this code is, taking one cell, taking each of the swimming
% episodes (taking the time when this this happening). Then we will take
% one this times and we will check if is inside the times where any episode
% of swimming is ocurring. if is inside it will put 1 otherwise 0. then it
% will take the next timepoint of the calcium response and it will check if
% is inside any episode of swimming. to not overwrite the previous vector
% created we will keep the ones. at the end, for each cell you will have a
%vector with length equal to the number of swimming episodes, and in each
%index will be a 1 if that cell responded to that specific swimming
%episode.



%
%
% %% CREATING DATA SHEET
%
% %The first row will be the name of the cell
% writecell(name_perloc',fullfile(path.data...
%     ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range','A1')
% % The second row will be the position normalize
% writematrix(CALCIUM.Normalized_x_distance,fullfile(path.data...
%     ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range','A2')
%
% sum1 = 0;
% for i = index_order
%     temp1 = fieldnames(WINDOW_Y);
%     myfield1 = char(temp1(i));
%     sum1 = length(i) + sum1; % for natural counting (1,2,3,...)
%     longitud = 0;
%     sum = 0;
%     for h=1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
%         temp2 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
%         myfield2 = char(temp2(h)); % extracting the name of each episode
%
%
%         if isfield(WINDOW_Y.(myfield1),(myfield2))
%             T([1+sum:length(WINDOW_Y.(myfield1).(myfield2))+sum],1) =...
%                 WINDOW_Y.(myfield1).(myfield2)
%             sum = sum + length(WINDOW_Y.(myfield1).(myfield2))
%         end
%     end
%     T = T(find(T)); % Removinf the zeros
%     writematrix(T,fullfile(path.data...
%         ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),'Sheet',1,'Range',...
%         strcat(char('A'+(sum1-1)),'3'))
%     clear longitud sum T
%
% end
%
%
% % Connect to Excel
% Excel = actxserver('excel.application');
% % Get Workbook object
% WB = Excel.Workbooks.Open(fullfile(path.data...
%     ,'dataCALCIUM',list(nfish),'nMLF.xlsx'),0,false);
%
% % Set the color of cell "A1" of Sheet 1 to RED
% WB.Worksheets.Item(1).Range(strcat('A2:',strcat(char('A'+(length(name_perloc)-1)),'2'))).Interior.ColorIndex = 3;
% % Save Workbook
% WB.Save();
% % Close Workbook
% WB.Close();
% % Quit Excel
% Excel.Quit();
