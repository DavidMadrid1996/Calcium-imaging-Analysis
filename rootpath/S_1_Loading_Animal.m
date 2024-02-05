%% David Madrid . Last Rev 13/09/2023
clear
close
user_settings

%% Select animal
% List with the numbers of the trial
temp_table(:,1) = table(list); % First column of the table --> list
temp_vector = 1:length(list); % Second column of the table --> 1:length(list)
temp_table(:,2) = table(temp_vector');
temp_table.Properties.VariableNames = {'List','nfish'};
disp(temp_table) % Display the table as info

prompt = 'Select number of animal from the list that you want to analyze: ';
str = str2double(input(prompt,'s'));
if isnumeric(str)
    nfish = str;
    display(str)% Number of the fish to load (this is the trial for the same fish)
    pause(1)
else
    error('You must introduced a number')
end

[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
sr=1/CALCIUMroiTS.deeplabcut.sr;
wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D);


%% Checking what is already analyzed

if isfield(CALCIUM, 'EPISODE_Y_BASELINE_AMP_SIGN') % Swimming
    disp('Swimming is analyzed')
else
    disp('Swimming not analyzed')
end 

if isfield(CALCIUM, 'BASELINE_TO_PEAK') % Calcium signal
     disp('Calcium is analyzed')
else
    disp('Calcium not analyzed')
end 

if isfield(CALCIUM, 'WINDOW_Y') % Window
     disp('Window is analyzed')
else
    disp('Window not analyzed')
end 
CALCIUM.path = path;
CALCIUMimg('save',CALCIUM,[],list,nfish);