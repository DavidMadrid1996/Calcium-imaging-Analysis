
%% Creating a list of parameters for the PCA.

%{

------------------------------------------------------------

** Branching 

-----> Bilateral (if BilateralIndex > BilateralIndex_Threshold)

-------------> Swimming (if abs(Amplitude) < Amplitude_Threshold && if Mean Frequency > Frequency_Threshold) % low amplitude and high frequency 

----------------------------> Fast Swimming (if Mean Frequency > Frequency_Threshold)

----------------------------> Slow Swimming (if )

-------------> Strugle (if abs(Amplitude) > Amplitude_Threshold && if Mean Frequency < Frequency_Threshold) % high amplitude and low frequency 
        


-----> Unilateral (if BilateralIndex < BilateralIndex_Threshold)

 

------------------------------------------------------------

** PCA: Analysis

                @ Matrix Dimensions: 

             --> each ROW is one behaviour 

             --> each COLUMN is one parameter


               @ List of Parameters: 

               1. Max Freq
               2. Mean Freq
               3. Std Freq
               4. Max Amplitude
               5. Mean Amplitude
               6. Std Amplitude
               7. Number of Cycles
               8. Bilateral index 
%}

clear
user_settings
S_11_Group_analysis_files
close all

%% Threshold for Braching

% BilateralIndex_Threshold = 0.6 ; % Value that goes from 0 to 1
% Amplitude_Threshold = 60 ; % Degrees
% Frequency_Threshold = 20 ; % Hz
% Fast_Frequency_Threshold = ; %
% Slow_Frequency_Threshold = ; %



for i = 1:size(Animal_List_unique,1)
    load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
        Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
        Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)
    myfield1 = char(Animal_List_unique(i,3)); % extracting the fieldnames of each cell of the loop



    % Manual lablelling of episodes (Branching)
    Manual_Label = [1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3];  % Manual label of the behaviour to color code the points

    % Next will be the billateral index, also in the bilaterallity
    % measurement in the script: "S_17_Group_analysis_Swimming"
    Variable = CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1); 
    Variable2 = CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1); 
    BilateralIndex= Bilateral_Index(Variable,Variable2); 
% 
%     if BilateralIndex > BilateralIndex_Threshold % The behaviour will be consider bilateral
%         if 
%         else
%         end 
% 
%     elseif
%     elseif
%     end











    % List of parameters
    
    if isempty(CALCIUM.FREQ_SWIM_Y.(myfield1)) % Because turning will have an empty frequency
    Parameters_list(1,i) = 0; % Max Freq
    Parameters_list(2,i) = 0; % Mean Freq
    Parameters_list(3,i) = 0; % Std Freq
    else 
    Parameters_list(1,i) = max(CALCIUM.FREQ_SWIM_Y.(myfield1)); % Max Freq
    Parameters_list(2,i) = mean(CALCIUM.FREQ_SWIM_Y.(myfield1)); % Mean Freq
    Parameters_list(3,i) = std(CALCIUM.FREQ_SWIM_Y.(myfield1)); % Std Freq
    end 
    Parameters_list(4,i) = max(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))); % Max Amplitude
    Parameters_list(5,i) = mean(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))); % Mean Amplitude
    Parameters_list(6,i) = std(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))); % Std Amplitude
    Parameters_list(7,i) = length(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)); % Number of Cycles

    Parameters_list(8,i) = BilateralIndex; % Bilateral index

    
end

% Transposing and PCA
Parameters_list =  Parameters_list'; %figure, imagesc(Parameters_list)
[u,z] = pca(Parameters_list,'NumComponents',3); 

% Ploting PCA 3D
figure
scatter3(z(:,1),z(:,2),z(:,3),40,Manual_Label,'filled'); 
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

% Ploting PCA 2D
figure
scatter(z(:,1),z(:,2),40,Manual_Label,'filled'); 
xlabel('PC1'); ylabel('PC2');





