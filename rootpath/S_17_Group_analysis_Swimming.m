%% This script is used to 

%% General view all recordings. This code will plot all the behavioural recordings in the Animal_List_unique variable 
clear
user_settings
S_11_Group_analysis_files
close all

for i = 1:size(Animal_List_unique,1)
    load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)

    figure(i)
    plot(CALCIUMroiTS.deeplabcut.time_scaled,wave),hold on
    yline(0,'LineWidth',2)


        myfield1 = char(Animal_List_unique(i,3)); % extracting the fieldnames of each cell of the loop
        plot(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),'LineWidth',2)
        scatter(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),30,'filled')
         text(CALCIUMroiTS.deeplabcut.time_scaled(CALCIUM.EPISODE_X.(myfield1)),CALCIUM.EPISODE_Y.(myfield1),string(CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.( myfield1)))

    title(strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),Animal_List_unique(i,3)))
end


%% Amplitude
close all
S_11_Group_analysis_files
suma2 = 1; % This is for boxplot
suma3 = 0;  % This is for boxplot
Group = table; 

for i = 1:size(Animal_List_unique,1)
    load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)
    %% Variable
    Variable = CALCIUM.EPISODE_Y_BASELINE_AMP; 
    Variable2 = CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN; 



        myfield1 = char(Animal_List_unique(i,3)); 
        CellName {1,i} = char(strcat(strcat(Animal_List_unique(i,1),Animal_List_unique(i,2),'_',myfield1)));
        Group(1:length(Variable.(myfield1)),i) = array2table(abs(Variable.(myfield1)).*Variable2.(myfield1)'); 

        suma2 = length(Variable.(myfield1)); 
        DataBoxplot_label([1+suma3:suma2+suma3],1) = string(repmat(CellName(1,i), 1,suma2));
        DataBoxplot_Data([1+suma3:suma2+suma3],1) = abs(Variable.(myfield1)).*Variable2.(myfield1)'; 
        suma3 = suma2 + suma3; 




end

Group.Properties.VariableNames(1:length(CellName)) = CellName ;
figure(1) % Amplitudes including positive and negatives values
boxplot(DataBoxplot_Data,char(DataBoxplot_label))
yline(0)
ylabel('Amplitude [degrees]')
figure(2) % Amplitudes including just positive values
boxplot(abs(DataBoxplot_Data),char(DataBoxplot_label))
yline(0)
ylabel(' Absoluyte value Amplitude [degrees]')

clear DataBoxplot_Data DataBoxplot_label


%% Frequency
clear
close all
user_settings
S_11_Group_analysis_files

suma2 = 1; % This is for boxplot
suma3 = 0;  % This is for boxplot
Group = table; 
figure(2)

for i = 1:size(Animal_List_unique,1)

 load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)

    %% Variable
       Variable = CALCIUM.GENERAL_FREQ_Y; %CALCIUM.FREQ_SWIM_Y; 

        myfield1 = char(Animal_List_unique(i,3));
        CellName {1,i} = char(strcat(Animal_List_unique(i,1),Animal_List_unique(i,2),'_',myfield1)); 
        if isempty(Variable.(myfield1)')
            Group(1:1,i) = array2table(0);
        else
            Group(1:length(Variable.(myfield1)),i) = array2table(Variable.(myfield1)');
        end
         suma2 = length(Variable.(myfield1)); 
        DataBoxplot_label([1+suma3:suma2+suma3],1) = string(repmat(CellName(1,i), 1,suma2));
        DataBoxplot_Data([1+suma3:suma2+suma3],1) = Variable.(myfield1); 
        suma3 = suma2 + suma3; 





end
Group.Properties.VariableNames(1:length(CellName)) = CellName ;
boxplot(DataBoxplot_Data,char(DataBoxplot_label))
ylabel('Frequency [Hz]')
clear DataBoxplot_Data DataBoxplot_label

%% Number of cycles

S_11_Group_analysis_files
Group = table; 
figure(3)

for i = 1:size(Animal_List_unique,1)

 load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)

    %% Variable
    Variable = CALCIUM.EPISODE_Y_BASELINE_AMP;

        myfield1 = char(Animal_List_unique(i,3));
        CellName {1,i} = char(strcat(Animal_List_unique(i,1),Animal_List_unique(i,2),'_',myfield1)); 

            Group(1,i) = array2table(length(Variable.(myfield1)));

        
        DataBoxplot_label(i,1) = CellName(1,i);
        DataBoxplot_Data(i,1) = length(Variable.(myfield1)); 


end
Group.Properties.VariableNames(1:length(CellName)) = CellName ;
plot(1:size(DataBoxplot_Data,1),DataBoxplot_Data)
newStr = strrep(DataBoxplot_label,'_','.'); 
set(gca,'xticklabel',newStr);
ylabel('Number of episodes')



%% Bilaterality
clear
close all
clear BilateralIndex
user_settings
S_11_Group_analysis_files
Group = table; 
figure(6)

for i = 1:size(Animal_List_unique,1)

 load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)
    myfield1 = char(Animal_List_unique(i,3));

%%   Control of errors: CHECKING THAT THE SINGS OF THE AMPLITUDES ARE OKEY

% but carreful, some times if the value of CALCIUM.EPISODE_Y is too close
% to 0 then when the normalization is applied the value of
% CALCIUM.EPISODE_Y_BASELINE_AMP can change signs, not because it is wrong,
% just because mabe the baseline was slightly up or down. 
for ii = 1:length(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))
    if CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)(ii)<0 && CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1)(ii)<0 %&& CALCIUM.EPISODE_Y.(myfield1)(ii)<0 
    elseif CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)(ii)>0  && CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1)(ii)>0 %&& CALCIUM.EPISODE_Y.(myfield1)(ii)>0
    else 
      
        error
    end 
end 


    %% Variable

    Variable = CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1); 
    Variable2 = CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1); 

    BilateralIndex(i) = Bilateral_Index(Variable,Variable2); 
   

       
       CellName {1,i} = char(strcat(Animal_List_unique(i,1),Animal_List_unique(i,2),'_',myfield1));




end
scatter(1:length(BilateralIndex),BilateralIndex)
ylabel('Bilaterallity index')


%% Struggle index

clear
close all
clear StruggleIndex
user_settings
S_11_Group_analysis_files
Group = table; 
figure(7)

for i = 24:size(Animal_List_unique,1)

 load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)
    myfield1 = char(Animal_List_unique(i,3));
   
    %% Variable

    TraceAngle = squeeze(CALCIUM.Ang3D); 

    matrix = TraceAngle(1:4,CALCIUM.EPISODE_X.(myfield1)(1)-50 : CALCIUM.EPISODE_X.(myfield1)(end)+50 )
    trace1 = TraceAngle(3,CALCIUM.EPISODE_X.(myfield1)(1)-50 : CALCIUM.EPISODE_X.(myfield1)(end)+50 )
    trace4 = TraceAngle(4,CALCIUM.EPISODE_X.(myfield1)(1)-50 : CALCIUM.EPISODE_X.(myfield1)(end)+50 )

     close all
     imagesc(matrix)
     plot(trace1),hold on
     plot(trace4)
     legend('Rostral','Caudal')

     plot(TraceAngle(1:4,CALCIUM.EPISODE_X.(myfield1)(1)-50 : CALCIUM.EPISODE_X.(myfield1)(end)+50 ))
       
    CellName {1,i} = char(strcat(Animal_List_unique(i,1),Animal_List_unique(i,2),'_',myfield1));




end
scatter(1:length(BilateralIndex),BilateralIndex)
ylabel('Bilaterallity index')




%% Swimming struggle comparison

S_11_Group_analysis_files
close all
figure(1)

for i = 1:size(Animal_List_unique,1)
    load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)

%     g{i} = CALCIUM.LABELLINGBH

myfield1 = char(Animal_List_unique(i,3));

scatter(CALCIUM.FREQ_SWIM_Y.(myfield1),CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)(2:end),40,'MarkerEdgeColor','b',...
              'MarkerFaceColor','b',...
              'LineWidth',1.5),hold on
xlim([0 50])
ylim([0 80])
xlabel('Freq [Hz]')
ylabel('Amplitude [º]')
end



%% Creating a continous vector with all behaviour recordings (Not animal dependent, which means it will
% not distinguis btw animals, it will take all behaviours ever recorded and it will put them together)

S_11_Group_analysis_files
close all

     suma = 0; 
     HowMuch = 20; % How much extra before and after the first and last peak are we taking
     zeropading = 10; % How uch will zeropad before and after (not including HowMuch)

for i = 1:size(Animal_List_unique,1)
    load(char(fullfile(ValidPathdataCalcium(i,:),strcat('CALCIUM_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = load(char(fullfile(ValidPathdatawaves(i,:),strcat('CALCIUMRoiTS_',strcat(Animal_List_unique(i,1),...
     Animal_List_unique(i,2),'.mat')))));
    [CALCIUMroiTS] = CALCIUMroiTS.VSDroiTS;
    sr=1/CALCIUMroiTS.deeplabcut.sr;
    wave = squeeze(CALCIUM.AngHeadTail3D);  %squeeze(CALCIUM.SumAng3D)
    

     myfield = char(Animal_List_unique(i,3));
     Selection =  ((CALCIUM.EPISODE_X.(myfield)(1))-HowMuch):((CALCIUM.EPISODE_X.(myfield)(end))+HowMuch);
     Trace = wave((Selection(1)-HowMuch):(Selection(end)+HowMuch)); plot(Trace)
     TraceFiltered = medfilt1(Trace,5); plot(TraceFiltered)
     Zero_TraceFiltered = [zeros(1,zeropading) TraceFiltered' zeros(1,zeropading)]; 
     continuousTrace(1+suma:length(Zero_TraceFiltered)+suma) = Zero_TraceFiltered; 
     suma = suma + length(Zero_TraceFiltered); 

end

plot(continuousTrace)
