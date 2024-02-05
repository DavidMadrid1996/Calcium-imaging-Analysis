%% RASTER. 
%
% For a specific colormap.
% For finding a specific color go to https://htmlcolorcodes.com/es/

% For Bilateral 
color1 = [0 0.4470 0.7410]; 
color2 = [0.9290 0.6940 0.1250]; 

lengthH = 1000;
neutral = [243, 220, 218]./255;
blue = color1;
colors_p = [linspace(neutral(1),blue(1),lengthH)', linspace(neutral(2),blue(2),lengthH)',...
    linspace(neutral(3),blue(3),lengthH)'];

%% Creating group 
load('Selection_Long_delay_bilateral.mat'); % Selection_Long_delay BILATERAL
load('Selection_Short_delay_bilateral.mat'); % Selection_Short_delay BILATERAL
load('Selection_Long_delay_turn_left.mat'); % Selection_Long_delay turn left
load('Selection_Short_delay_turn_left.mat'); % Selection_Short_delay turn left
load('Selection_Long_delay_turn_right.mat'); % Selection_Long_delay turn left
load('Selection_Short_delay_turn_right.mat'); % Selection_Short_delay turn left

% 
Group_Data = [Data_ShortDelay(Selection_Short_delay,:);Data_LongDelay(Selection_Long_delay,:)]; %;Data_NonActive
Group_coordinates = [ShortDelay_COORDINATES(Selection_Short_delay,:);LongDelay_COORDINATES(Selection_Long_delay,:)]; %;NonActive_COORDINATES
Group_names = [Name__ShortDelay(Selection_Short_delay,:); Name__LongDelay(Selection_Long_delay,:)]; %; Name__NonActive

% Group_Data = [Data_ShortDelay(Selection_Short_delay,:);Data_LongDelay(Selection_Long_delay,:) ; Data_NonActive([1 2 3],:)]; 

%Turn left
Group_Data = [Data_ShortDelay(Selection_Short_delay_turn_left,:);Data_LongDelay(Selection_Long_delay_turn_left,:);Data_NonActive]; %;Data_NonActive %;Data_LongDelay
Group_coordinates = [ShortDelay_COORDINATES(Selection_Short_delay_turn_left,:);LongDelay_COORDINATES(Selection_Long_delay_turn_left,:);NonActive_COORDINATES]; %;NonActive_COORDINATES %;LongDelay_COORDINATES
Group_names = [Name__ShortDelay(Selection_Short_delay_turn_left,:); Name__LongDelay(Selection_Long_delay_turn_left,:);Name__NonActive]; %; Name__NonActive % ; Name__LongDelay
Group_time = [Time_ShortDelay(Selection_Short_delay_turn_left,:); Time_LongDelay(Selection_Long_delay_turn_left,:)]; 



%Turn right
Group_Data = [Data_ShortDelay(Selection_Short_delay_turn_right,:);Data_LongDelay(Selection_Long_delay_turn_right,:);Data_NonActive]; %;Data_NonActive %;Data_LongDelay 
Group_coordinates = [ShortDelay_COORDINATES(Selection_Short_delay_turn_right,:);LongDelay_COORDINATES(Selection_Long_delay_turn_right,:);NonActive_COORDINATES]; %;NonActive_COORDINATES %;LongDelay_COORDINATES 
Group_names = [Name__ShortDelay(Selection_Short_delay_turn_right,:); Name__LongDelay(Selection_Long_delay_turn_right,:);Name__NonActive]; %; Name__NonActive % ; Name__LongDelay ;Name__NonActive
Group_time = [Time_ShortDelay(Selection_Short_delay_turn_right,:); Time_LongDelay(Selection_Long_delay_turn_right,:)]; 


%plot(Group_Data(4,:))
%% Normalization

close all
clear tempfilter
temp = -Group_Data;
prompt = 'Do you want to normalize it?';
str = input(prompt,'s');


% 24
%plot(-Group_Data(19,:))
%plot(Data_NonActive(:,1))
suma = 1; 
for i =20:40
    
    if str=='Y'

        if max(temp(i,:))<0

            tempfilter(suma,:) = medfilt1(temp(i,:)./min(temp(i,:)),8);
        else
            tempfilter(suma,:) = medfilt1(temp(i,:)./max(temp(i,:)),8);
        end
    elseif str=='N'
        tempfilter(suma,:) = medfilt1(temp(i,:));
    end
   plot(tempfilter(suma,:),'G'),hold on
    %         uiwait
    %     tempfilter = temp
    % plot(tempfilter(i,:)),hold on
%     pause(2)
    % disp(i)
    suma = 1 + suma; 
end
 plot(mean(tempfilter),'G','LineWidth',5)
% % 
% clear sum newGroup_Data
% vec = [8 1 9 10 11 2 3 4 5 6 7 17 16 15 13 14 18]; 
% vec = [18 14 13 15 16 17 7 6 5 4 3 2 11 10 9 1 8]
% sum = 1
% for i = vec
% newGroup_Data(sum,:) = Group_Data(i,:)
% sum = sum + 1 ;
% end 
% 
% Group_Data = newGroup_Data; 
% imagesc(tempfilter)


clear temp

%% Raster

    % If medio lateral raster--> raster = 1
    % If caudal rostral raster--> raster = 2
    raster = 2; 


    Excuding_Neurons = []; 
 Test = Group_coordinates(:,raster);
% Test = Group_time; 
%     Test(Excuding_Neurons)=[]; 
    if raster == 1
        LeftRigthOrder = sort(Test); 
        Title = 'Medio lateral normalized'; 
    elseif raster == 2
    LeftRigthOrder = sort(Test,"descend"); 
    Title = 'Caudal to rostral normalized'; 
    end 

    count = 1; 
    clear LeftRigthOrderMatrix Num
    for i = 1:length(LeftRigthOrder)
        temp = find(LeftRigthOrder(i)==Test); 
        LeftRigthOrderMatrix(count,:) = tempfilter(temp(1),:);
        LeftRigthOrderMatrix_Name(count,:) = Group_names (temp(1),:);
            count = count + 1; 
        temp1 = find(LeftRigthOrder(i)==Test); 
        Num(i) = temp1(1); 
    end 
figure(101)
    Raster = imagesc(LeftRigthOrderMatrix); 
    %plot(LeftRigthOrderMatrix(2,:))
%     imagesc(LeftRigthOrderMatrix(100,:))
%    colormap(colors_p), colorbar
colormap('turbo')

    hold on
   caxis([0 1])
   colorbar
    yline(length(find( CALCIUM.ZLeftOrRigth<0))+0.5,'LineWidth',3)
    xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
    title(Title)
    yticks([1:length(LeftRigthOrder)])
    yticklabels(num2cell(Num) )



% %% Next plot will generate a raster plot with neurons order or mediolateral o rostrocaudal. The raster plot
% % will include all the positions and just eh position which have neurons in
% % this trials will be filled with data
% 
% clear
% clc
% user_settings
% list
% 
%     nfish =8; % Number of the fish to load (this is the trial for the same fish)
%     % If medio lateral raster--> raster = 1
%     % If caudal rostral raster--> raster = 2
%     raster = 2; 
% 
%     [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
%     CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
%     CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data'; 
%     wave = squeeze(CALCIUM.SumAng3D)'; 
%     if raster == 1
%         LeftRigthOrder = sort(CALCIUM.ZAbsoluteCoordinates(:,raster)); 
%         Title = 'Medio lateral normalized'; 
%     elseif raster == 2
%     LeftRigthOrder = sort(CALCIUM.ZAbsoluteCoordinates(:,raster),"descend"); 
%     Title = 'Caudal to rostral normalized'; 
%     %CALCIUM.ZAbsoluteCoordinates(:,2) = CALCIUM.ZAbsoluteGeneralRostralDistance' + (rand(15,1)./100000000)
%     end 
% 
% 
%     linspace(0,2,0.1)
%     count = 1; 
%     for i = 1:length(LeftRigthOrder)
%         LeftRigthOrderMatrix(count,:) =  CALCIUMroiTS.diff_perc03.data(find(LeftRigthOrder(i)==CALCIUM.ZAbsoluteCoordinates(:,raster)),:);  %- min(CALCIUMroiTS.diff_perc03.data(find(LeftRigthOrder(i)==CALCIUM.ZAbsoluteCoordinates(:,raster)),:)))...
%             %./max(CALCIUMroiTS.diff_perc03.data(find(LeftRigthOrder(i)==CALCIUM.ZAbsoluteCoordinates(:,raster)),:)- min(CALCIUMroiTS.diff_perc03.data(find(LeftRigthOrder(i)==CALCIUM.ZAbsoluteCoordinates(:,raster)),:))); 
%         count = count + 1; 
%         Num(i) = find(LeftRigthOrder(i)==CALCIUM.ZAbsoluteCoordinates(:,raster)); 
% 
%     end 
% 
%     Raster = imagesc(-LeftRigthOrderMatrix), colormap('turbo'), hold on
%    caxis([-50 1000])
%    colorbar
%     yline(length(find( CALCIUM.ZLeftOrRigth<0))+0.5,'LineWidth',3)
%     xline(nonzeros(unique(CALCIUM.CaMovementSeconds))./CALCIUM.srate,'r','LineWidth',3)
%     title(Title)
%     yticks([1:length(LeftRigthOrder)])
%     yticklabels(num2cell(Num) )
% 
% 
% 
%    
% 
% %     for i = 1:length(fieldnames(CALCIUM.AMPLITUDE_X))
% %         temp1 = fieldnames(CALCIUM.AMPLITUDE_X);
% %         myfield1 = char(temp1(i))
% %         CALCIUM.AMPLITUDE_X.(myfield1)
% %     end 
% % 
% %     LeftOrRigth
% % 
% 
% 
% 
%     % Ploting rescaling time 
% for i = 1:length(CALCIUMroiTS.roi.labels)
%     yyaxis left
%     r = -waves(:,i,1);
% figure(i)
% TimeCameraRescale = (1:length(x))./length(x);
% plot(TimeCameraRescale,y,'g','LineWidth',3); hold on
% yyaxis right
% TimeCalciumRescale  = (1:length(timev))./length(timev);
% plot(TimeCalciumRescale,r,'b','LineWidth',3); 
% 
% yyaxis left
% title(strcat('Calcium activity-Tail movement','-',CALCIUM.ref,'-',CALCIUMroiTS.roi.labels(i)));
% xlabel('Time [s]')
% ylabel('\Deltax Coordinates tail3 [px]')
% 
% yyaxis right
% ylabel('% \Delta F')
% saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),...
%     strcat(char(CALCIUMroiTS.roi.labels(i)),'_','TimeRescaled',char(CALCIUM.list(nfish,1)),'.fig')))
% end 
