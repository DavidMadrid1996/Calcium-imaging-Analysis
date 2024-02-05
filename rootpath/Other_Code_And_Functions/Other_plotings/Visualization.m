%% 3D RECONSTRUCTION PER TRIAL
    user_settings
    list

    nfish =4; % Number of the fish to load (this is the trial for the same fish)

    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;
     wave = squeeze(CALCIUM.SumAng3D); 
%      plot(squeeze(CALCIUM.AngHeadTail3D)); 


 % BINARY ACTIVATION PLOT. IS THE CELL ACTIVE DURING THE EPISODE OR NOT? 
% 
%  scatter3(CALCIUM.ZAbsoluteCoordinates(:,1),CALCIUM.ZAbsoluteCoordinates(:,2),CALCIUM.ZAbsoluteCoordinates(:,2),'k'), hold on % Coordinates of the lateral reference
% % plot3(CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,1), CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,2),CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,3), 'color', 'k', 'LineWidth', 1); hold on

 TEMP = CALCIUM.BASELINE_TO_PEAK; 
temp1 = fieldnames(TEMP);


numEpisodeSwimming = 1; 
   
for i = 1:size(fieldnames(TEMP),1)
    myfield1 = char(temp1(i));

    if CALCIUM.RESPONSE.(myfield1)(numEpisodeSwimming)==0
    scatter(CALCIUM.ZAbsoluteCoordinates(i,1),CALCIUM.ZAbsoluteCoordinates(i,2),30,'r','filled'), hold on % Coordinates of the lateral reference
 
    else
        scatter(CALCIUM.ZAbsoluteCoordinates(i,1),CALCIUM.ZAbsoluteCoordinates(i,2),30,'g','filled'), hold on % 
    end 

end 



% Select trial

 % NON BINARY ACTIVATION PLOT. HOW MUCH THE CELL IS ACTIVE DURING THE EPISODE? 
clear col
   user_settings



    [CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
    CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);
    sr=1/CALCIUMroiTS.deeplabcut.sr;
     wave = squeeze(CALCIUM.SumAng3D); 
%      plot(squeeze(CALCIUM.AngHeadTail3D)); 

 TEMP = CALCIUM.BASELINE_TO_PEAK; 
temp1 = fieldnames(TEMP);
         
for i = 1:size(fieldnames(TEMP),1)
    myfield1 = char(temp1(i));
    col(i) = TEMP.(myfield1)(1)
        
end 

scatter(CALCIUM.ZAbsoluteCoordinates(:,1),CALCIUM.ZAbsoluteCoordinates(:,2),30,col,'filled')
caxis([-50 1000])



%%%%%%%%%% 

    for i = 1:length(CALCIUM.Latreference.manual_poly) % carefull now
        % is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.Latreference.manual_poly{i}; % will contain the x and y values of all the
        %point that contains the midline ROI. This temp variable will be
        %use later. dont erase it.
    end
    scatter3(temp(:,1),temp(:,2),CALCIUM.Zposition.CurrentRef.*ones(1,length(temp(:,1)))), hold on


scatter3([Coordinates(:,2);temp(:,1)],[Coordinates(:,3);temp(:,2)],[Coordinates(:,4);CALCIUM.Zposition.CurrentRef.*ones(1,length(temp(:,1)))'],[],'filled')

%% HISTOGRAM ROSTRO CAUDAL AXIS

Type = table2cell(ANIMAL(:,6));
for i = 1:length(Type)
    A(i) = isequal(Type{i},'SW'); % which are swimming
end

Change = table2array(ANIMAL(:,4)); % change to 0 or 1
for i = 1:length(Change)
    if Change(i) == 0
        Change2(i) = 0
    else
        Change2(i) = 1
    end
end

Z = table2array(ANIMAL(:,3));
c = sort(unique(Z));

for i=1:length(c)
    Zcoordinates (:,i) = c(i)==Z; % index of diferents z cordinates
end

col = ['b','g','r']; 
for l = 1:size(Zcoordinates ,2)

    clear temp
    temp = table2array(ANIMAL(:,2)).*Change2'.*A'.* Zcoordinates (:,l); % Y coordinates of the cells that are RESPONDING, SW, 
    temp = nonzeros(temp);
    nbins = 10;

    %%%%%
     scatter(temp, nonzeros(table2array(ANIMAL(:,5)).*Change2'.*A'.* Zcoordinates (:,l)),col(l)), hold on
    %%%%



    vec =  0:1/nbins:1;
    CountingCells = NaN(1,size(vec,2)-1);
    ActiveCells = NaN(1,size(vec,2)-1);
    for i = 1:length(vec)-1
        temp2 = find((table2array(ANIMAL(:,2))>vec(i)).*(table2array(ANIMAL(:,2))<vec(i+1)));
        temp3 = find((temp > vec(i)).*(temp < vec(i+1)));
        CountingCells(i) = length(temp2); % Number of cell in each subsegment (bin)
        ActiveCells (i) = length(temp3); % Number of active cells in each subsegment (bin)
    end

    PerActiveCell = ActiveCells./CountingCells;
    % PerActiveCell = medfilt1(PerActiveCell);
%     area(PerActiveCell,'FaceAlpha',0.2,'FaceColor','r'), hold on
    % bar(PerActiveCell,'FaceAlpha',0.5,'FaceColor','b')
    plot(PerActiveCell,col(l)), hold on
    % bar(CountingCells,'FaceAlpha',0.5), hold on
    % bar(ActiveCells,'FaceAlpha',0.5), hold on
    set(gca,'xtick',[1 length(CountingCells)],'xticklabel',({'SP/MEDULLA BOUNDARY','MEDULA/CEREBLLUM BOUNDARY'}))

 

    ylabel('% Active Cells')
    title(table2cell(CALCIUM.DATA(1,5)))


    % This plot is not taking in considerations de number of cell that are in
    % each subsegment of th HB. Is possible that some subsegments have more
    % cells so the absolute number of cells active will be bigger and other
    % subsegment there are few cells active but also few cells in general.
    subplot(122)
    plot(CountingCells); hold on
    xlabel('0 SPINAL/HB BOUNDARY ------ 1 MAUTHNER BODIES')
    title('Distribution of CHX10+ Cells in HB')


end
legend('Dorsal','middle','Ventral')
end
