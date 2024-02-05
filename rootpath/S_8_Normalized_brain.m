%% Normalization will be

% ROSTRO-CAUDAL--> (0)JUNCTION HINDBRAIN-SPINAL CORD --- (1)MAUTHNER AXON CROSS
% MEDIO-LATERAL--> (0)MIDLINE --- (1) THE WIDEST POINT IN THE RIGHT CAUDAL HINDBRAIN

if exist('WHOLE_ZSTACK', 'dir')

else
    mkdir (fullfile(path.data,'WHOLE_ZSTACK',CALCIUM.ref)); % Creating the directory for dataCALCIUM

    prompt = {'Where the mth axon cross','Juntion spinal cord-hindbrain','Widest point caudal hindbrain','midline'};

    % Here you have to select the medial point

    dlgtitle = 'Absolute Z stack normalization';
    dims = [1 35];

    definput ={sprintf('%.0f',0),sprintf('%.0f',0),sprintf('%.0f',0),sprintf('%.0f',0)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);


    % Upload the zstack
    data2 = bfopen(char(CALCIUM.zstack_directory_loop));

    for j = 1:length(data2{1,1})
        movies4D2 (:,:,j) = data2{1,1}{j,1};
    end
    movies4D2 = im2double(movies4D2);
    zstack = movies4D2;

    CALCIUM.zstackMost.Rostral =  zstack(:,:,str2double(answer(1))) ;
    CALCIUM.zstackMost.Caudal = zstack(:,:,str2double(answer(2))) ;
    CALCIUM.zstackMost.Lateral =  zstack(:,:,str2double(answer(3))) ;
    CALCIUM.zstackMost.Medial =  zstack(:,:,str2double(answer(4))) ;

    CALCIUMimg('save',CALCIUM,[],list,nfish);

end

% For this we dont need to draw a complete ROI, we just need the most
% rostral, cadudal and lateral point of the zstack

%% For the most Lateral
imagesc(histeq(CALCIUM.zstackMost.Lateral)); colormap('bone')
title('Mark the lateral reference/s')
disp('Mark the lateralreference/s.')
clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
CALCIUM.zstackMost.LateralCoord(1) = num_neurons_x;
CALCIUM.zstackMost.LateralCoord(2) = num_neurons_y;
close all
clear num_neurons_x num_neurons_y


%% For the most medial

clear temp2 temp3 temp temp8 temporal

imagesc(histeq(CALCIUM.zstackMost.Medial)); colormap('bone'),hold on
scatter(CALCIUM.zstackMost.LateralCoord(1),CALCIUM.zstackMost.LateralCoord(2),100,'filled','red')
title('Mark the Medial reference/s')
disp('Mark the Medial reference/s.')
clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
CALCIUM.zstackMost.MedialCoord(1) = num_neurons_x;
CALCIUM.zstackMost.MedialCoord(2) = num_neurons_y;
close all
clear num_neurons_x num_neurons_y


temp2 = sqrt((CALCIUM.zstackMost.MedialCoord(1)-CALCIUM.zstackMost.LateralCoord(1))^2 +...
    (CALCIUM.zstackMost.MedialCoord(2)-CALCIUM.zstackMost.LateralCoord(2))^2);

ZAbsoluteMedial2LateralDistance = temp2 ; % Absolute distance btw the most lateral and the most medial borders in the hind brain

ZAbsoluteGeneralLateralDistance = (CALCIUM.ZDistance2Midline./ZAbsoluteMedial2LateralDistance).*CALCIUM.ZLeftOrRigth; % Normalize with respect to the cube
% where 0 is the midline and 1 is the most lateral border of the whole
% hindbrain.

ZAbsoluteCoordinates(:,1) = ZAbsoluteGeneralLateralDistance;







%% For the most caudal
imagesc(histeq(CALCIUM.zstackMost.Caudal)); colormap('bone')
title('Mark the Caudal reference/s')
disp('Mark the Caudal reference/s.')
clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
CALCIUM.zstackMost.CaudalCoord(1) =  num_neurons_x;
CALCIUM.zstackMost.CaudalCoord(2) =  num_neurons_y;
close all
clear num_neurons_x num_neurons_y


% For this we calculate the distance from each midpoint of each ith neuron
% to the most caudal edge of the hindbrain
for i = 1:size(CALCIUM.ZMedialReference ,1)
    temp2(i) = sqrt(( CALCIUM.ZMedialReference (i,1)-CALCIUM.zstackMost.CaudalCoord(1))^2 + ...
        (CALCIUM.ZMedialReference  (i,2)-CALCIUM.zstackMost.CaudalCoord(2))^2); % Absolute Distance for each midpoint to the most caudal
    % border
end
ZAbsoluteDistance2Caudal =temp2;


clear temp2 temp3 temp temp8 temporal


%% For the most rostral
imagesc(histeq(CALCIUM.zstackMost.Rostral)); colormap('bone')
title('Mark the Rostral reference/s')
disp('Mark the Rostral reference/s.')
clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
CALCIUM.zstackMost.RostralCoord(1) = num_neurons_x;
CALCIUM.zstackMost.RostralCoord(2) = num_neurons_y;
close all
clear num_neurons_x num_neurons_y


temp2= sqrt((CALCIUM.zstackMost.CaudalCoord(1)-CALCIUM.zstackMost.RostralCoord(1))^2 + ...
    (CALCIUM.zstackMost.CaudalCoord(2)-CALCIUM.zstackMost.RostralCoord(2))^2); % Absolute distance from the most caudal and
% most rostral borders of the hinbrain

ZAbsoluteCaudal2RostralDistance =temp2;


ZAbsoluteGeneralRostralDistance = ZAbsoluteDistance2Caudal./ZAbsoluteCaudal2RostralDistance;


ZAbsoluteCoordinates(:,2) = ZAbsoluteGeneralRostralDistance;  % + (rand(size(ROI_point_pos,1),1)/1000000);
%CALCIUM.ZAbsoluteCoordinates(:,2) = CALCIUM.ZAbsoluteGeneralRostralDistance' + (rand(15,1)./100000000)
% that 2 neurons have the exact same value for the rostrocaudal. this is
% addind variability without affecting to the actual position beacuse is
% adding a decimal of 0.0000

% The dorso ventral normalization is already done
ZAbsoluteCoordinates(:,3) = CALCIUM.ZCoordinates(:,3);

% Saving
CALCIUM.ZAbsoluteCoordinates = ZAbsoluteCoordinates;
CALCIUM.ZAbsoluteCaudal2RostralDistance = ZAbsoluteCaudal2RostralDistance ;
CALCIUM.ZAbsoluteGeneralRostralDistance  = ZAbsoluteGeneralRostralDistance ;
CALCIUMimg('save',CALCIUM,[],list,nfish);



%% Final test ploting for the rostro caudal axis

%{

% Finding ZmidlineReference  for each corresponding ZlateralReference
clear temp2 temporal
    for i = 1:length(CALCIUM.ZLatreference.manual_poly{1,1}(:,1))
      temp2= sqrt((CALCIUM.ZMIDreference.manual_poly{1,1}(:,1)-CALCIUM.ZLatreference.manual_poly{1,1}(i,1)).^2 + ...
    (CALCIUM.ZMIDreference.manual_poly{1,1}(:,2)-CALCIUM.ZLatreference.manual_poly{1,1}(i,2)).^2);
    temporal =  CALCIUM.ZMIDreference.manual_poly{1,1}(temp2 == min(temp2),:); 
            ZAbsoluteMedialReferenceFORLATERAL(i,[1 2]) =temporal(1,:); % closest point in the midline to the ith neuron
            ZAbsoluteDistance2MidlineFORLATERAL(i) = (min(temp2));
        clear temp2 temporal
        ZLeftOrRigthFORLATERAL(i) = sign( ZAbsoluteMedialReferenceFORLATERAL(i,1)-CALCIUM.ZLatreference.manual_poly{1,1}(i,1))
%         figure(i)
%         scatter(CALCIUM.ZLatreference.manual_poly{1,1}(:,1),CALCIUM.ZLatreference.manual_poly{1,1}(:,2)), hold on
% scatter(CALCIUM.ZLatreference.manual_poly{1,1}(i,1),CALCIUM.ZLatreference.manual_poly{1,1}(i,2),'r')
% scatter(CALCIUM.ZMIDreference.manual_poly{1,1}(:,1),CALCIUM.ZMIDreference.manual_poly{1,1}(:,2),'k')
% scatter( ZMedialReferenceFORLATERAL(i,1), ZMedialReferenceFORLATERAL(i,2),'r')
    end 
 % Absolute distance btw the most lateral and the most medial borders in the hind brain

ZAbsoluteGeneralLateralDistanceFORLATERAL = ( ZAbsoluteDistance2MidlineFORLATERAL./ZAbsoluteMedial2LateralDistance).*ZLeftOrRigthFORLATERAL; % Normalize with respect to the cube
% where 0 is the midline and 1 is the most lateral border of the whole
% hindbrain.


% For the rostral
clear temp2
for i = 1:size( ZAbsoluteMedialReferenceFORLATERAL,1)
temp2(i) = sqrt(( ZAbsoluteMedialReferenceFORLATERAL (i,1)-CALCIUM.zstackMost.CaudalCoord(1))^2 + ...
    (ZAbsoluteMedialReferenceFORLATERAL(i,2)-CALCIUM.zstackMost.CaudalCoord(2))^2); % Absolute Distance for each midpoint to the most caudal
% border
end 
 ZAbsoluteDistance2CaudalFORLATERAL = temp2; 
 
ZAbsoluteGeneralRostralDistanceFORLATERAL = ZAbsoluteDistance2CaudalFORLATERAL./ZAbsoluteCaudal2RostralDistance;

scatter(ZAbsoluteGeneralLateralDistanceFORLATERAL,ZAbsoluteGeneralRostralDistanceFORLATERAL), hold on % Coordinates of the lateral reference
scatter(ZAbsoluteCoordinates(:,1),ZAbsoluteCoordinates(:,2)) % Coordinates of the neurons

%tO HAVE BETTER RESULT PUT MORE POINTS IN THE ZMIDLINE OTHERWISE IT WILL
%LOCK REALLY STRATIFICATE BECAUSE A LOT OF LATERAL REFERENCES WILL HAVE THE
%SAME MID LINE REFERNCE


CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,1) = ZAbsoluteGeneralLateralDistanceFORLATERAL;
CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,2) = ZAbsoluteGeneralRostralDistanceFORLATERAL; 
CALCIUM.ZAbsoluteCoordinatesFORLATERAL (:,3) = CALCIUM.ZAbsoluteCoordinates(:,3); 
CALCIUMimg('save',CALCIUM,[],list,nfish);





%% Automatic Finding reference with respect to the zstack

         clear temp2 temp3 temp temp8 temporal

 imagesc(CALCIUM.mean); colormap('bone'); hold on

for nroi = 1:size(CALCIUM.roi.manual_poly,1)
    coord = CALCIUM.roi.manual_poly{nroi};
    plot(coord(:,1), coord(:,2), 'LineWidth', 1); hold on
    text(mean(coord(:,1)), mean(coord(:,2)), CALCIUM.roi.labels(nroi),'Color','red','FontSize',14);
end

 Bright = 0.2; 
    figure(2)
  
%      imagesc(CALCIUM.zstackImage), colormap('bone'); % Use if you have the SNAP
%     imagesc(adapthisteq(CALCIUM.zstackImage+ Bright)), colormap('bone'); % 
    imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'), hold on
pause(0.2)

    prompt = 'Which neuron do you want as referene?';
    str = input(prompt,'s');
    
    close all
    figure(1)
            imagesc(CALCIUM.mean), colormap('bone');
%     imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone'); % Use if you have the SNAP
    hold on

    jjjj =  str2double(str); % # Number of the cell that will be used us a reference
    scatter(ROI_point_pos(jjjj,1),ROI_point_pos(jjjj,2),'w','filled')
    
    Bright = 0.2; 
    figure(2)
  
%      imagesc(CALCIUM.zstackImage), colormap('bone'); % Use if you have the SNAP
    imagesc(adapthisteq(CALCIUM.zstackImage+ Bright)), colormap('bone'); % Use if you have the SNAP
    hold on

    % Finding neuron of reference
    CALCIUMREFERENCENEURON.roi.labels = {'R1'};
    figure(2)
    clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
    [manual_poly, manual_mask] = roi_draw(adapthisteq(CALCIUM.zstackImage+ Bright),CALCIUMREFERENCENEURON.roi.labels,num_neurons_x,num_neurons_y); %in each column, one set of rois



    temp = manual_poly{1}; % for each x and y coordinate of one neuron
    xpos_REF = mean (temp(:,1)); % the first column is the x coordinates of each point
    % of one ROI of one neuron
    ypos_REF= mean (temp(:,2)); % the second column is the y coordinates of each point
    % of one ROI of one neuron
    clear temp

    % Each neuron will be treated like a individual point.
    REF_point_pos = [xpos_REF' ypos_REF'];  % Matrix containing the average x and y position for each neuron
    %   scatter( ROI_point_pos(:,1), -ROI_point_pos(:,2))



    %%%% Midreference

    for i = 1:length(CALCIUM.ZMIDreference.manual_poly) % carefull now
            % is made just for one ROI (MID REFERENCE), the midline
            temp = CALCIUM.ZMIDreference.manual_poly{i}; % will contain the x and y values of all the
            %point that contains the midline ROI. This temp variable will be
            %use later. dont erase it.
    end

     % Calculating the reference point in the midline for each cell
        % cell in a loop
        ZMedialReference = NaN(size(REF_point_pos,1),2);
        ZDistance2Midline = NaN(size(REF_point_pos,1),1);
        ZLeftOrRigth= NaN(size(REF_point_pos,1),1);
%         Ztemp2= NaN(1, size(temp,1));

        for i = 1:size(REF_point_pos,1)
            for ii = 1:size(temp,1)
                temp2(ii) = sqrt((temp(ii,1)-REF_point_pos(i,1))^2 + (temp(ii,2)-REF_point_pos(i,2))^2);
            end
            
            ZMedialReference(i,[1 2]) = temp(temp2 == min(temp2),:); % closest point in the midline to the ith neuron 
            ZDistance2Midline(i) = unique(min(temp2)); % Distance from closest point to the ith neuron

            % to figure out if the ith neuron is to the left or to the right we
            % will create a vector form the MedialReference to one with a
            % higher y value. this is creating a direction (up and maybe tilted if the brain is tilted)
            % After that we will create a vector for the MedialReference tot
            % the ith neuron. if the angle is created clock wise means that is
            % to the rigth, if is created counterclockwise is to the left.
            temp4 = find(temp(:,2)>ZMedialReference(i,2));
            if isempty(temp4)
                greater(i,[1,2]) = greater(i-1,[1,2]); 
            else
                greater(i,[1,2]) = temp(temp4(1),:);
            end 
            
            UpVector = [greater(i,1)-ZMedialReference(i,1) greater(i,2)-ZMedialReference(i,2)];
            IthVector = [REF_point_pos(i,1)-ZMedialReference(i,1) REF_point_pos(i,2)-ZMedialReference(i,2)];
            %figure, scatter(ROI_point_pos(i,1),ROI_point_pos(i,2)), hold on, scatter(MedialReference(i,1),MedialReference(i,2))
            %  figure,arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
            %  arrow([MedialReference(i,1),MedialReference(i,2)],[ROI_point_pos(i,1),ROI_point_pos(i,2)],1,30,60)

            temp3 = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));

            if temp3 > 0
                ZLeftOrRigth(i) = -1; % to the left of the midline
            elseif temp3 < 0

                ZLeftOrRigth(i) = 1;  % to the rigth of the midline
            else
                ZLeftOrRigth = NaN;
            end

           
        end

         clear temp2 temp3 temp temp8

        % Final Variables --> MedialReference, Distance2Midline, LeftOrRigth


         %%%% Lateralreference

       for i = 1:length(CALCIUM.ZLatreference.manual_poly) % carefull now
            % is made just for one ROI (MID REFERENCE), the midline
            temp = CALCIUM.ZLatreference.manual_poly{i}; % will contain the x and y values of all the
            %point that contains the midline ROI. This temp variable will be
            %use later. dont erase it.
        end

        % Calculating the reference point in the midline for each cell
        % cell in a loop
        close all
        ZLateralReference = NaN(size(REF_point_pos,1),2);
        ZDistance2Lateral = NaN(size(REF_point_pos,1),1);
        ZRostralReference = NaN(size(REF_point_pos,1),2);
        ZDistance2Rostral = NaN(size(REF_point_pos,1),1);
        ZCaudalReference = NaN(size(REF_point_pos,1),2);
        ZDistance2Caudal = NaN(size(REF_point_pos,1),1);



        for i = 1:size(REF_point_pos,1)

            UpVector = [greater(i,1)-ZMedialReference(i,1) greater(i,2)-ZMedialReference(i,2)];

            % For the lateral

            for ii = 1:size(temp,1)
                temp2= [temp(ii,1)-ZMedialReference(i,1) temp(ii,2)-ZMedialReference(i,2)];
                temp3(ii,1) = atan2d(UpVector(1)*temp2(2)-UpVector(2)*temp2(1),UpVector(1)*temp2(1)+UpVector(2)*temp2(2));
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
            end


            if ZLeftOrRigth (i) == 1 % is to the right angles will be negative
                temp4 = temp3 + 90;
                temp4(temp4<0) = max(temp4);
                ZLateralReference(i,[1 2]) = temp(temp4 == min(temp4),:);
                %         scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
                %         scatter(LateralReference(i,1),LateralReference(i,2))
                ZDistance2Lateral(i) = sqrt((ZLateralReference(i,1) -REF_point_pos(i,1))^2 + (ZLateralReference(i,2) -REF_point_pos(i,2))^2);
            elseif ZLeftOrRigth (i) == -1 % is to the left  angles will be positive
                temp4 = temp3 - 90;
                temp4(temp4<0) = max(temp4);
                ZLateralReference(i,[1 2]) = unique(temp(temp4 == min(temp4),:));
                %           scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
                %           scatter(LateralReference(i,1),LateralReference(i,2),'b')
                ZDistance2Lateral(i) = sqrt((ZLateralReference(i,1) -REF_point_pos(i,1))^2 + (ZLateralReference(i,2) -REF_point_pos(i,2))^2);
            end




            % For the rostrocaudal

            for ii = 1:size(temp,1)
                IthVector = [temp(ii,1)-REF_point_pos(i,1) temp(ii,2)-REF_point_pos(i,2)];
                temp8(ii,1) = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
                %                 arrow([ROI_point_pos(i,1),ROI_point_pos(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
            end
            temp8 = abs(temp8);

            bbb = temp(temp8==min(temp8),:);
            if size(bbb,1)>1
                ZCaudalReference(i,[1 2]) =  bbb(1,:);
            else
                ZCaudalReference(i,[1 2]) = bbb;
                ZDistance2Caudal(i) = sqrt((ZCaudalReference(i,1) -REF_point_pos(i,1))^2 + (ZCaudalReference(i,2) -REF_point_pos(i,2))^2);
            end

            clear bbb

            bbb = temp(temp8==max(temp8),:);
            if size(bbb,1)>1

                ZRostralReference (i,[1 2]) = bbb(1,:);
            else
                ZRostralReference (i,[1 2]) = bbb;
                ZDistance2Rostral(i) = sqrt((ZRostralReference(i,1) -REF_point_pos(i,1))^2 + (ZRostralReference(i,2) -REF_point_pos(i,2))^2);
            end

            clear bbb
        end





%     Testing zmedial a zlateral zreference.
%     close all
    figure(3)
    imagesc(CALCIUM.zstackImage)
    colormap('jet')
    hold on

    for i = 1:length(CALCIUM.ZLatreference.manual_poly) % carefull now
%         is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.ZLatreference.manual_poly{i}; % will contain the x and y values of all the
%         point that contains the midline ROI. This temp variable will be
%         use later. dont erase it.
    end
    scatter(temp(:,1),temp(:,2),'w')

    for i = 1:length(CALCIUM.ZMIDreference.manual_poly) % carefull now
%         is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.ZMIDreference.manual_poly{i}; % will contain the x and y values of all the
%         point that contains the midline ROI. This temp variable will be
%         use later. dont erase it.
    end
    scatter(temp(:,1),temp(:,2),'w')

    i = 1; % # Number of the cell
    scatter(ZLateralReference(i,1),ZLateralReference(i,2),'r','filled')
    scatter(ZMedialReference(i,1),ZMedialReference(i,2),'r','filled')
    scatter(ZRostralReference(i,1),ZRostralReference(i,2),'g','filled')
    scatter(ZCaudalReference(i,1),ZCaudalReference(i,2),'g','filled')

    scatter(REF_point_pos(i,1),REF_point_pos(i,2),'w','filled')
    
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('ZstackReferenceCell_',char(CALCIUM.list(nfish,1)),'.fig')))
uiwait



    % Lateromedial Normalization.
    %Distance from the medial reference to the lateral reference for each
    %cell
    ZMedial2LateralDistance = sqrt((ZLateralReference(:,1) -ZMedialReference(:,1)).^2 + (ZLateralReference(:,2) -ZMedialReference(:,2)).^2);
    ZRostral2CauldaDistance = sqrt((ZRostralReference(:,1) -ZCaudalReference(:,1)).^2 + (ZRostralReference(:,2) -ZCaudalReference(:,2)).^2);

    % Normalizing each cell with respect to its own lateral edge of the tisue
    ZOwnLateralDistance = (ZDistance2Midline./ ZMedial2LateralDistance).*ZLeftOrRigth; %-1 left lateral 0 medial 1 right lateral

     % Normalizing each cell with respect to its own Rostral edge of the tisue
    ZOwnRostralDistance = (ZDistance2Caudal./ZRostral2CauldaDistance); % 1 Rostral 0 Caudal

    % Normalazing each cell with respect to the most lateral edge of the
    % tissue. (Imagine HB like a rectangular box.this is the distance to the lateral of the box
    % which is in the most lateral point)
    ZGeneralLateralDistance = (ZDistance2Midline./ max(ZMedial2LateralDistance)).*ZLeftOrRigth;

      % Normalazing each cell with respect to the most Rostral edge of the
    % tissue 
    ZGeneralRostralDistance = ZDistance2Caudal./ max(ZRostral2CauldaDistance);

    
    
    ZCoordinates(:,1) = ZGeneralLateralDistance; 
    ZCoordinates(:,2) =  ZGeneralRostralDistance; 
    ZCoordinates(:,3) = CALCIUM.Zposition.CurrentRef.*ones(1,size(REF_point_pos,1)); 
    CALCIUM.ZCoordinates = ZCoordinates; 
    

    %% Comparison and extrapolation
% 
%     cf = ZMedial2LateralDistance / Medial2LateralDistance(jjjj);
   cf = ZDistance2Midline/ Distance2Midline(jjjj);
%         cfX = ZDistance2Midline/Distance2Midline(jjjj);  % Base on the resolution (Num px in the zstack/ Num px in the Snap) in the X direction
%        cfY = 1012/3788; % Base on the resolution (Num px in the zstack/ Num px in the Snap) in the Y direction


%     test = sqrt((ROI_point_pos(1,1)-ROI_point_pos(2,1))^2 + (ROI_point_pos(1,2)-ROI_point_pos(2,2))^2);

 for i = 1:size(ROI_point_pos,1)

    Distance2RefX(i) = (ROI_point_pos(i,1)-ROI_point_pos(jjjj,1));
    Distance2RefY(i) = (ROI_point_pos(i,2)-ROI_point_pos(jjjj,2)); 


   Distance2RefX = Distance2RefX*cf;
     Distance2RefY = Distance2RefY*cf;


   ROI_NewLoc (i,1)= REF_point_pos(1,1) + Distance2RefX(i);
   ROI_NewLoc (i,2) = REF_point_pos(1,2) + Distance2RefY(i); 

   
 end 


 %Test 

 %     Testing zmedial a zlateral zreference.
    figure(3)

   for i = 1:size(ROI_point_pos,1)

    scatter(ROI_NewLoc (i,1),ROI_NewLoc (i,2),'w','filled')

   end 




   % Now finding the references and normalization


    %%%% Midreference

         clear temp2 temp3 temp temp8 temporal

    for i = 1:length(CALCIUM.ZMIDreference.manual_poly) % carefull now
            % is made just for one ROI (MID REFERENCE), the midline
            temp = CALCIUM.ZMIDreference.manual_poly{i}; % will contain the x and y values of all the
            %point that contains the midline ROI. This temp variable will be
            %use later. dont erase it.
    end


     % Calculating the reference point in the midline for each cell
        % cell in a loop
        ZMedialReference = NaN(size(ROI_NewLoc,1),2);
        ZDistance2Midline = NaN(size(ROI_NewLoc,1),1);
        ZLeftOrRigth= NaN(size(ROI_NewLoc,1),1);
%         Ztemp2= NaN(1, size(temp,1));

        for i = 1:size(ROI_NewLoc,1)
            for ii = 1:size(temp,1)
                temp2(ii) = sqrt((temp(ii,1)-ROI_NewLoc(i,1))^2 + (temp(ii,2)-ROI_NewLoc(i,2))^2);
            end
            temporal =  temp(temp2 == min(temp2),:); 
            ZMedialReference(i,[1 2]) =temporal(1,:); % closest point in the midline to the ith neuron
            ZDistance2Midline(i) = unique(min(temp2)); % Distance from closest point to the ith neuron

            % to figure out if the ith neuron is to the left or to the right we
            % will create a vector form the MedialReference to one with a
            % higher y value. this is creating a direction (up and maybe tilted if the brain is tilted)
            % After that we will create a vector for the MedialReference tot
            % the ith neuron. if the angle is created clock wise means that is
            % to the rigth, if is created counterclockwise is to the left.
            temp4 = find(temp(:,2)>ZMedialReference(i,2));
            if isempty(temp4)
                greater(i,[1,2]) = greater(i-1,[1,2]); 
            else
                greater(i,[1,2]) = temp(temp4(1),:);
            end 
            
            UpVector = [greater(i,1)-ZMedialReference(i,1) greater(i,2)-ZMedialReference(i,2)];
            IthVector = [ROI_NewLoc(i,1)-ZMedialReference(i,1) ROI_NewLoc(i,2)-ZMedialReference(i,2)];
            %figure, scatter(ROI_point_pos(i,1),ROI_point_pos(i,2)), hold on, scatter(MedialReference(i,1),MedialReference(i,2))
            %  figure,arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
            %  arrow([MedialReference(i,1),MedialReference(i,2)],[ROI_point_pos(i,1),ROI_point_pos(i,2)],1,30,60)

            temp3 = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));

            if temp3 > 0
                ZLeftOrRigth(i,:) = -1; % to the left of the midline
            elseif temp3 < 0

                ZLeftOrRigth(i,:) = 1;  % to the rigth of the midline
            else
                ZLeftOrRigth = NaN;
            end

           
        end

         clear temp2 temp3 temp temp8 temporal

        % Final Variables --> MedialReference, Distance2Midline, LeftOrRigth


         %%%% Lateralreference

  
       for i = 1:length(CALCIUM.ZLatreference.manual_poly) % carefull now
            % is made just for one ROI (MID REFERENCE), the midline
            temp = CALCIUM.ZLatreference.manual_poly{i}; % will contain the x and y values of all the
            %point that contains the midline ROI. This temp variable will be
            %use later. dont erase it.
        end

        % Calculating the reference point in the midline for each cell
        % cell in a loop
        close all
        ZLateralReference = NaN(size(ROI_NewLoc,1),2);
        ZDistance2Lateral = NaN(size(ROI_NewLoc,1),1);
        ZRostralReference = NaN(size(ROI_NewLoc,1),2);
        ZDistance2Rostral = NaN(size(ROI_NewLoc,1),1);
        ZCaudalReference = NaN(size(ROI_NewLoc,1),2);
        ZDistance2Caudal = NaN(size(ROI_NewLoc,1),1);



        for i = 1:size(ROI_NewLoc,1)

            UpVector = [greater(i,1)-ZMedialReference(i,1) greater(i,2)-ZMedialReference(i,2)];

            % For the lateral

            for ii = 1:size(temp,1)
                temp2= [temp(ii,1)-ZMedialReference(i,1) temp(ii,2)-ZMedialReference(i,2)];
                temp3(ii,1) = atan2d(UpVector(1)*temp2(2)-UpVector(2)*temp2(1),UpVector(1)*temp2(1)+UpVector(2)*temp2(2));
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
            end


            if ZLeftOrRigth (i) == 1 % is to the right angles will be negative
                temp4 = temp3 + 90;
                temp4(temp4<0) = max(temp4);
                ZLateralReference(i,[1 2]) = unique(temp(temp4 == min(temp4),:));
                %         scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
                %         scatter(LateralReference(i,1),LateralReference(i,2))
                ZDistance2Lateral(i) = sqrt((ZLateralReference(i,1) -ROI_NewLoc(i,1))^2 + (ZLateralReference(i,2) -ROI_NewLoc(i,2))^2);
            elseif ZLeftOrRigth (i) == -1 % is to the left  angles will be positive
                temp4 = temp3 - 90;
                temp4(temp4<0) = max(temp4);
                ZLateralReference(i,[1 2]) = unique(temp(temp4 == min(temp4),:));
                %           scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
                %           scatter(LateralReference(i,1),LateralReference(i,2),'b')
                ZDistance2Lateral(i) = sqrt((ZLateralReference(i,1) -ROI_NewLoc(i,1))^2 + (ZLateralReference(i,2) -ROI_NewLoc(i,2))^2);
            end




            % For the rostrocaudal

            for ii = 1:size(temp,1)
                IthVector = [temp(ii,1)-ROI_NewLoc(i,1) temp(ii,2)-ROI_NewLoc(i,2)];
                temp8(ii,1) = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));
                %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
                %                 arrow([ROI_point_pos(i,1),ROI_point_pos(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
            end
            temp8 = abs(temp8);

            bbb = temp(temp8==min(temp8),:);
            if size(bbb,1)>1
                ZCaudalReference(i,[1 2]) =  bbb(1,:);
            else
                ZCaudalReference(i,[1 2]) = bbb;
                ZDistance2Caudal(i) = sqrt((ZCaudalReference(i,1) -ROI_NewLoc(i,1))^2 + (ZCaudalReference(i,2) -ROI_NewLoc(i,2))^2);
            end

            clear bbb

            bbb = temp(temp8==max(temp8),:);
            if size(bbb,1)>1

                ZRostralReference (i,[1 2]) = bbb(1,:);
            else
                ZRostralReference (i,[1 2]) = bbb;
                ZDistance2Rostral(i) = sqrt((ZRostralReference(i,1) -ROI_NewLoc(i,1))^2 + (ZRostralReference(i,2) -ROI_NewLoc(i,2))^2);
            end

            clear bbb
        end


         clear temp2 temp3 temp temp8 temporal


   

    

   
 
% TESTING !!!!!!!!!!!! PLOT ALL OF THEM
figure(1)
 imagesc(CALCIUM.mean); colormap('bone'); hold on

for nroi = 1:size(CALCIUM.roi.manual_poly,1)
    coord = CALCIUM.roi.manual_poly{nroi};
    plot(coord(:,1), coord(:,2), 'LineWidth', 1); hold on
    text(mean(coord(:,1)), mean(coord(:,2)), CALCIUM.roi.labels(nroi),'Color','red','FontSize',14);
end

    for j = 1:CALCIUM.neurons
        figure(j+1)
           imagesc(CALCIUM.zstackImage)
    colormap('jet')
    hold on


    for i = 1:length(CALCIUM.ZLatreference.manual_poly) % carefull now
%         is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.ZLatreference.manual_poly{i}; % will contain the x and y values of all the
%         point that contains the midline ROI. This temp variable will be
%         use later. dont erase it.
    end
    scatter(temp(:,1),temp(:,2),'w')

    for i = 1:length(CALCIUM.ZMIDreference.manual_poly) % carefull now
%         is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.ZMIDreference.manual_poly{i}; % will contain the x and y values of all the
%         point that contains the midline ROI. This temp variable will be
%         use later. dont erase it.
    end
    scatter(temp(:,1),temp(:,2),'w')

    
    scatter(ZLateralReference(j,1),ZLateralReference(j,2),'r','filled')
    scatter(ZMedialReference(j,1),ZMedialReference(j,2),'r','filled')
    scatter(ZRostralReference(j,1),ZRostralReference(j,2),'g','filled')
    scatter(ZCaudalReference(j,1),ZCaudalReference(j,2),'g','filled')

    scatter(ROI_NewLoc(j,1),ROI_NewLoc(j,2),100,'w','filled')
%     i = 2
%     scatter(ROI_point_pos(i,1),ROI_point_pos(i,2),'w','filled')
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('ZCoordinatesRefernce_','CEEL',num2str(j),'-',char(CALCIUM.list(nfish,1)),'.fig')))
uiwait

    end 
    close all





    % Lateromedial Normalization.
    %Distance from the medial reference to the lateral reference for each
    %cell
    ZMedial2LateralDistance = sqrt((ZLateralReference(:,1) -ZMedialReference(:,1)).^2 + (ZLateralReference(:,2) -ZMedialReference(:,2)).^2);
    ZRostral2CauldaDistance = sqrt((ZRostralReference(:,1) -ZCaudalReference(:,1)).^2 + (ZRostralReference(:,2) -ZCaudalReference(:,2)).^2);

    % Normalizing each cell with respect to its own lateral edge of the tisue
    ZOwnLateralDistance = (ZDistance2Midline./ ZMedial2LateralDistance).*ZLeftOrRigth; %-1 left lateral 0 medial 1 right lateral

     % Normalizing each cell with respect to its own Rostral edge of the tisue
    ZOwnRostralDistance = (ZDistance2Caudal./ZRostral2CauldaDistance); % 1 Rostral 0 Caudal

    % Normalazing each cell with respect to the most lateral edge of the
    % tissue. (Imagine HB like a rectangular box.this is the distance to the lateral of the box
    % which is in the most lateral point)
    ZGeneralLateralDistance = (ZDistance2Midline./ max(ZMedial2LateralDistance)).*ZLeftOrRigth;

      % Normalazing each cell with respect to the most Rostral edge of the
    % tissue 
    ZGeneralRostralDistance = ZDistance2Caudal./ max(ZRostral2CauldaDistance);


    clear ZCoordinates 
    CALCIUM.ZDistance2Midline = ZDistance2Midline; 
    CALCIUM.ZMedialReference = ZMedialReference; 
    ZCoordinates(:,1) = ZGeneralLateralDistance; 
    ZCoordinates(:,2) =  ZGeneralRostralDistance; 
    ZCoordinates(:,3) = CALCIUM.Zposition.CurrentRef.*ones(1,size(REF_point_pos,1)); 
    CALCIUM.ZCoordinates = ZCoordinates; 
    CALCIUM.ZLeftOrRigth = ZLeftOrRigth; 

CALCIUMimg('save',CALCIUM,[],list,nfish);



%}
%%%%%%%%