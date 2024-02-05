

%% Manual Finding reference with respect to zstack

close all
clear temp2 temp3 temp temp8 temporal
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(CALCIUM.mean); colormap('bone'); hold on

Bright = 0.2;

figure('units','normalized','outerposition',[0 0 1 1])

%      imagesc(CALCIUM.zstackImage), colormap('bone'); % Use if you have the SNAP
%     imagesc(adapthisteq(CALCIUM.zstackImage+ Bright)), colormap('bone'); %
imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'), hold on

for nroi = 1:size(CALCIUM.roi.manual_poly,1)
    figure(1)
    coord = CALCIUM.roi.manual_poly{nroi};
    plot(coord(:,1), coord(:,2), 'LineWidth', 1); hold on
    text(mean(coord(:,1)), mean(coord(:,2)), CALCIUM.roi.labels(nroi),'Color','red','FontSize',14);
    plot(CALCIUMroiTS.roi.manual_poly{nroi}(:,1),CALCIUMroiTS.roi.manual_poly{nroi}(:,2))

    figure(2)
    clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;


    % Each neuron will be treated like a individual point.
    REF_point_pos(nroi,:) = [num_neurons_x num_neurons_y];  % Matrix containing the average x and y position for each neuron
    figure(2)
    scatter( REF_point_pos(nroi,1),REF_point_pos(nroi,2),100,'filled','red')
    text(REF_point_pos(nroi,1), REF_point_pos(nroi,2), CALCIUM.roi.labels(nroi),'Color','red','FontSize',14);
    figure(1)
    cla
    imagesc(CALCIUM.mean); colormap('bone'); hold on
end
close all


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

    unice =  temp(temp2 == min(temp2),:);
    ZMedialReference(i,[1 2]) = unice(1,:); % closest point in the midline to the ith neuron
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
for ii = 1:size(CALCIUM.roi.manual_poly,1)
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


    scatter(ZLateralReference(ii,1),ZLateralReference(ii,2),'r','filled')
    scatter(ZMedialReference(ii,1),ZMedialReference(ii,2),'r','filled')
    scatter(ZRostralReference(ii,1),ZRostralReference(ii,2),'g','filled')
    scatter(ZCaudalReference(ii,1),ZCaudalReference(ii,2),'g','filled')

    scatter(REF_point_pos(ii,1),REF_point_pos(ii,2),'w','filled')

     saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('ZstackReferenceCell_',num2str(ii),char(CALCIUM.list(nfish,1)),'.fig')))
    uiwait
end



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
CALCIUM.ZLeftOrRigth = ZLeftOrRigth;
CALCIUM.ZDistance2Midline = ZDistance2Midline; 
CALCIUM.ZMedialReference = ZMedialReference;
cCALCIUM.ZDistance2Caudal = ZDistance2Caudal; 


CALCIUMimg('save',CALCIUM,[],list,nfish);