
%% Center of mass of each neuron

xpos_N = NaN(1,size(CALCIUM.roi.manual_poly,1)); % Prelocate the vector for the mean x values
ypos_N = NaN(1,size(CALCIUM.roi.manual_poly,1)); % Prelocate the vector for the mean y values

for i = 1:size(CALCIUM.roi.manual_poly,1) % the size will be equal to the number of neurons
    temp = CALCIUM.roi.manual_poly{i}; % for each x and y coordinate of one neuron
    xpos_N(i) = mean (temp(:,1)); % the first column is the x coordinates of each point
    % of one ROI of one neuron
    ypos_N(i)= mean (temp(:,2)); % the second column is the y coordinates of each point
    % of one ROI of one neuron
    clear temp
end

% Each neuron will be treated like a individual point.
ROI_point_pos = [xpos_N' ypos_N'];  % Matrix containing the average x and y position for each neuron
%   scatter( ROI_point_pos(:,1), -ROI_point_pos(:,2))



if isfield(CALCIUM,'Midreference') % checking the condition of the existence of a mid reference

    for i = 1:length(CALCIUM.Midreference.manual_poly) % carefull now
        % is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.Midreference.manual_poly{i}; % will contain the x and y values of all the
        %point that contains the midline ROI. This temp variable will be
        %use later. dont erase it.
    end
    % Test
    imagesc(CALCIUM.mean)
    %         imagesc(CALCIUM.SNAPImage)
    colormap('jet')
    hold on
    scatter(temp(:,1),temp(:,2),'w')

    % Calculating the reference point in the midline for each cell
    % cell in a loop
    MedialReference = NaN(size(ROI_point_pos,1),2);
    Distance2Midline = NaN(size(ROI_point_pos,1),1);
    LeftOrRigth= NaN(size(ROI_point_pos,1),1);
    temp2= NaN(1, size(temp,1));

    for i = 1:size(ROI_point_pos,1)
        for ii = 1:size(temp,1)
            temp2(ii) = sqrt((temp(ii,1)-ROI_point_pos(i,1))^2 + (temp(ii,2)-ROI_point_pos(i,2))^2);
        end
        temporal =  temp(temp2 == min(temp2),:);
        MedialReference(i,[1 2]) = temporal(1,:); % closest point in the midline to the ith neuron
        Distance2Midline(i) = unique(min(temp2)); % Distance from closest point to the ith neuron

        % to figure out if the ith neuron is to the left or to the right we
        % will create a vector form the MedialReference to one with a
        % higher y value. this is creating a direction (up and maybe tilted if the brain is tilted)
        % After that we will create a vector for the MedialReference tot
        % the ith neuron. if the angle is created clock wise means that is
        % to the rigth, if is created counterclockwise is to the left.
        temp4 = find(temp(:,2)>MedialReference(i,2));
        if isempty(temp4)
            greater(i,[1,2]) = greater(i-1,[1,2]);
        else
            greater(i,[1,2]) = temp(temp4(1),:);
        end

        UpVector = [greater(i,1)-MedialReference(i,1) greater(i,2)-MedialReference(i,2)]; % midline vector
        IthVector = [ROI_point_pos(i,1)-MedialReference(i,1) ROI_point_pos(i,2)-MedialReference(i,2)];
        %figure, scatter(ROI_point_pos(i,1),ROI_point_pos(i,2)), hold on, scatter(MedialReference(i,1),MedialReference(i,2))
        %  figure,arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
        %  arrow([MedialReference(i,1),MedialReference(i,2)],[ROI_point_pos(i,1),ROI_point_pos(i,2)],1,30,60)

        temp3 = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));

        if temp3 > 0
            LeftOrRigth(i) = -1; % to the left of the midline
        elseif temp3 < 0

            LeftOrRigth(i) = 1;  % to the rigth of the midline
        else
            LeftOrRigth = NaN;
        end


    end

    clear temp2 temp3 temp temporal

    % Final Variables --> MedialReference, Distance2Midline, LeftOrRigth
end



if isfield(CALCIUM,'Latreference') % checking the condition of the existence of a mid reference

    for i = 1:length(CALCIUM.Latreference.manual_poly) % carefull now
        % is made just for one ROI (MID REFERENCE), the midline
        temp = CALCIUM.Latreference.manual_poly{i}; % will contain the x and y values of all the
        %point that contains the midline ROI. This temp variable will be
        %use later. dont erase it.
    end

    % Calculating the reference point in the midline for each cell
    % cell in a loop
    close all
    LateralReference = NaN(size(ROI_point_pos,1),2);
    Distance2Lateral = NaN(size(ROI_point_pos,1),1);
    RostralReference = NaN(size(ROI_point_pos,1),2);
    Distance2Rostral = NaN(size(ROI_point_pos,1),1);
    CaudalReference = NaN(size(ROI_point_pos,1),2);
    Distance2Caudal = NaN(size(ROI_point_pos,1),1);

    for i = 1:size(ROI_point_pos,1)

        UpVector = [greater(i,1)-MedialReference(i,1) greater(i,2)-MedialReference(i,2)];

        % For the lateral

        for ii = 1:size(temp,1)
            temp2= [temp(ii,1)-MedialReference(i,1) temp(ii,2)-MedialReference(i,2)];
            temp3(ii,1) = atan2d(UpVector(1)*temp2(2)-UpVector(2)*temp2(1),UpVector(1)*temp2(1)+UpVector(2)*temp2(2));
            %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
            %                 arrow([MedialReference(i,1),MedialReference(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
        end


        if LeftOrRigth (i) == 1 % is to the right angles will be negative
            temp4 = temp3 + 90;
            temp4(temp4<0) = max(temp4);
            LateralReference(i,[1 2]) = temp(temp4 == min(temp4),:);
            %         scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
            %         scatter(LateralReference(i,1),LateralReference(i,2))
            Distance2Lateral(i) = sqrt((LateralReference(i,1) -ROI_point_pos(i,1))^2 + (LateralReference(i,2) -ROI_point_pos(i,2))^2);
        elseif LeftOrRigth (i) == -1 % is to the left  angles will be positive
            temp4 = temp3 - 90;
            temp4(temp4<0) = max(temp4);

            LateralReference(i,[1 2]) = unique(temp(temp4 == min(temp4),:));
            %           scatter(ROI_point_pos(i,1),ROI_point_pos(i,2))
            %           scatter(LateralReference(i,1),LateralReference(i,2),'b')
            Distance2Lateral(i) = sqrt((LateralReference(i,1) -ROI_point_pos(i,1))^2 + (LateralReference(i,2) -ROI_point_pos(i,2))^2);
        end




        % For the rostrocaudal

        for ii = 1:size(temp,1)
            IthVector = [temp(ii,1)-ROI_point_pos(i,1) temp(ii,2)-ROI_point_pos(i,2)];
            temp8(ii,1) = atan2d(UpVector(1)*IthVector(2)-UpVector(2)*IthVector(1),UpVector(1)*IthVector(1)+UpVector(2)*IthVector(2));
            %                 arrow([MedialReference(i,1),MedialReference(i,2)],[greater(i,1),greater(i,2)],1,30,60),hold on
            %                 arrow([ROI_point_pos(i,1),ROI_point_pos(i,2)],[temp(ii,1),temp(ii,2)],1,30,60)
        end
        temp8 = abs(temp8);

        bbb = temp(temp8==min(temp8),:);
        if size(bbb,1)>1
            CaudalReference(i,[1 2]) =  bbb(1,:);
        else
            CaudalReference(i,[1 2]) = bbb;
            Distance2Caudal(i) = sqrt((CaudalReference(i,1) -ROI_point_pos(i,1))^2 + (CaudalReference(i,2) -ROI_point_pos(i,2))^2);
        end

        clear bbb

        bbb = temp(temp8==max(temp8),:);
        if size(bbb,1)>1

            RostralReference (i,[1 2]) = bbb(1,:);
        else
            RostralReference (i,[1 2]) = bbb;
            Distance2Rostral(i) = sqrt((RostralReference(i,1) -ROI_point_pos(i,1))^2 + (RostralReference(i,2) -ROI_point_pos(i,2))^2);
        end

        clear bbb
    end

end

%     %     TESTING !!!!!!!! medial a lateral reference.
% for j = 1:size(CALCIUM.roi.manual_poly,1) % # Number of the cell
% 
% 
%     figure(j)
% %     close all
%      imagesc(CALCIUM.mean)
% %     imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone'); % Use if you have the SNAP
%     hold on
% 
%     for i = 1:length(CALCIUM.Latreference.manual_poly) % carefull now
% %         is made just for one ROI (MID REFERENCE), the midline
%         temp = CALCIUM.Latreference.manual_poly{i}; % will contain the x and y values of all the
% %         point that contains the midline ROI. This temp variable will be
% %         use later. dont erase it.
%     end
%     scatter(temp(:,1),temp(:,2),'w')
% 
%     for i = 1:length(CALCIUM.Midreference.manual_poly) % carefull now
% %         is made just for one ROI (MID REFERENCE), the midline
%         temp = CALCIUM.Midreference.manual_poly{i}; % will contain the x and y values of all the
% %         point that contains the midline ROI. This temp variable will be
% %         use later. dont erase it.
%     end
%     scatter(temp(:,1),temp(:,2),'w')
% 
% 
%     scatter(LateralReference(j,1),LateralReference(j,2),'r','filled')
%     scatter(MedialReference(j,1),MedialReference(j,2),'r','filled')
%     scatter(RostralReference(j,1),RostralReference(j,2),'g','filled')
%     scatter(CaudalReference(j,1),CaudalReference(j,2),'g','filled')
% 
%     scatter(ROI_point_pos(j,1),ROI_point_pos(j,2),100,'b','filled')
% %     i = 2
% %     scatter(ROI_point_pos(i,1),ROI_point_pos(i,2),'w','filled')
% 
% end

%     % Lateromedial Normalization.
%     %Distance from the medial reference to the lateral reference for each
%     %cell
%     Medial2LateralDistance = sqrt((LateralReference(:,1) -MedialReference(:,1)).^2 + (LateralReference(:,2) -MedialReference(:,2)).^2);
%     Rostral2CauldaDistance = sqrt((RostralReference(:,1) -CaudalReference(:,1)).^2 + (RostralReference(:,2) -CaudalReference(:,2)).^2);
%
%     % Normalizing each cell with respect to its own lateral edge of the tisue
%     OwnLateralDistance = (Distance2Midline./ Medial2LateralDistance).*LeftOrRigth; %-1 left lateral 0 medial 1 right lateral
%
%      % Normalizing each cell with respect to its own Rostral edge of the tisue
%     OwnRostralDistance = (Distance2Caudal./Rostral2CauldaDistance); % 1 Rostral 0 Caudal
%
%     % Normalazing each cell with respect to the most lateral edge of the
%     % tissue. (Imagine HB like a rectangular box.this is the distance to the lateral of the box
%     % which is in the most lateral point)
%     GeneralLateralDistance = (Distance2Midline./ max(Medial2LateralDistance)).*LeftOrRigth;
%
%       % Normalazing each cell with respect to the most Rostral edge of the
%     % tissue
%     GeneralRostralDistance = Distance2Caudal./ max(Rostral2CauldaDistance);
%
%
%
%     Coordinates(:,1) = GeneralLateralDistance;
%     Coordinates(:,2) =  GeneralRostralDistance;
%     Coordinates(:,3) = CALCIUM.Zposition.CurrentRef.*ones(1,size(ROI_point_pos,1));
%     CALCIUM.Coordinates = Coordinates;
%
%