
    %Normalize the traces 
%     tail0 = [CALCIUMroiTS.deeplabcut.Mouth.x-CALCIUMroiTS.deeplabcut.Mouth.x CALCIUMroiTS.deeplabcut.Mouth.y-CALCIUMroiTS.deeplabcut.Mouth.y]
    tail1 = [CALCIUMroiTS.deeplabcut.tail1.x-CALCIUMroiTS.deeplabcut.tail1.x CALCIUMroiTS.deeplabcut.tail1.y-CALCIUMroiTS.deeplabcut.tail1.y];
    tail2 = [CALCIUMroiTS.deeplabcut.tail2.x-CALCIUMroiTS.deeplabcut.tail1.x CALCIUMroiTS.deeplabcut.tail2.y-CALCIUMroiTS.deeplabcut.tail1.y];
    tail3 = [CALCIUMroiTS.deeplabcut.tail3.x-CALCIUMroiTS.deeplabcut.tail1.x CALCIUMroiTS.deeplabcut.tail3.y-CALCIUMroiTS.deeplabcut.tail1.y];
    tail4 = [CALCIUMroiTS.deeplabcut.tail4.x-CALCIUMroiTS.deeplabcut.tail1.x CALCIUMroiTS.deeplabcut.tail4.y-CALCIUMroiTS.deeplabcut.tail1.y];
%     tail5 = [CALCIUMroiTS.deeplabcut.tail5.x-CALCIUMroiTS.deeplabcut.tail1.x CALCIUMroiTS.deeplabcut.tail5.y-CALCIUMroiTS.deeplabcut.tail1.y]

for j = 1:length(fieldnames(CALCIUM.EPISODE_X ))
    temp = fieldnames(CALCIUM.EPISODE_X );
    myfield = char(temp (j));
    Selection =  ((CALCIUM.EPISODE_X.(myfield)(1))-30):((CALCIUM.EPISODE_X.(myfield)(end))+50);
    count = 1;

    Len = length(Selection);
    red = [0 0 1];
    pink = [1 1 0];
    colors_p = [linspace(red(1),pink(1),Len)', linspace(red(2),pink(2),Len)', linspace(red(3),pink(3),Len)'];

    theta = 90; % to rotate 90 counterclockwise
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    subplot(3,3,[1 4 ])
    cla

    for i = Selection

        X = [tail1(i,1) tail2(i,1) tail3(i,1) tail4(i,1)];
        Y = [tail1(i,2) tail2(i,2) tail3(i,2) tail4(i,2)];
        %     X = [tail1(i,1) tail2(i,1) tail3(i,1) tail4(i,1) tail5(i,1)];
        %     Y = [tail1(i,2) tail2(i,2) tail3(i,2) tail4(i,2) tail5(i,2)];
        for ii = 1:size(X,2)

            temp =R*[X(ii) Y(ii)]';
            X(ii) = temp(1);
            Y(ii) = temp(2);

        end

        scatter(X,-Y,20,colors_p(count,:,:),'filled'),hold on % 4 if just 4 points to plot, if 5 put 5
        plot(X,-Y,'color',colors_p(count,:,:),'LineWidth',2)
        count = count+1;

    end
    clear count
    ylim([-100 0])
    xlim([-60 60])

    colormap(colors_p)
    colorbar


    % How much do you want to take to have a clear trace
    HowMuch = 10;
    subplot(3,3,[7 8 9])
    cla
    plot(1:length(wave((Selection(1)-HowMuch):(Selection(end)+HowMuch))),wave((Selection(1)-HowMuch):(Selection(end)+HowMuch)),'k'),hold on
    plot((HowMuch+1):(HowMuch+length(Selection)),wave(Selection),'b','LineWidth',2)


    % Final subplot


    subplot(3,3,[2 3 5 6])
    cla
    count  = 1;
    for i = Selection

        X = [tail1(i,1) tail2(i,1) tail3(i,1) tail4(i,1)];
        Y = [tail1(i,2) tail2(i,2) tail3(i,2) tail4(i,2)];
        %     X = [tail1(i,1) tail2(i,1) tail3(i,1) tail4(i,1) tail5(i,1)];
        %     Y = [tail1(i,2) tail2(i,2) tail3(i,2) tail4(i,2) tail5(i,2)];
        for ii = 1:size(X,2)

            temp =R*[X(ii) Y(ii)]';
            X(ii) = temp(1);
            Y(ii) = temp(2);

        end

        scatter(X+count*10 ,-Y,20,colors_p(count,:,:),'filled'),hold on % 4 if just 4 points to plot, if 5 put 5
        plot(X+count*10,-Y,'color',colors_p(count,:,:),'LineWidth',2)
        count = count+1;

    end
    clear count
    %             ylim([-100 0])
    %     xlim([-60 60])

    colormap(colors_p)
    colorbar

    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat(myfield,'_detailed_',char(CALCIUM.list(nfish,1)))))
     uiwait

end


%% Parameters
for i = 1:length(fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP))
    figure(i)
    temp1 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
    myfield1 = char(temp1(i)); % extracting the fieldnames of each cell of the loop
    CellName {i} = myfield1;

% Amplitudes
  subplot(4,3,[1 2 4 5])

%   
%     temp1 = fieldnames(CALCIUM.EPISODE_Y_BASELINE_AMP);
%     myfield1 = char(temp1(i)); % extracting the fieldnames of each cell of the loop
%     CellName {i} = myfield1;

    boxplot(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))
    xticklabels(CellName) 
    title('Amplitude')
    ylabel('Degrees')
 



% Frequencies

  subplot(4,3,[7 8 10 11])

    boxplot(CALCIUM.GENERAL_FREQ_Y.(myfield1))

        xticklabels(CellName) 



    title('Frequencies')
    ylabel('[Hz]')

 
    %% Bilateral swimming
subplot(4,3,[3 6])

   plot(1:length(CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1)),...
       CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1))

title('Bilaterallity')

    %% Bilateral swimming
subplot(4,3,[9 12])


   plot(1:length(CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1)),...
       abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)).*CALCIUM.EPISODE_Y_BASELINE_AMP_SIGN.(myfield1)')

title('Bilaterallity')

end 





