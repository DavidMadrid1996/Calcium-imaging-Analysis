S_1_Loading_Animal

    NameDestin = strcat(char(CALCIUM.list(nfish,1)),'.xlsx');

    for i = 1:length(fieldnames(CALCIUM.AMPLITUDE_X)) % Loop across all the cells (number of cells)

        temp1 = fieldnames(CALCIUM.AMPLITUDE_X);
        myfield1 = char(temp1(i)); % extracting the fieldnames of each cell of the loop
        CellName {i} = myfield1;

    end

    % Writning the Calcium traces

    writematrix(CALCIUMroiTS.diff_perc03.data,fullfile(CALCIUM.path.data,'dataCALCIUM',...
        char(CALCIUM.list(nfish,1)), NameDestin),'Sheet','Calcium Traces','Range','A2')

    writecell(CellName,fullfile(CALCIUM.path.data,'dataCALCIUM',...
        char(CALCIUM.list(nfish,1)), NameDestin),'Sheet','Calcium Traces','Range','A1')

    % Writing the swimming trace

    writematrix(wave,fullfile(CALCIUM.path.data,'dataCALCIUM',...
        char(CALCIUM.list(nfish,1)), NameDestin),'Sheet','Swimming Trace','Range','A2')

    writecell({'Swimming'},fullfile(CALCIUM.path.data,'dataCALCIUM',...
        char(CALCIUM.list(nfish,1)), NameDestin),'Sheet','Swimming Trace','Range','A1')


    for i = 1:length(fieldnames(CALCIUM.EPISODE_X)) % Loop across swm episodes

        temp1 = fieldnames(CALCIUM.EPISODE_X);
        myfield1 = char(temp1(i)); % extracting the fieldnames of each cell of the loop
        Episode{i} = myfield1;
        Swm = table;
        Swm(:,1) = array2table(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)));  % Amplitude in degrees
        Swm(:,2) = array2table([CALCIUM.FREQ_SWIM_Y.(myfield1)' ; 0]); % Freq in degrees
        Swm(1,3) = array2table(mean(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)))); % Mean Amplitude in degrees
        Swm(1,4) = array2table(std(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)))); % std Amplitude in degrees
        Swm(1,5) = array2table(max(abs(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1)))); % Max Amplitude in degrees
        Swm(1,6) = array2table(mean(CALCIUM.FREQ_SWIM_Y.(myfield1))); % Mean Freq in degrees
        Swm(1,7) = array2table(std(CALCIUM.FREQ_SWIM_Y.(myfield1))); % std Freq in degrees
        Swm(1,8) = array2table(max(CALCIUM.FREQ_SWIM_Y.(myfield1))); % Max Freq in degrees
        Swm(1,9) = array2table(CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield1)(end)-CALCIUM.FINAL_POINTS_SWIMMING_X.(myfield1)(1)); % Duration Episode Seconds
        Swm(1,10) = array2table(length(CALCIUM.EPISODE_Y_BASELINE_AMP.(myfield1))); % Number Half Cycles



        Swm.Properties.VariableNames(1:10) = {'Amp [Degrees] ', 'Freq [Hz]',...
            'Mean Amp [Degrees]','std Amp [Degrees]','Max Amp [Degrees]','Mean Freq [Hz]','std Freq [Hz]','Max Freq [Hz]',...
            'Duration Epiosde [s]', 'Number Half Cycles' };

        for ii = 1:length(fieldnames(CALCIUM.AMPLITUDE_X)) % Loop across all the cells (number of cells)

            temp1 = fieldnames(CALCIUM.AMPLITUDE_X);
            myfield = char(temp1(ii)); % extracting the fieldnames of each cell of the loop

            if CALCIUM.RESPONSE.(myfield)(i)~=0 % If is not equal to zero means that the neuron is responding to the swimming episode
                CalOfEpi(ii) = CALCIUM.BASELINE_TO_PEAK.(myfield)(CALCIUM.RESPONSE.(myfield)(i));
%                 DelOfEpi(ii) = CALCIUM.DELAY.(myfield)(i); 
                Binary(ii) = 1;
            else
                CalOfEpi(ii) = 0;
%                 DelOfEpi(ii) = CALCIUM.DELAY.(myfield)(i);  
                Binary(ii) = 0;
            end


        end

        Cal = table;
        Cal(1,:) = array2table(Binary)
        Cal(2,:) = array2table(CalOfEpi);
%         Cal(3,:) = array2table(DelOfEpi); 
        Cal([3 4 5],:) = array2table( CALCIUM.ZAbsoluteCoordinates')

        Cal.Properties.VariableNames(1:length(CellName)) = CellName ;


        writetable(Swm,fullfile(CALCIUM.path.data,'dataCALCIUM',...
            char(CALCIUM.list(nfish,1)), NameDestin),'Sheet',myfield1,'Range','A1')

        writetable(Cal,fullfile(CALCIUM.path.data,'dataCALCIUM',...
            char(CALCIUM.list(nfish,1)), NameDestin),'Sheet',myfield1,'Range','L1')

        clear CalOfEpi Binary
    end


