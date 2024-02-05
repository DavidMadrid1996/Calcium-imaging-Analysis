
    %%%%%%%%%% Finding onset of the calcium wave
    for kk = 1:length(pksCAL) % Loop across waves (peaks)

        % FINDING ONSET REFERING TO THE A % OF THE PEAK
        %{
          rising_threshold = 0.2; % Amplitud must be 20% of the original peak
          rising_points = pksCAL(kk)*0.2; % value of the amplitude of the rising point of the first wave (kk)
          for the first neuron (i)
        %}


        % Because there can be several waves, there will be several rising_points
        % that will be equal to the rising_points of one cell. To avoid this, we
        % have to restrict the time points where we will find this rising_points to
        % points that just precede the peak that we are analysing

        time_segment(1) = locsCAL(kk)*CALCIUM.srate-50*CALCIUM.srate; % start of the segment
        if time_segment(1)<0
            time_segment(1)=CALCIUM.srate;
        end
        % ntimes*sampling rate before the peak
        time_segment(2) = locsCAL(kk)*CALCIUM.srate; % end of the segment is the
        % moment in seconds of the peak
        clear index_segment

        new_vector_time = time_segment(1):CALCIUM.srate:time_segment(2);


        for kkk = 1:length(new_vector_time)
            index_segment (kkk)= find(round(CALCIUMroiTS.diff_perc03.times,4)==...
                round(new_vector_time(kkk),4));  % indexs of the time CALCIUMroiTS.diff_perc03.times
            % that correspond with ntimes*sampling rate before the peak. we will
            % use this to index the CALCIUMroiTS.diff_perc03.data
        end


        % FINDING ONSET REFERING TO THE A % OF THE PEAK
        %{
          temp_rising = abs(rising_points-abs(CALCIUMroiTS.diff_perc03.data(index_segment,i)));
          index_rising = find(temp_rising==min(temp_rising));
          time_rising = new_vector_time(index_rising(end)); % these will be the temporal indexs
        %}

        segment_data = CALCIUMroiTS.diff_perc03.data(i,index_segment);

        for jj = 1:length(CALCIUMroiTS.diff_perc03.data(i,index_segment))-1% rate of change
            slope(jj)=segment_data(jj)...
                -segment_data(jj+1);
        end

        thresholh_segment = 6;
        sumation = 1;
        maxim_threshold = 0;
        index = min(find(slope>thresholh_segment));
        while size(index,2) == 0 % Can be that the original thershold was too big and there is not
            % any nuber above it. It will create a problem in the next part
            % of the code. For this reason
            thresholh_segment = thresholh_segment-sumation;
            index = min(find(slope>thresholh_segment));
            maxim_threshold = 1;


            % This while loop wil continue until there is an aceptable value
            % for the threshold

        end
        if maxim_threshold==1 % To inform to the user
            msgbox('Highest possible value for threshold reached')
            pause(2)
            clear maxim_threshold
        end

        CALCIUM.ONSET_X.(my_field)(kk) = find(round(CALCIUMroiTS.diff_perc03.times,4)==...
            round(new_vector_time(index),4));

        scatter(CALCIUM.ONSET_X.(my_field)(kk),-segment_data(index),30,'filled','r')

        prompt = {'Do you wanna set a new threshold','Which one?','Do you prefer using the value of the peak?'};
        dlgtitle = 'Input';
        dims = [1 35];

        definput ={'N',num2str(thresholh_segment),'N'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        answer1 =  answer(1);
        answer3 = answer(3);

        if isequal(answer3,{'Y'}) %if your wave comes after a second one, maybe
            % you dont have a proper onset, in this case, take the peak
            CALCIUM.ONSET_X.(my_field)(kk) = locsCAL(kk);
            scatter(CALCIUM.ONSET_X.(my_field)(kk),pksCAL(kk),30,'filled','b')

        elseif isequal(answer1,{'Y'})
            while isequal(answer1,{'Y'})
                thresholh_segment = str2double(cell2mat(answer(2)));
                index = min(find(slope>thresholh_segment))
                maxim_threshold = 0;
                while size(index,2) == 0 % Can be that the original thershold was too big and there is not
                    % any nuber above it. It will create a problem in the next part
                    % of the code. For this reason
                    thresholh_segment = thresholh_segment-sumation;
                    index = min(find(slope>thresholh_segment));
                    maxim_threshold = 1;


                    % This while loop wil continue until there is an aceptable value
                    % for the threshold

                end
                if maxim_threshold==1 % To inform to the user
                    msgbox('Highest possible value for threshold reached')
                    pause(2)
                    clear maxim_threshold

                end

                CALCIUM.ONSET_X.(my_field)(kk) = (find(round(CALCIUMroiTS.diff_perc03.times,4)==...
                    round(new_vector_time(index),4)));

                scatter(CALCIUM.ONSET_X.(my_field)(kk),-segment_data(index),30,'filled','b')

                prompt = {'Do you wanna set a new threshold','Which one?'};
                dlgtitle = 'Input';
                dims = [1 35];

                definput ={'Y',num2str(thresholh_segment)};
                answer = inputdlg(prompt,dlgtitle,dims,definput);
                answer1 =  answer(1);

            end
        end

        CALCIUM.ONSET_X.(my_field)(kk) = CALCIUM.ONSET_X.(my_field)(kk)*CALCIUM.srate;


        saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',...
            char(CALCIUM.list(nfish,1)),strcat('Peaks_and_onset',my_field, char(CALCIUM.list(nfish,1)))))

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Finally we clear the pksCAL and locsCAL
    clear locsCAL pksCAL sumation
    w = waitforbuttonpress; % Once the scatter is plotted just click any key to go to next celL
    close all

