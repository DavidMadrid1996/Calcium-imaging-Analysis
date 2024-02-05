%% This will plot the variable Data_XXX_Blank for the
%  type of neuron that we want

%% Introduce variable to analize
 temp = - [Data_LongDelay_Blank ; Data_ShortDelay_Blank]; 
% temp = - Data_LongDelay_Blank(:,[2400:2600]); 
%imagesc(temp)

%% Normalization

close all
clear tempfilter
prompt = 'Do you want to normalize it?';
str = input(prompt,'s');


% 24
%plot(-Group_Data(19,:))
%plot(Data_NonActive(:,1))
suma = 1; 
for i =1:size(temp)-2
    
    if str=='Y'

        if max(temp(i,:))<0

            tempfilter(suma,:) = medfilt1(temp(i,:)./min(temp(i,:)),8);
        else
            tempfilter(suma,:) = medfilt1(temp(i,:)./max(temp(i,:)),8);
        end
    elseif str=='N'
        tempfilter(suma,:) = medfilt1(temp(i,:));
    end
   plot(tempfilter(suma,:),'B'),hold on
    %         uiwait
    %     tempfilter = temp
    % plot(tempfilter(i,:)),hold on
%     pause(2)
    % disp(i)
    suma = 1 + suma; 
end

imagesc(tempfilter)


%% Removing zeros rows
suma = 1; 
for i = 1:size(tempfilter,1)
    if sum(tempfilter(i,:))>0
        tempfilter_No_Zeros (suma,:)= tempfilter(i,:);  
        suma = 1 + suma; 
    else
    end 
end 
imagesc( tempfilter_No_Zeros)
