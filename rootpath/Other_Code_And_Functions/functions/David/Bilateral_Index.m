function [BilateralIndex] = Bilateral_Index(Amplitude, Amplitude_Sign)
% Amplitude = Variable; 
% Amplitude_Sign = Variable2;  

  limit = 5; % how many degrees ayou are willing to accept for not crossing the midline
  


  if length(Amplitude) == 1 % Just one movement, 
      Bilateral = 0; % so it is unilateral by definition
  else
     
      for i = 2:length(Amplitude)


          if  Amplitude_Sign(i-1)~=Amplitude_Sign(i) % If the sign of the previous is not the same than the sign of the current,
              Bilateral(i-1) = 1;  % then it means that there was a cross of the midline, so it is bilateral
          elseif Amplitude_Sign(i-1)==Amplitude_Sign(i) && abs(Amplitude(i))<limit  % If the sign of the previous is the same than
              % the sign of the current (there is no cross of the midline), but
              % is not further away than X limit degrees (giving this space beacuse maybe there are errors in the labelling)
              Bilateral(i-1) = 1; % then it is also bilateral
          else
              Bilateral(i-1) = 0; % otherwise is not bilateral
          end

      end


  end
  BilateralIndex = mean(Bilateral);
end

