function [ROImatrix] = roi_TSave2(Y,roimask)
% [outwave] = roi_TSaverage(Y,roimask). Extracts timeserie of the input
% movie (excluding last frame, in case that the background was included)

% INPUT: 'Y' 3D data matrix (movie: x*y*frames); 'roimask', 2D logic;

% OUTPUT: 'wave' of the ROI timeserie, that is the average value of all the pixels
% in the ROI (roimask)

Nframes= size(Y,3)-1; % substract last frame
ROImatrix=NaN(size(Y,1),size(Y,2),Nframes); % vector of length the nº of frames in data (Y) 

for frame=1:Nframes

        ROIvalues = Y(:,:,frame).*roimask;
        ROIvalues = fillmissing(ROIvalues,'constant', 0); 
        ROImatrix(:,:,frame)=ROIvalues;
           
end

end

%%  Created: 06/12/21 (from old code)
% old function, last debugged on: 21/10/20

% Updated: 07/02/21
