
function [coordinateSelected, minIdx] = showZValueFcn(hObj, event)
%  FIND NEAREST (X,Y,Z) COORDINATE TO MOUSE CLICK
% Inputs:
%  hObj (unused) the axes
%  event: info about mouse click
% OUTPUT
%  coordinateSelected: the (x,y,z) coordinate you selected
%  minIDx: The index of your inputs that match coordinateSelected. 
x = hObj.XData; 
y = hObj.YData; 
z = hObj.ZData; 
pt = event.IntersectionPoint;       % The (x0,y0,z0) coordinate you just selected
coordinates = [x(:),y(:),z(:)];     % matrix of your input coordinates
dist = pdist2(pt,coordinates);      %distance between your selection and all points
[~, minIdx] = min(dist);            % index of minimum distance to points
coordinateSelected = coordinates(minIdx,:); %the selected coordinate
% from here you can do anything you want with the output.  This demo
% just displays it in the command window.  
fprintf('[x,y,z] = [%.5f, %.5f, %.5f]\n', coordinateSelected)

end % <--- optional if this is embedded into a function