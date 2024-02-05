function selectPoint(~, eventData)
    % Get the current point clicked
    pt = eventData.IntersectionPoint;
    x = pt(1);
    y = pt(2);
    z = pt(3);
   
    % Retrieve the existing points from the workspace
    xdata = evalin('base', 'xdata');
    ydata = evalin('base', 'ydata');
    zdata = evalin('base', 'zdata');
   
    % Find the closest point
    distances = sqrt((xdata-x).^2 + (ydata-y).^2 + (zdata-z).^2);
    [~, idx] = min(distances);
    closest_point = [xdata(idx), ydata(idx), zdata(idx)];
   
    % Add the selected point to the matrix
    selectedPoints = evalin('base', 'selectedPoints');
    selectedPoints = [selectedPoints; closest_point];
   
    % Display the selected point
    fprintf('Selected point: (%.2f, %.2f, %.2f)\n', closest_point);
   
    % Update the matrix in the workspace
    assignin('base', 'selectedPoints', selectedPoints);

end