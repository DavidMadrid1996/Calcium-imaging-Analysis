function [] = Hist(Data,BinsMin,BinsMax, NumberMins,color,Axis)
% Generate random data
data = Data;

% Create the histogram
binEdges = linspace(BinsMin, BinsMax, NumberMins); % Define the bin edges
[counts, binEdges] = histcounts(data, binEdges);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Plot the histogram bars
% barh(binCenters, counts, 'FaceColor', 'blue', 'EdgeColor', 'none');
hold on;

% Fit a PDF to the data
pd = fitdist(data, 'Kernel'); % Example: normal distribution
y = linspace(min(binCenters), max(binCenters), 100);
pdfValues = pdf(pd, y);

% Calculate the scaling factor
scalingFactor = max(counts) / max(pdfValues);

% Plot the fitted PDF (smooth curve)
if Axis=='V'
plot(pdfValues * scalingFactor, y, 'color',color, 'LineWidth', 2);
% Add labels and title
xlabel('Count');
ylabel('Data');
elseif Axis=='H'
plot( y,pdfValues * scalingFactor, 'color',color, 'LineWidth', 2);
% Add labels and title
xlabel('Data');
ylabel('Count');
else
    error('Valid inputs are H or V')
end 

title('Histogram with Fitted PDF');

% Add legend
end

