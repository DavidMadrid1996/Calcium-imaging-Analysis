close all
% Example data (replace this with your actual neuron trace matrix)
% Assuming your data is stored in a variable named 'neuron_traces'
% Replace this with your actual data, ensuring it is in the correct format
neuron_traces = tempfilter; % Example matrix of 100 neurons with 500 time points each

% Compute DTW distances between neuron traces
dtw_distance_matrix = zeros(size(neuron_traces, 1));

for i = 1:size(neuron_traces, 1)
    for j = 1:size(neuron_traces, 1)
        dtw_distance_matrix(i, j) = dtw(neuron_traces(i, :), neuron_traces(j, :));
    end
end

% Perform hierarchical clustering on the DTW distance matrix
linkage_method = 'average'; % You can try different linkage methods
cluster_tree = linkage(dtw_distance_matrix, linkage_method);

% Visualize the dendrogram
figure;
dendrogram(cluster_tree);
title('Hierarchical Clustering Dendrogram (DTW)');

% Determine cluster assignments based on a threshold or criterion
% Here, 'cutoff' represents the threshold for forming clusters
cutoff = 10; % Replace this with an appropriate threshold
cluster_indices = cluster(cluster_tree, 'cutoff', cutoff, 'criterion', 'distance');

% Visualize the clustered data (for demonstration purposes)
% Replace this visualization with your preferred method
figure;
scatter3(neuron_traces(:, 1), neuron_traces(:, 2), neuron_traces(:, 3), 50, cluster_indices, 'filled');
xlabel('Dimension 1');
ylabel('Dimension 2');
zlabel('Dimension 3');
title('Clustered Neuron Traces (DTW)');
