clc;
clear;
close all;

%% 1. Import Data
roundNumber = 10;
filePath = sprintf('./pfResults/pfResults_round%d.csv', roundNumber);
particles = readmatrix(filePath);

image = imread("./map/map_II.pgm");
image = flipud(image(1:250, 1:300));
resolution = 0.05;
origin = [-4.000000, -5.000000];

srcGroundTruth_world = [
    [1.25, 0.6]; [1.25, -0.6]; [2.50, 0.6]; [2.50, -0.6]; 
    [3.74, 0.6]; [3.74, -0.6]; [4.99, 0.6]; [4.99, -0.6]; 
    [6.23, 0.6]; [6.23, -0.6]
    ];

srcGroundTruth = (srcGroundTruth_world - origin) / resolution;

%% 2. Define Grid Search Parameters
epsilon_values = 1:1:15; % Epsilon range
min_samples_values = 5:1:20; % MinPts range
best_epsilon = NaN;
best_min_samples = NaN;
best_metric = Inf;

%% 3. Grid Search for Best Hyper-Parameters
for epsilon = epsilon_values
    for min_samples = min_samples_values
        metric = compute_clustering_metrics(particles, srcGroundTruth_world, origin, resolution, epsilon, min_samples);
        if metric < best_metric
            best_metric = metric;
            best_epsilon = epsilon;
            best_min_samples = min_samples;
        end
    end
end

fprintf('Best epsilon: %f, Best min_samples: %d, Best metric: %f\n', best_epsilon, best_min_samples, best_metric);

%% 4. Plot the Best Clustering Results
labels = dbscan_lite(particles, best_epsilon, best_min_samples);
uniqueLabels = unique(labels);

figure;
imshow(image, 'XData', [0, size(image,2)], 'YData', [0, size(image,1)],'InitialMagnification', 300);
set(gca, 'YDir', 'normal');
hold on;
simple_gscatter(particles(:,1), particles(:,2), labels);
scatter(srcGroundTruth(:,1), srcGroundTruth(:,2), 100, 'rp', 'filled','DisplayName', 'Source ground truth');

% Plot estimated source positions as triangles
validClusters = uniqueLabels(uniqueLabels > 0);

legend_added = false;
for i = 1:length(validClusters)
    cluster_points = particles(labels == validClusters(i), :);
    cluster_center = mean(cluster_points, 1);
    estimated_source = cluster_center;
    x = estimated_source(1);
    y = estimated_source(2);
    
    % Add the 'DisplayName' property only once
    if ~legend_added
        patch([x-2 x+2 x], [y-2 y-2 y+2], 'g', 'FaceAlpha', 0.6, 'DisplayName', 'Estimated Source');
        legend_added = true; % Update the flag
    else
        patch([x-2 x+2 x], [y-2 y-2 y+2], 'g', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
    end
end

title(sprintf('DBSCAN Clustering Results (ε=%f, MinPts=%d)', best_epsilon, best_min_samples));
legend('show');
axis equal;
set(gca, 'YDir', 'normal');
axis on;
axis image;
hold off;

% attention: only support matlab 2020a or later
% exportgraphics(gcf, sprintf('./figures/figures_round%d.png', roundNumber), 'Resolution', 600);

%% 5. Estimate Source Positions and Calculate Distances with Best Parameters
clusterCenters = [0 0];
distances = [0];
clusterCounts = [];

validClusters = uniqueLabels(uniqueLabels ~= -1);

for i = 1:length(validClusters)
    cluster_points = particles(labels == validClusters(i), :);
    if i ~= 1
        cluster_center = mean(cluster_points, 1);
        clusterCenters = [clusterCenters; cluster_center];
        estimated_source = cluster_center * resolution + origin;
        dists = sqrt(sum((srcGroundTruth_world - estimated_source).^2, 2));
        distances = [distances; min(dists)];
    end
    clusterCounts = [clusterCounts; numel(cluster_points) / 2];
end

resultsTable = table(validClusters, clusterCenters, clusterCounts, distances, 'VariableNames', {'Cluster', 'Cluster Center (x, y)', 'Number of Points', 'Distance to Nearest Ground Truth (m)'});
disp(resultsTable);

%% dbscan lite
function labels = dbscan_lite(X, epsilon, minPts)
    % Initialize labels for all points to -1 (unvisited)
    labels = -ones(size(X, 1), 1);
    C = 0;  % Cluster counter
    
    % Iterate through each point
    for i = 1:size(X, 1)
        if labels(i) ~= -1
            continue;
        end
        
        % Find neighbors of the point
        neighbors = regionQuery(X, i, epsilon);
        
        % Check if the point is noise
        if numel(neighbors) < minPts
            labels(i) = 0;  % Mark as noise
            continue;
        end
        
        % Create a new cluster
        C = C + 1;
        labels(i) = C;
        
        % Expand the cluster
        k = 1;
        while true
            if k > numel(neighbors)
                break;
            end
            
            j = neighbors(k);
            
            if labels(j) == 0
                labels(j) = C;  % Change noise to border point
            end
            
            if labels(j) ~= -1
                k = k + 1;
                continue;
            end
            
            labels(j) = C;
            new_neighbors = regionQuery(X, j, epsilon);
            
            if numel(new_neighbors) >= minPts
                neighbors = [neighbors; new_neighbors];  % Add new neighbors to the cluster
            end
            
            k = k + 1;
        end
    end
end

function neighbors = regionQuery(X, idx, epsilon)
    % Calculate distances from the point to all other points
    distances = sqrt(sum((X - X(idx, :)).^2, 2));
    % Find neighbors within the epsilon radius
    neighbors = find(distances <= epsilon);
end


%% Function to Compute Clustering Metrics
function metric = compute_clustering_metrics(particles, srcGroundTruth_world, origin, resolution, epsilon, min_samples)
    labels = dbscan_lite(particles, epsilon, min_samples);
    uniqueLabels = unique(labels);
    validClusters = uniqueLabels(uniqueLabels > 0);
    if isempty(validClusters)
        metric = Inf; % Return a high value if no valid clusters found
        return;
    end
    
    distances = [];
    for i = 1:length(validClusters)
        cluster_points = particles(labels == validClusters(i), :);
        cluster_center = mean(cluster_points, 1);
        estimated_source = cluster_center * resolution + origin;
        dists = sqrt(sum((srcGroundTruth_world - estimated_source).^2, 2));
        distances = [distances; min(dists)];
    end
    
    % average without weight
    metric = mean(distances);
end

%% simple gscatter
% I know 'gscatter' is embedded in matlab
% but I do not want to install extra add-in
function simple_gscatter(x, y, group, color_map, markers)
    if nargin < 4
        color_map = lines(max(group) + 1);
    end
    if nargin < 5
        markers = 'o';
    end
    
    hold on;
    unique_groups = unique(group);
    for i = 1:length(unique_groups)
        g = unique_groups(i);
        idx = group == g;
        if g == 0
            scatter(x(idx), y(idx), 16, color_map(g+1, :), '*', 'DisplayName', 'Noise');
        else
            scatter(x(idx), y(idx), 20, color_map(g+1, :), markers, 'filled', 'DisplayName', sprintf('Group %d', g));
        end
    end
end

