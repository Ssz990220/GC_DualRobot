function [ball_labels,labels] = cluster_balls(balls,k)
    minDistance = 2;
    [labels,numClusters] = pcsegdist(balls,minDistance);
    point_counts = zeros(1,numClusters);
    for i = 1:numClusters
        point_counts(i) = sum((labels == i));
    end
    n_balls = sum(point_counts>80000);         % Valid point cloud should posses more than 1e6 points
    if n_balls<3
        fprintf("Only find %d points on dataset %d \n", n_balls,k)
    end
    [~,ball_labels] = sort(point_counts,'descend');
    ball_labels = ball_labels(1:3);
end

