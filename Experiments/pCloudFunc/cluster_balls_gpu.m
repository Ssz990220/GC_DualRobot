function labels = cluster_balls_gpu(balls_Location) %# gpu codegen
    % Add kernelfun pragma to trigger kernel creation
%     coder.gpu.kernelfun;
    balls = pointCloud(balls_Location);
    minDistance = 3;
    [labels,~] = pcsegdist(balls,minDistance);
%     point_counts = zeros(1,numClusters);
%     for i = 1:numClusters
%         point_counts(i) = sum((labels == i));
%     end
%     [~,ball_labels] = gpucoder.sort(point_counts,'descend');
%     ball_labels = ball_labels(1:3);
end

