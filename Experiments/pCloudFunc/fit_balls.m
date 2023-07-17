function [ball_collection,ib_collection,Centers,model_b1,model_b2,model_b3] = fit_balls(balls,ball_labels,labels)
maxDistance = 0.5;
index = (1.0:1.0:balls.Count)';
ball1 = select(balls,index(labels==ball_labels(1)));
ball2 = select(balls,index(labels==ball_labels(2)));
ball3 = select(balls,index(labels==ball_labels(3)));
[model_b1,ib1,~] = pcfitsphere(ball1,maxDistance,Confidence=99.9);
[model_b2,ib2,~] = pcfitsphere(ball2,maxDistance,Confidence=99.9);
[model_b3,ib3,~] = pcfitsphere(ball3,maxDistance,Confidence=99.9);
Centers = zeros(3,3);
Centers(1,:) = model_b1.Center;
Centers(2,:) = model_b2.Center;
Centers(3,:) = model_b3.Center;
ball_collection = {ball1,ball2,ball3}; ib_collection = {ib1,ib2,ib3};
end

