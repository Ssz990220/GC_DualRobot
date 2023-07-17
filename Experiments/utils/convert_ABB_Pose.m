function Ts = convert_ABB_Pose(last_link_pos)
%CONVERT_REAL_ROBOT_POS Summary of this function goes here
%   This function convert last joint posture read on the pendant to a
%   standard SE3 matrix
Ts = zeros([4,4,size(last_link_pos,2)]);
for i = 1:size(last_link_pos,2)
    T = last_link_pos(1:3,i);
    q = last_link_pos(4:7,i)';
%     R =[2*q(1)^2-1+2*q(2)^2, 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
%             2*(q(2)*q(3)+q(1)*q(4)), 2*q(1)^2-1+2*q(3)^2, 2*(q(3)*q(4)-q(1)*q(2));
%             2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 2*q(1)^2-1+2*q(4)^2];
    R = quat2rotm(q);
    Ts(:,:,i) = [R,T;zeros(1,3),1];
end
end