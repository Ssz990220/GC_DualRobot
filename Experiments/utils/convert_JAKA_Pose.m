function Ts = convert_JAKA_Pose(Poses)
%CONVERT_JAKA_POSE 此处显示有关此函数的摘要
%   此处显示详细说明
Ts = zeros([4,4,size(Poses,2)]);

for i = 1:size(Poses,2)
    T = Poses(1:3,i);
    q = (Poses(4:6,i)/180*pi)';
    R = eul2rotm([q(3),q(2),q(1)]);
    Ts(:,:,i) = [R,T;zeros(1,3),1];
end
end

