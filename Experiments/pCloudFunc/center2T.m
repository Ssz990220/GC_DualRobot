function [Ts,orders,dis] = center2T(centers, normals,k)
    % Three balls are enough to define a coordinate. 
    % Say point 3 is the point on y axis, point 1 and 2 sits on x axis
    % The Origin O can be found as the foot on line 12'
    % Z axis is line 12 cross line O3

    % Noted that the order of point 1,2,3 can be flipped so you have to set
    % some rules to determine which point sits on Y.

    if nargin == 2
        k=0;
    end
    N = size(centers,3);
    dis_local = zeros(3,1);
    Ts = zeros(4,4,N);
    for i = 1:N
        dis_local(1) = norm(centers(1,:,i)-centers(2,:,i));
        dis_local(2) = norm(centers(3,:,i)-centers(2,:,i));
        dis_local(3) = norm(centers(3,:,i)-centers(1,:,i));
        [~,order] = sort(dis_local);
        if order(1) == 1
            % centers3 is the Y point
            [Ts(:,:,i),flip] = xyz2T(centers(1,:,i),centers(2,:,i),centers(3,:,i),normals(i,:));
            if flip 
                orders = [2,1,3];
                dis = [dis_local(1);dis_local(3);dis_local(2)];
            else
                orders = [1,2,3];
                dis = [dis_local(1);dis_local(2);dis_local(3)];
            end
        elseif order(1) == 2
            % centers1 is the Y point
            [Ts(:,:,i),flip] = xyz2T(centers(2,:,i),centers(3,:,i),centers(1,:,i),normals(i,:));
            if flip 
                orders = [3,2,1];
                dis = [dis_local(2);dis_local(1);dis_local(3)];
            else
                orders = [2,3,1];
                dis = [dis_local(2);dis_local(3);dis_local(1)];
            end
        else
            % centers2 is the Y point
            [Ts(:,:,i),flip] = xyz2T(centers(1,:,i),centers(3,:,i),centers(2,:,i),normals(i,:));
            if flip % centers3 is point 1, centers1 is point 2, centers2 is point 3
                orders = [3,1,2];
                dis = [dis_local(3);dis_local(1);dis_local(2)];
            else % centers1 is point 1, centers3 is point 2, centers2 is point 3
                orders = [1,3,2];
                dis = [dis_local(3);dis_local(2);dis_local(1)];
            end
        end
        [dis_local,~] = sort(dis_local);
        if dis_local(1) > 61 | dis_local(1) < 59 | dis_local(2) > 77 | dis_local < 75 | dis_local(3) > 77 | dis_local(3) < 75
            % line 12 is around 60 (if line 12 is x axis), line 13 and line
            % 23 should have length around 76 in my case.
            fprintf("Point Cloud %d error, please Modify it Manually or remove it from dataset...\n",k)
        end
    end
end

function T = xyzO2T(x,y,z,O)
    if isequal(size(x),[3,1])
        T = [x,y,z,O;0,0,0,1];
    else
        T = [x',y',z',O';0,0,0,1];
    end
end

function [T,flip] = xyz2T(p1,p2,p3,normal)
    vec12 = p2 - p1;
    vec12 = vec12/norm(vec12);
    vec13 = p3 - p1;
    O = p1 + vec12 *(vec13*vec12');
    x = -vec12; y = p3-O; y = y/norm(y);
    z = cross(x,y);
    z = z/norm(z);
    if z*double(normal)' < 0
        % use facial normal to decide the direction of Z axis.
        z = -z;
        x = -x;
        flip = true;
    else
        flip = false;
    end
    T = xyzO2T(x,y,z,O);
end