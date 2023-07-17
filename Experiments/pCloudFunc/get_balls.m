function [balls,normal] = get_balls(ptCloud)
%% Filter out the big surface at the bottom of the balls.
maxDistance = 0.5;
[model1,inlierIndices,outlierIndices,meanError] = pcfitplane(ptCloud,...
        maxDistance,'Confidence',99,'MaxNumTrials',1000);
[warnMsg, ~] = lastwarn;
if ~isempty(warnMsg) & meanError > 0.2
    % if meanError is less than 0.2, most likely it has found the plane but
    % either it's so large to find the best fit parameter.
    balls = ptCloud;
    normal = [0,0,1]; % always face the camera so the default normal will work as well
    lastwarn('')
else
    ptCloud_noPlane = select(ptCloud,outlierIndices);
    Locations = [ptCloud_noPlane.Location,ones(size(ptCloud_noPlane.Location,1),1)];
    indicator = Locations*model1.Parameters';
    il0_Points = indicator<-5;                  %% index of points on one side, but filter out what remains by the plane
    ig0_Points = indicator>5;                  %% index of points on other side, but filter out what remains by the plane
    nl0_Points = sum(il0_Points);                                   %% Number of points on one side
    ng0_Points = sum(ig0_Points);                                   %% Number of points on other side
    index = 1.0:1.0:ptCloud_noPlane.Count;
    if nl0_Points > ng0_Points                  % select the side with more points
        balls = select(ptCloud_noPlane,index(il0_Points));
    else
        balls = select(ptCloud_noPlane,index(ig0_Points));
    end
    mballs = mean(balls.Location);
    mplanes = mean(select(ptCloud,inlierIndices).Location);
    if (mballs-mplanes)*model1.Normal'>0        % set normal in right direction (faces "up" )
        normal = model1.Normal;
    else
        normal = -model1.Normal;
    end
end

end
