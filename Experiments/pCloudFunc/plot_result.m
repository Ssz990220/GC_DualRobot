function [ax,ax2] = plot_result(Ts,ball1, ball2, ball3, display,normview)
    if nargin == 4
        ax = figure();
    else
        if display
            ax = figure();
        else
            ax = figure('visible','off');
        end
    end
    balls = pointCloud([ball1.Location;ball2.Location;ball3.Location]);
    color = [ones(1,ball1.Count),ones(1,ball2.Count)*2,ones(1,ball3.Count)*3];
    color = color(1:100:balls.Count);
    balls = select(balls,1:100:balls.Count);
    pcshow(balls.Location,color);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    if normview
        view(Ts(1:3,3));
    else
        view(30,60);
    end
    hold on
    % plot3(centers(1,1),centers(1,2),centers(1,3),"-o","MarkerSize",15)
    % plot3(centers(2,1),centers(2,2),centers(2,3),"-o","MarkerSize",15)
    % plot3(centers(3,1),centers(3,2),centers(3,3),"-o","MarkerSize",15)
    scale = 10;
    quiver3(Ts(1,4),Ts(2,4),Ts(3,4),Ts(1,1)*scale,Ts(2,1)*scale,Ts(3,1)*scale,'-b','LineWidth',2)
    quiver3(Ts(1,4),Ts(2,4),Ts(3,4),Ts(1,2)*scale,Ts(2,2)*scale,Ts(3,2)*scale,"-g",'LineWidth',2)
    quiver3(Ts(1,4),Ts(2,4),Ts(3,4),Ts(1,3)*scale,Ts(2,3)*scale,Ts(3,3)*scale,"-r",'LineWidth',2)

    if nargin == 4
        ax2 = figure();
    else
        if display
            ax2 = figure();
        else
            ax2 = figure('visible','off');
        end
    end
    subplot(1,3,1); pcshow(select(ball1, 1:100:ball1.Count));
    subplot(1,3,2); pcshow(select(ball2, 1:100:ball2.Count));
    subplot(1,3,3); pcshow(select(ball3, 1:100:ball3.Count));
end