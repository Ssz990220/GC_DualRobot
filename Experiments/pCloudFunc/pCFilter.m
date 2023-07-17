function [Centers,normal,T,R,dis] = pCFilter(ptCloud,k,cfg)
%PCFILTER 处理三球点云以获取球心位置
%% delete tha big plane (if there is one)
[balls,normal] = get_balls(ptCloud);
%% Filters--Find Target Balls by looking for the biggest three clusters
[ball_labels,labels] = cluster_balls(balls,k);
%% Sphere Registry--Ransac
[ball_collection, ib_collection,Centers,model1,model2,model3] = fit_balls(balls,ball_labels,labels);
R = [model1.Radius, model2.Radius, model3.Radius]';
%% get T
[T,order,dis] = center2T(Centers,normal,k);
%% reorder Balls
ball1 = ball_collection{order(1)};ball2 = ball_collection{order(2)};ball3 = ball_collection{order(3)};
ib1 = ib_collection{order(1)};ib2 = ib_collection{order(2)};ib3 = ib_collection{order(3)};
R = [R(order(1));R(order(2));R(order(3))];
%% save
if isfield(cfg,"save_ball")
    if cfg.save_ball
        pcwrite(select(ball1,ib1),sprintf('.\\meshes\\%s\\balls\\%d_1.ply',cfg.name,k));
        pcwrite(select(ball2,ib2),sprintf('.\\meshes\\%s\\balls\\%d_2.ply',cfg.name,k));
        pcwrite(select(ball3,ib3),sprintf('.\\meshes\\%s\\balls\\%d_3.ply',cfg.name,k));
    end
end
if cfg.savefig
    [ax,ax2] = plot_result(T,ball1, ball2, ball3, cfg.show,cfg.normview);
    savefig(ax,sprintf(".\\result\\%s\\%d.fig",cfg.name,k))
    saveas(ax,sprintf(".\\result\\%s\\%d.png",cfg.name,k))
    if ~cfg.show
        close(ax)
    end
    savefig(ax2,sprintf(".\\result\\%s\\%d-3b.fig",cfg.name,k))
    saveas(ax2,sprintf(".\\result\\%s\\%d-3b.png",cfg.name,k))
    if ~cfg.show
        close(ax2)
    end
end
end

