clc;
clear;
close all;
addpath(genpath('.'));

NAME = "Validation_New";
ALGS = ["Wang","G3","Wu","Fu","Ma","GA"];
%% Load ABB DATA
abb_file = fopen(sprintf("./traj/%s/ABBTraj.txt",NAME),'r');
formatSpec = '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n';
ABB_Raw = fscanf(abb_file,formatSpec,[13,Inf]);
fclose(abb_file);
Ts_ABB = convert_ABB_Pose(ABB_Raw(7:13,:));
%% Load JAKA DATA
jaka_file = fopen(sprintf("./traj/%s/JAKATraj.txt",NAME),'r');
formatSpec = '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n';
JAKA_Raw = fscanf(abb_file,formatSpec,[12,Inf]);
fclose(abb_file);
Ts_JAKA = convert_JAKA_Pose(JAKA_Raw(1:6,:));
%% Load XYZ and derive the tfmatrix
load("./result/XYZ/XYZ_ALL.mat")
N = size(Ts_JAKA,3);
N_alg = 6;
Tfs = zeros(4,4,N,N_alg);
for i = 1:N
    Tfs(:,:,i,1) = (Ym*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*Xm);
    Tfs(:,:,i,2) = (Yg3*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*Xg3);
    Tfs(:,:,i,3) = (Yl*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*Xl);
    Tfs(:,:,i,4) = (YFu*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*XFu);
    Tfs(:,:,i,5) = (YMa*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*XMa);
    Tfs(:,:,i,6) = (YGA*Ts_ABB(:,:,i))\(Ts_JAKA(:,:,i)*XGA);
end
%% Parameters
N_Pose = 8;
N_Shots = ones(1,N_Pose)*5;
assert(size(N_Shots,2)==N_Pose,"Not Enough Shots are provided.")
%% Load Filtered Point Cloud
figsave = false;
radius = zeros(3,N_Pose,N_alg);
idx = cumsum(N_Shots);
idx = [0,idx(1:end-1)]+1;
for s = 1:N_alg
    parfor i = 1:N_Pose
        c = idx(i);
        if figsave
            ax = figure("Visible","off");
            hold on
        end
        for k = 1:3
            counter = c;
            ptball = pcread(sprintf("./meshes/%s/balls/%d_%d.ply",NAME,counter,k));
            ptball = pctransform(ptball,rigid3d(Tfs(1:3,1:3,counter,s)',Tfs(1:3,4,counter,s)'));
            locations = ptball.Location;
            coutner = counter + 1;
            for j = 2:N_Shots(i)
                 ptC = pcread(sprintf("./meshes/%s/balls/%d_%d.ply",NAME,counter,k));
                 ptC = pctransform(ptC,rigid3d(Tfs(1:3,1:3,counter,s)',Tfs(1:3,4,counter,s)'));
                 locations = cat(1,locations,ptC.Location);
                 counter = counter + 1;
            end
            maxDistance = 100;
            ptball = pointCloud(locations);
            if figsave
                pcshow(select(ptball,1:100:ptball.Count));
%                 pcshow(ptball);
            end
    %         balls = pcfitsphere(ptball,maxDistance,Confidence=99.9);
            [C,R] = sphereFit(ptball.Location);
            radius(k,i,s) = R;
        end
        if figsave
            savefig(ax,sprintf(".\\result\\%s\\match_pose_balls_%s_%d",NAME,ALGS(s),i))
            saveas(ax,sprintf(".\\result\\%s\\match_pose_balls_%s_%d.png",NAME,ALGS(s),i))
            fprintf("%s pose %d finished\n",ALGS(s),i)
            close(ax);
        end
    end
end
%% RESULT
R = [30.0060,30.0075,30.0085]'/2;
err = zeros(N_alg,N_Pose);
for i = 1:N_alg
    err(i,:) = sum(abs(radius(:,:,i)-R),1)/3;
end
% err1(3) = []; err2(3) = [];
fprintf("Wang's Method has error: %f,\t with var %f, \t max %f\n" ,mean(err(1,:)),var(err(1,:)),max(err(1,:)))
fprintf("G3's Method has error: %f,\t with var %f, \t max %f\n",mean(err(2,:)),var(err(2,:)),max(err(2,:)))
fprintf("Liao's Method has error: %f,\t with var %f, \t max %f\n" ,mean(err(3,:)),var(err(3,:)),max(err(3,:)))
fprintf("Fu's Method has error: %f,\t with var %f, \t max %f\n" ,mean(err(4,:)),var(err(4,:)),max(err(4,:)))
fprintf("Ma's Method has error: %f,\t with var %f, \t max %f\n" ,mean(err(5,:)),var(err(5,:)),max(err(5,:)))
% fprintf("PGA's Method has error: %f,\t with var %f, \t max %f\n" ,mean(err(6,:)),var(err(6,:)),max(err(6,:)))