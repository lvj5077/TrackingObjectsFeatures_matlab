load('obj_stat.mat')
load('vpts.mat')
I1_pdic = imread('/Users/jin/Desktop/results/predictions125.jpg');
I2_pdic = imread('/Users/jin/Desktop/results/predictions530.jpg');
colorstyle = ["g*","b*","c*","m*","r*","mx","b<","r<"];
%%
close all
figure,showMatchedFeatures(I1_pdic,I2_pdic,vpts1(obj_stat>0,:),vpts2(obj_stat>0,:),'montage','PlotOptions',{'yo','yo','y-'});
for i = 1:length(vpts1)
    ccstyle = "r<";
    if (obj_stat(i,1)>0)
        ccstyle = colorstyle( obj_stat(i,1) );
    end
    hold on
    plot(vpts1(i,1),vpts1(i,2),ccstyle);
    hold on
    plot(vpts2(i,1)+640,vpts2(i,2),ccstyle);
end

%%
clc
% [T,inlierPts,inlierObjs] = RANSACpose2d2d_obj(vpts1,vpts2,ones(length(vpts1),1))
[T,inlierPts,inlierObjs] = RANSACpose2d2d_obj(vpts1,vpts2,obj_stat);
inlierObjs'
vpts1_inlier = vpts1(inlierPts==1,:);
vpts2_inlier = vpts2(inlierPts==1,:);
obj_stat_inler = obj_stat(inlierPts==1,:);%(inlierPts==1&obj_stat(:)~=5&error_Pts(:)<0.2,:);
figure,showMatchedFeatures(I1_pdic,I2_pdic,vpts1_inlier(obj_stat_inler>0,:),vpts2_inlier(obj_stat_inler>0,:),'montage','PlotOptions',{'yo','yo','y-'});

% obj_stat(obj_stat(:,1)==0) = 8;
for i = 1:length(vpts1_inlier)
    ccstyle = "r<";
    if (obj_stat_inler(i,1)>0)
        ccstyle = colorstyle( obj_stat_inler(i,1) );
    end
    hold on
    plot(vpts1_inlier(i,1),vpts1_inlier(i,2),ccstyle);
    hold on
    plot(vpts2_inlier(i,1)+640,vpts2_inlier(i,2),ccstyle);
end