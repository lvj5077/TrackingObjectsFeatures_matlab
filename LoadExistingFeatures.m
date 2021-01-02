clear
close all
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
    else
        obj_stat(i,1) = 1;
        ccstyle = colorstyle( obj_stat(i,1) );
    end
    hold on
    plot(vpts1(i,1),vpts1(i,2),ccstyle);
    hold on
    plot(vpts2(i,1)+640,vpts2(i,2),ccstyle);
end

%%
clc
% [T,inlierPts,inlierObjs] = RANSACpose2d2d_obj(vpts1,vpts2,ones(length(vpts1),1));
[T,inlierPts,inlierObjs] = RANSACpose2d2d_obj(vpts1,vpts2,obj_stat);
inlierObjs'

error_Pts = sampsonErrf(T, vpts1, vpts2);
[idx1,~] = find(inlierPts==1);
[idx2,val] = find(error_Pts'<3e-4);
idx = idx1(inlierObjs(obj_stat(idx1)) >0);
%%
% [idx,~]=find(obj_stat==5);
% test1 = vpts1(idx,:);
% test2 = vpts2(idx,:);
% test_error = sampsonErrf(T, test1, test2);
% [testIdx,~] = find(test_error'<3e-4);
% rr = length(testIdx)/length(test1)
% testIdxinlier = idx(testIdx);
% 
% idx = testIdxinlier;
% [idx,~] = find(obj_stat(idx2)>0);
vpts1_inlier = vpts1(idx,:);%(inlierPts==1&error_Pts(:)<5e-4,:);
vpts2_inlier = vpts2(idx,:);%(inlierPts==1&error_Pts(:)<5e-4,:);
obj_stat_inler = obj_stat(idx,:);%(inlierPts==1&error_Pts(:)<5e-4,:);%(inlierPts==1&obj_stat(:)~=5&error_Pts(:)<0.2,:);
figure,showMatchedFeatures(I1_pdic,I2_pdic,vpts1_inlier,vpts2_inlier,'montage','PlotOptions',{'yo','yo','y-'});

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