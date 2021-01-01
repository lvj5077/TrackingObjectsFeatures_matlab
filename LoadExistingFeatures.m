load('obj_stat.mat')
load('vpts.mat')
I1_pdic = imread('/Users/jin/Desktop/results/predictions125.jpg');
I2_pdic = imread('/Users/jin/Desktop/results/predictions530.jpg');
colorstyle = ["g*","b*","c*","m*","r*","mx","b<","r<"];

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
