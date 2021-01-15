clear
close all
clc
addpath('/Users/jin/Q_Mac/mexopencv');
%%
data_path = "/Volumes/BlackSSD/rotateIP12/rpy/y/";
color_path = data_path+"summary/color/";
I1 = imread(color_path+num2str(1)+".png");
I1 = rgb2gray(I1);

K = [1460.110474,0,956.812561;0,1460.110474,652.749084;0,0,1];
secNum = 7;

Ipre = I1;
Inxt = I1;

tempCell = cv.goodFeaturesToTrack(Ipre, 'MaxCorners', 800, 'QualityLevel', 0.01, 'MinDistance', 20);
prevPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
firstPts = prevPts;
nextPts = prevPts; 
for i = 2:secNum
    I_cur = imread(color_path+num2str(i)+".png");
    I_cur = rgb2gray(I_cur);

    Ipre = Inxt;
    Inxt = I_cur;

    prevPts = nextPts;

    tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
    nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
    % C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
    idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
    nextPts = nextPts(idx,:);
    prevPts = prevPts(idx,:);

    [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
    nextPts = nextPts(mask==1,:);
    prevPts = prevPts(mask==1,:);
    firstPts = firstPts(idx,:);
    firstPts = firstPts(mask==1,:);
%     figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');
end
prevPts = firstPts;
nextPts = prevPts; 
checkPts = prevPts;
Inxt = I1;
for i = 2:secNum
    I_cur = imread(color_path+num2str(i)+".png");
    I_cur = rgb2gray(I_cur);

    Ipre = Inxt;
    Inxt = I_cur;

    prevPts = checkPts;

    tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
    nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
    % C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
    idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
    nextPts = nextPts(idx,:);
    prevPts = prevPts(idx,:);
    
    checkPts = nextPts;

    [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
    nextPts = nextPts(mask==1,:);
    prevPts = prevPts(mask==1,:);
%     h = figure;
%     showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
%     outDir = data_path+"summary/"+num2str(i-1)+"to"+num2str(i)+".png";
%     exportgraphics(h,outDir) 
    
    prevPts = firstPts;
    nextPts = checkPts;
    [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
    nextPts = nextPts(mask==1,:);
    prevPts = prevPts(mask==1,:);
    close all
    h = figure;
    showMatchedFeatures(I1,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
%     outDir = data_path+"summary/1"+"to"+num2str(i)+".png";
%     exportgraphics(h,outDir) 
    E = K' * F * K;
    [R, t, good, mask, triangulatedPoints] = cv.recoverPose(E,prevPts,nextPts);
    rotm2eul(R)*180/pi
end