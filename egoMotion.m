clear
close all
clc
addpath('/Users/jin/Q_Mac/mexopencv');
I1 = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/1.png');
I1 = rgb2gray(I1);
I2 = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/2.png');
I2 = rgb2gray(I2);


Ipre = I1;
Inxt = I2;

tempCell = cv.goodFeaturesToTrack(Ipre, 'MaxCorners', 800, 'QualityLevel', 0.01, 'MinDistance', 20);
prevPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';

tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
% C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
nextPts = nextPts(idx,:);
prevPts = prevPts(idx,:);
[F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
nextPts = nextPts(mask==1,:);
prevPts = prevPts(mask==1,:);
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');

firstPts = prevPts;
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/3.png');
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
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/4.png');
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

figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/5.png');
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
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');

firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/6.png');
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
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');

firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/7.png');
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
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');

firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);
%%
I_first = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/1.png');
I_first = rgb2gray(I_first);

Ipre = I_first;
% Inxt = I_cur;
% % 
% % prevPts = nextPts;
% % 
% % tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
% % nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
% % % C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
% % idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
% % nextPts = nextPts(idx,:);
% % prevPts = prevPts(idx,:);
[F, mask] = cv.findFundamentalMat(firstPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
firstPts = firstPts(mask==1,:);
nextPts = nextPts(mask==1,:);
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Inxt,Ipre,nextPts,firstPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'7to1.png') 

%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/2.png');
I_cur = rgb2gray(I_cur);

Ipre = I_first;
Inxt = I_cur;

prevPts = firstPts;

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
figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'1to2.png') 
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/3.png');
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
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'2to3.png') 
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/4.png');
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

h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'3to4.png') 
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/5.png');
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
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'4to5.png') 

firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/6.png');
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
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'5to6.png') 

firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);
%%
I_cur = imread('/Volumes/BlackSSD/rotateIP12/rpy/y/summary/color/7.png');
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
h = figure('units','normalized','outerposition',[0 0 1 1]);
showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
exportgraphics(h,'6to7.png') 
firstPts = firstPts(idx,:);
firstPts = firstPts(mask==1,:);