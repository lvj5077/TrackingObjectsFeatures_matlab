clear
close all
% clc
addpath('/Users/jin/Q_Mac/mexopencv');
%%
testNum = 200;
stdIdx = 1;
gt_d = [    3.3100
    6.3700
    9.4400
   12.5000
   15.4400
   18.6300
   21.6900
   24.7500
   27.8200
   30.7600]';
secNum = length(gt_d);
results = zeros(testNum,secNum,3);
resultsR = zeros(testNum,secNum);
gt = [zeros(1,secNum);zeros(1,secNum);gt_d];
% secNum = secNum+1;


for testId = 1:testNum
    imgId = 200+testId;
    data_path = "/Volumes/BlackSSD/rotateIP12/rpy/x/rot_roll_";
    I1 = imread(data_path+num2str(stdIdx)+"/color/"+num2str(imgId)+".png");
    I1 = rgb2gray(I1);

    I_dpt = imread(data_path+num2str(stdIdx)+"/depth/"+num2str(imgId)+".png");

    Ipre = I1;
    Inxt = I1;

    tempCell = cv.goodFeaturesToTrack(Ipre, 'MaxCorners', 801, 'QualityLevel', 0.01, 'MinDistance', 20);
    prevPts = double(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
    firstPts = prevPts;
    nextPts = prevPts; 
    
    for i = 1:secNum
%         K = [1460.110474,0,956.812561;0,1460.110474,652.749084;0,0,1];
        imgId = 200+testId;
        I_cur = imread(data_path+num2str(i+stdIdx)+"/color/"+num2str(imgId)+".png");
%         I_cur = imread(color_path+num2str(i)+".png");
        I_cur = rgb2gray(I_cur);

        Ipre = Inxt;
        Inxt = I_cur;

        prevPts = nextPts;

        tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
        nextPts = double(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
        % C2 = double(reshape(tempMat',[2,length(tempMat)/2]))';
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
    for i = 1:secNum
        imgId = 200+testId;
%         I_cur = imread(color_path+num2str(i)+".png");
        I_cur = imread(data_path+num2str(i+stdIdx)+"/color/"+num2str(imgId)+".png");
        I_cur = rgb2gray(I_cur);

        
        
        Ipre = Inxt;
        Inxt = I_cur;

        prevPts = checkPts;

        tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
        nextPts = double(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
        % C2 = double(reshape(tempMat',[2,length(tempMat)/2]))';
        idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
        nextPts = nextPts(idx,:);
        prevPts = prevPts(idx,:);

        checkPts = nextPts;

%         [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
%         nextPts = nextPts(mask==1,:);
%         prevPts = prevPts(mask==1,:);
%         h = figure;
%         showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
    %     outDir = data_path+"summary/"+num2str(i-1)+"to"+num2str(i)+".png";
    %     exportgraphics(h,outDir) 

        prevPts = firstPts;
        nextPts = checkPts;
        [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',1);
        nextPts = nextPts(mask==1,:);
        prevPts = prevPts(mask==1,:);
%         close all
%         h = figure;
%         showMatchedFeatures(I1,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
    %     outDir = data_path+"summary/1"+"to"+num2str(i)+".png";
    %     exportgraphics(h,outDir) 
        imgId = 200+testId;
        dd = importdata(data_path+num2str(1)+"/Frames.txt");
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K1 = [fx,0,cx;0,fy,cy;0,0,1];
        
%         dd = importdata(data_path+num2str(i)+"/Frames.txt");
%         fx = dd(imgId,3);
%         fy = dd(imgId,4);
%         cx = dd(imgId,5);
%         cy = dd(imgId,6);
%         K2 = [fx,0,cx;0,fy,cy;0,0,1];
        mvImg = [];
        fixObj = [];
        for j = 1:length(prevPts)
            ui = double(prevPts(j,1));
            vi = double(prevPts(j,2));
            xi = round(ui/7.5);
            yi = round(vi/7.5);
            xi = max(xi,1);
            yi = max(yi,1);
%             I_resize = imresize(I_dpt,[1440,1920]);
% 
%             d = double(I_resize(vi, ui)) * 0.001;
%             closeDepth = 1;
            
        	d = double(I_dpt(yi, xi)) * 0.001;
            closeDepth = 0;

%             dist_x = ui - (7.5*xi-3.75);
%             dist_y = vi - (7.5*yi-3.75);
%             dis_thred = 2;

            dist_x = xi - ui/7.5;
            dist_y = yi - vi/7.5;
            dis_thred = 0.4;
            if (abs(dist_x)<dis_thred && abs(dist_y)<dis_thred)
                closeDepth = 1;
            end


            if(d>= 0.1 && d <= 5 && closeDepth>0) 
                pt_x = d*(ui-cx)/fx;
                pt_y = d*(vi-cy)/fy;
                pt = [pt_x,pt_y,d];
                fixObj = [pt;fixObj];
                mvImg = [double(nextPts(j,:));mvImg];
            end
        end
        [rvec, tvec, success, inliers] = cv.solvePnPRansac(fixObj, mvImg, K1);
%         K = (K1+K2)/2;
%         E = K1' * F * K1;
%         [R, t, good, mask, triangulatedPoints] = cv.recoverPose(E,prevPts,nextPts);

%         tvec'
        rotm = cv.Rodrigues(rvec);
        results(testId,i,:) = rotm2eul(rotm)*180/pi;
        axang = rotm2axang(rotm);
        resultsR(testId,i) = axang(4)*180/pi;
    end
end
%%

fff=abs(resultsR(:,:)-gt_d.*ones(testNum,1));
[x,y] = find(fff>0.5);
length(x)
resultsR(x,:)=[];
results(x,:,:)=[];


Meanresults = zeros(secNum,3);
MeanError = zeros(secNum,3);
STDresults = zeros(secNum,3);
MeanresultsM = zeros(secNum,1);
MeanresultsMm = zeros(secNum,1);
MeanresultsMstd = zeros(secNum,1);
for i = 1:secNum
    Meanresults(i,1) = mean( ( results(:,i,1) )  );
    Meanresults(i,2) = mean( ( results(:,i,2)  )  );
    Meanresults(i,3) = mean( ( results(:,i,3)  )  );
    
    MeanresultsM(i) = mean( resultsR(:,i) );
    MeanresultsMm(i) = mean( abs( resultsR(:,i)-gt_d(i) ));
    MeanresultsMstd(i) = std(  resultsR(:,i)-gt_d(i) );

    MeanError(i,1) = mean( abs( results(:,i,1)-gt(1,i)  )  );
    MeanError(i,2) = mean( abs( results(:,i,2)-gt(2,i)  )  );
    MeanError(i,3) = mean( abs( results(:,i,3)-gt(3,i)  )  );
    STDresults(i,:) = reshape(std(results(:,i,:)),[1,3]);
    
    figure,plot(resultsR(:,i))
end
Meanresults
MeanError
STDresults
MeanresultsM
MeanresultsMm
MeanresultsMstd

figure,
subplot(3,1,1)
plot(results(:,5,1))
subplot(3,1,2)
plot(results(:,5,2))
subplot(3,1,3)
plot(results(:,5,3))