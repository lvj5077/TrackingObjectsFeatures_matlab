clear
close all
% clc
addpath('/Users/jin/Q_Mac/mexopencv');
%%
testNum = 200;
gt = [3.1300    6.1900    9.2500   12.3200   15.3200   18.4400];
secNum = length(gt);
results = zeros(testNum,secNum,3);
gt = [zeros(1,secNum);gt;zeros(1,secNum)];
secNum = secNum+1;


for testId = 1:testNum
    imgId = 200+testId;
    data_path = "/Volumes/BlackSSD/rotateIP12/rpy/y/rot_pitch_";
    I1 = imread(data_path+"1/color/"+num2str(imgId)+".png");
    I1 = rgb2gray(I1);


    Ipre = I1;
    Inxt = I1;

    tempCell = cv.goodFeaturesToTrack(Ipre, 'MaxCorners', 801, 'QualityLevel', 0.01, 'MinDistance', 20);
    prevPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
    firstPts = prevPts;
    nextPts = prevPts; 
    
    for i = 2:secNum
%         K = [1460.110474,0,956.812561;0,1460.110474,652.749084;0,0,1];
        imgId = 200+testId;
        I_cur = imread(data_path+num2str(i)+"/color/"+num2str(imgId)+".png");
%         I_cur = imread(color_path+num2str(i)+".png");
        I_cur = rgb2gray(I_cur);

        Ipre = Inxt;
        Inxt = I_cur;

        prevPts = nextPts;

        tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
        nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
        % C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
        imgId = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
        nextPts = nextPts(imgId,:);
        prevPts = prevPts(imgId,:);

        [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
        nextPts = nextPts(mask==1,:);
        prevPts = prevPts(mask==1,:);
        firstPts = firstPts(imgId,:);
        firstPts = firstPts(mask==1,:);
    %     figure,showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage');
    end
    prevPts = firstPts;
    nextPts = prevPts; 
    checkPts = prevPts;
    Inxt = I1;
    for i = 2:secNum
        imgId = 200+testId;
%         I_cur = imread(color_path+num2str(i)+".png");
        I_cur = imread(data_path+num2str(i)+"/color/"+num2str(imgId)+".png");
        I_cur = rgb2gray(I_cur);

        Ipre = Inxt;
        Inxt = I_cur;

        prevPts = checkPts;

        tempCell = cv.calcOpticalFlowPyrLK(Ipre, Inxt, prevPts);
        nextPts = single(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
        % C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
        imgId = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
        nextPts = nextPts(imgId,:);
        prevPts = prevPts(imgId,:);

        checkPts = nextPts;

%         [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
%         nextPts = nextPts(mask==1,:);
%         prevPts = prevPts(mask==1,:);
    %     h = figure;
    %     showMatchedFeatures(Ipre,Inxt,prevPts,nextPts,'montage','PlotOptions',{'g.','r.','y-'});
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
        
        dd = importdata(data_path+num2str(i)+"/Frames.txt");
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K2 = [fx,0,cx;0,fy,cy;0,0,1];
        E = K2' * F * K1;
        [R, t, good, mask, triangulatedPoints] = cv.recoverPose(E,prevPts,nextPts);
        results(testId,i-1,:) = rotm2eul(R)*180/pi;
    end
end
%%
fff=abs(results(:,:,3));
[x,y] = find(fff>100);
results(x,:,:)=[];
Meanresults = zeros(secNum-1,3);
MeanError = zeros(secNum-1,3);
STDresults = zeros(secNum-1,3);

for i = 1:secNum-1
    Meanresults(i,1) = mean( ( results(:,i,1) )  );
    Meanresults(i,2) = mean( ( results(:,i,2)  )  );
    Meanresults(i,3) = mean( ( results(:,i,3)  )  );
    
    MeanError(i,1) = mean( abs( results(:,i,1)-gt(1,i)  )  );
    MeanError(i,2) = mean( abs( results(:,i,2)-gt(2,i)  )  );
    MeanError(i,3) = mean( abs( results(:,i,3)-gt(3,i)  )  );
    STDresults(i,:) = reshape(std(results(:,i,:)),[1,3]);
end
Meanresults
MeanError
STDresults
%%
figure,
subplot(3,1,1)
plot(results(:,2,1))
subplot(3,1,2)
plot(results(:,2,2))
subplot(3,1,3)
plot(results(:,2,3))