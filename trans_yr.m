clear
close all
format short g
clc
addpath('/Users/jin/Q_Mac/mexopencv');

testNum = 250;
stdIdx = 10;
gt_d = -10*[    10
    20
    30
   40
   50]';
secNum = length(gt_d);
results = zeros(testNum,secNum,4);
gt = [gt_d;zeros(1,secNum);zeros(1,secNum)];
% secNum = secNum+1;

%%
for testId = 1:testNum
%     testId
    imgId = 250+testId;
    data_path = "/Volumes/BlackSSD/rotateIP12/trans_y/trans_y_";
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
        imgId = 250+testId;
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
        imgId = 250+testId;
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
        imgId = 250+testId;
        dd = importdata(data_path+num2str(i)+"/Frames.txt");
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K2 = [fx,0,cx;0,fy,cy;0,0,1];
        
        dd = importdata(data_path+num2str(1)+"/Frames.txt");
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K1 = [fx,0,cx;0,fy,cy;0,0,1];
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
        
        
        [rvec, tvec, success, inliers] = cv.solvePnPRansac(fixObj, mvImg, K2);
%         K = (K1+K2)/2;
%         E = K2' * F * K1;
%         [R, t, good, mask, triangulatedPoints] = cv.recoverPose(E,prevPts,nextPts);

%         tvec'
        results(testId,i,1:3) = 1000*tvec';
        results(testId,i,4) = 1000*norm(tvec);
    end
end

%%
results_bkp = results;

%%
results = results_bkp;
%%
% results = importdata("trans_y_pnp.csv");
% results = reshape(results,[250,5,4]);
% fff=abs(results(:,:,3));
% [x,y] = find(fff>100);
% results(x,:,:)=[];
Meanresults = zeros(secNum,4);
MeanError = zeros(secNum,4);
STDresults = zeros(secNum,4);
MeanresultsM = zeros(secNum,1);
MeanresultsMm = zeros(secNum,1);
MeanresultsMstd = zeros(secNum,1);
h = figure('units','normalized','outerposition',[0 0 1 1]);
t = tiledlayout(1,secNum,'TileSpacing','Compact','units','normalized','outerposition',[0 0 1 1]);
for i = 1:secNum
    Meanresults(i,1) = mean( ( results(:,i,1) )  );
    Meanresults(i,2) = mean( ( results(:,i,2)  )  );
    Meanresults(i,3) = mean( ( results(:,i,3)  )  );
    Meanresults(i,4) = mean( ( results(:,i,4)  )  );

    MeanError(i,1) = mean( abs( results(:,i,1)-gt(1,i)  )  );
    MeanError(i,2) = mean( abs( results(:,i,2)-gt(2,i)  )  );
    MeanError(i,3) = mean( abs( results(:,i,3)-gt(3,i)  )  );
    MeanError(i,4) = mean( abs( results(:,i,4)-abs(gt_d(i))  )  );
    STDresults(i,:) = reshape(std(results(:,i,:)),[1,4]);
    
    nexttile
%     subplot(1,secNum,i)
    pt1 = plot(results(:,i,4),'b.','Markersize',10);
    hold on
    pt3 = plot(abs(gt_d(i))*ones(1,length(results)),'r','Linewidth',5);
    hold on
    pt2 = plot(Meanresults(i,4)*ones(1,length(results)),'g','Linewidth',5);
    hold on
    pt4 = plot((Meanresults(i,4)-2*STDresults(i,4))*ones(1,length(results)),'k--','Linewidth',3);
    pt4 = plot((Meanresults(i,4)+2*STDresults(i,4))*ones(1,length(results)),'k--','Linewidth',3);
    ylim([Meanresults(i,4)-20,Meanresults(i,4)+20])
    grid minor
    set(gca,'FontSize',20)
    set(gcf,'color','w');
    if i ==1
        legend([pt1 pt2 pt3 pt4],{'Measurement','Mean','Groundtruth','2-sigma'},'FontSize',10)
    end
    
end
title(t,'translation Y-axis','FontSize',40)
xlabel(t,'trail number','FontSize',30)
ylabel(t,'Measurement (mm)','FontSize',30)

exportgraphics(h,'/Users/jin/Library/Mobile Documents/com~apple~CloudDocs/iphone_ego/trans_x.png') 
csvwrite("/Users/jin/Library/Mobile Documents/com~apple~CloudDocs/iphone_ego/trans_x_pnp.csv",results)

Meanresults
MeanError
STDresults