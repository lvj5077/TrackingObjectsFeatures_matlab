clear
close all
format short g
clc
addpath('/Users/jin/Q_Mac/mexopencv');
%%
testNum = 250;
stdIdx = 1;
gt_d = [3.1300    6.1900    9.2500   12.3200   15.3200   18.4400];
secNum = length(gt_d);
results = zeros(testNum,secNum,4);
results_pnp = zeros(testNum,secNum,4);
gt = [zeros(1,secNum);gt_d;zeros(1,secNum)];

for testId = 1:testNum
%     testId
    imgId = 250+testId;
    data_path = "/Volumes/BlackSSD/rotateIP12/rpy/y/rot_pitch_";
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
        K = (K1+K2)/2;
        E = K2' * F * K1;
        [R, t, good, mask, triangulatedPoints] = cv.recoverPose(E,prevPts,nextPts);
        results(testId,i,1:3) = rotm2eul(R)*180/pi;
        axang = rotm2axang(R);
        results(testId,i,4) = axang(4)*180/pi;
        
        
%         tvec'
        rotm = cv.Rodrigues(rvec);
        results_pnp(testId,i,1:3) = rotm2eul(rotm)*180/pi;
        axang = rotm2axang(rotm);
        results_pnp(testId,i,4) = axang(4)*180/pi;
    end
end

%%
results_pnp_bkp = results_pnp;
results_bkp = results;
csvwrite("/Users/jin/Library/Mobile Documents/com~apple~CloudDocs/iphone_ego/rot_y_8p.csv",results)
csvwrite("/Users/jin/Library/Mobile Documents/com~apple~CloudDocs/iphone_ego/rot_y_pnp.csv",results_pnp)

%%
results = results_bkp;
results_pnp = results_pnp_bkp;
outlierThes = 0.5;
sigmaTimes = 5;
Meanresults = zeros(secNum,4);
MeanError = zeros(secNum,4);
STDresults = zeros(secNum,4);
h = figure('units','normalized','outerposition',[0 0 1 1]);
t = tiledlayout(2,secNum,'TileSpacing','Compact','units','normalized','outerposition',[0 0 1 1]);
for i = 1:secNum

    Meanresults(i,4) = mean( results(:,i,4) );
    MeanError(i,4) = mean( abs( results(:,i,4)-gt_d(i) ));
    STDresults(i,:) = reshape(std(results(:,i,:)),[1,4]);
    
    checkOutlier=abs(results(:,i,4)-median( results(:,i,4) ).*ones(length(results),1));
    [x,y] = find(checkOutlier>min(1,max(outlierThes,STDresults(i,4)*sigmaTimes)));
%     [x,y] = find(checkOutlier>STDresults(i,4)*sigmaTimes|checkOutlier>0.5);
    results(x,:,:)=[];
%     [x,y] = find(checkOutlier>1);
%     results(x,:,:)=[];
    

    Meanresults(i,1) = mean( ( results(:,i,1) )  );
    Meanresults(i,2) = mean( ( results(:,i,2)  )  );
    Meanresults(i,3) = mean( ( results(:,i,3)  )  );

    MeanError(i,1) = mean( abs( results(:,i,1)-gt(1,i)  )  );
    MeanError(i,2) = mean( abs( results(:,i,2)-gt(2,i)  )  );
    MeanError(i,3) = mean( abs( results(:,i,3)-gt(3,i)  )  );
    
    Meanresults(i,4) = mean( results(:,i,4) );
    MeanError(i,4) = mean( abs( results(:,i,4)-gt_d(i) ));
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
    ylim([Meanresults(i,4)-.5,Meanresults(i,4)+.5])
    grid minor
    set(gca,'FontSize',20)
    set(gcf,'color','w');
    if i ==1
%         legend([pt1 pt2 pt3 pt4],{'Measurement','Mean','Groundtruth','2-sigma'},'FontSize',10)
%         set(get(gca,'ylabel'),'rotation',270)
        y = ylabel('8Point','FontSize',30);
%         set(y, 'position', get(y,'position')-[0.01,0.1,0]); 
    end
    
end


Meanresults
MeanError
STDresults



MeanresultsPnP = zeros(secNum,4);
MeanErrorPnP = zeros(secNum,4);
STDresultsPnP = zeros(secNum,4);
% h = figure('units','normalized','outerposition',[0 0 1 0.3]);
for i = 1:secNum
    MeanresultsPnP(i,4) = mean( results_pnp(:,i,4) );
    MeanErrorPnP(i,4) = mean( abs( results_pnp(:,i,4)-gt_d(i) ));
    STDresultsPnP(i,:) = reshape(std(results_pnp(:,i,:)),[1,4]);
    
    checkOutlierPnP=abs(results_pnp(:,i,4)-median( results_pnp(:,i,4) ).*ones(length(results_pnp),1));
    
    [x,y] = find(checkOutlierPnP(:,1)>min(1,max(outlierThes,STDresultsPnP(i,4)*sigmaTimes)));
    results_pnp(x,:,:)=[];
    

    MeanresultsPnP(i,1) = mean( ( results_pnp(:,i,1) )  );
    MeanresultsPnP(i,2) = mean( ( results_pnp(:,i,2)  )  );
    MeanresultsPnP(i,3) = mean( ( results_pnp(:,i,3)  )  );

    MeanErrorPnP(i,1) = mean( abs( results_pnp(:,i,1)-gt(1,i)  )  );
    MeanErrorPnP(i,2) = mean( abs( results_pnp(:,i,2)-gt(2,i)  )  );
    MeanErrorPnP(i,3) = mean( abs( results_pnp(:,i,3)-gt(3,i)  )  );
    
    MeanresultsPnP(i,4) = mean( results_pnp(:,i,4) );
    MeanErrorPnP(i,4) = mean( abs( results_pnp(:,i,4)-gt_d(i) ));
    STDresultsPnP(i,:) = reshape(std(results_pnp(:,i,:)),[1,4]);
    
    nexttile
%     subplot(1,secNum,i)
    pt1 = plot(results_pnp(:,i,4),'b.','Markersize',10);
    hold on
    pt3 = plot(abs(gt_d(i))*ones(1,length(results_pnp)),'r','Linewidth',5);
    hold on
    pt2 = plot(MeanresultsPnP(i,4)*ones(1,length(results_pnp)),'g','Linewidth',5);
    hold on
    pt4 = plot((MeanresultsPnP(i,4)-2*STDresultsPnP(i,4))*ones(1,length(results_pnp)),'k--','Linewidth',3);
    pt4 = plot((MeanresultsPnP(i,4)+2*STDresultsPnP(i,4))*ones(1,length(results_pnp)),'k--','Linewidth',3);
    ylim([MeanresultsPnP(i,4)-.5,MeanresultsPnP(i,4)+.5])
    grid minor
    set(gca,'FontSize',20)
    set(gcf,'color','w');
    if i ==1
%         set(get(gca,'ylabel'),'rotation',0)
        y = ylabel('PnP','FontSize',30);
%         set(y, 'position', get(y,'position')-[0.01,0.1,0]); 
    end
    
    
end

title(t,'rotation Y-axis','FontSize',40)
xlabel(t,'trail number','FontSize',30)
ylabel(t,'Measurement (degree)','FontSize',30)

exportgraphics(h,'/Users/jin/Library/Mobile Documents/com~apple~CloudDocs/iphone_ego/rot_y.png') 


MeanresultsPnP
MeanErrorPnP
STDresultsPnP


length(results)/testNum
length(results_pnp)/testNum