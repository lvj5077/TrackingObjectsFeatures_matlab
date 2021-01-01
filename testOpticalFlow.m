clear
close all
clc
addpath('/Users/jin/Q_Mac/mexopencv');
I1 = imread('/Users/jin/Desktop/results/rgb125.png');
I1 = rgb2gray(I1);
I2 = imread('/Users/jin/Desktop/results/rgb530.png');
I2 = rgb2gray(I2);

fx = 486.9024;
fy = 486.9024;
cx = 316.8537;
cy = 219.8356;

figure,
I1_pdic = imread('/Users/jin/Desktop/results/predictions125.jpg');
imshow(I1_pdic)

% % points = detectMinEigenFeatures(I1);
% points = detectORBFeatures(I1);
% [idx,C1] = kmeans(points.Location,500);
% % % prevPts = num2cell(C,[2,6])';
% hold on
% plot(C1(:,1),C1(:,2),'r*')


colorstyle = ["g*","b*","c*","m*","r*","mx","b<","r<"];

obj_dec1 = [1,480,146,480,4;56,280,1,125,3;30,155,500,640,5;56,200,410,525,2;35,125,90,130,1];

objsize = (obj_dec1(:,2)-obj_dec1(:,1)).*(obj_dec1(:,4)-obj_dec1(:,3));
obj_dec1 = [obj_dec1,objsize];

[~,idx] = sort(objsize); 
obj_dec1 = obj_dec1(idx,:);



figure,
I2_pdic = imread('/Users/jin/Desktop/results/predictions530.jpg');
imshow(I2_pdic)
masks_obj1 = zeros(480,640);
prevPts = [];


[w,h] = size(obj_dec1);
for i = 1:w
    mask_tmp = zeros(480,640);
    
    mask_tmp(obj_dec1(i,1):obj_dec1(i,2),obj_dec1(i,3):obj_dec1(i,4)) = 255;
%     subplot(5,2,i*2-1)
%     imshow(mask_tmp)
%     i
%     obj_dec1(i,5)
%     
    masks_obj1(masks_obj1(:,:)==0 & mask_tmp(:,:)>0) = obj_dec1(i,5);
%     figure, imshow(masks_obj)
    mask = zeros(480,640);
    mask(masks_obj1(:,:)==obj_dec1(i,5))=255;
    tempPts = cv.goodFeaturesToTrack(I1, 'MaxCorners', 70, 'QualityLevel', 0.01, 'MinDistance', length(find(mask(:,:)==255))*30/(640*480), 'Mask', mask);
    prevPts = [prevPts,tempPts];
    
    length(tempPts)
    tempMat = cell2mat(tempPts);
    C1 = single(reshape(tempMat',[2,length(tempMat)/2]))';
    hold on
    plot(C1(:,1),C1(:,2),colorstyle( obj_dec1(i,5)+1 ))
end


mask = zeros(480,640);
mask(masks_obj1(:,:)==0)=255;

tempPts = cv.goodFeaturesToTrack(I1, 'MaxCorners', 150, 'QualityLevel', 0.01, 'MinDistance', 0.5*length(find(mask(:,:)==255))*30/(640*480), 'Mask', mask);
length(tempPts)
prevPts = [prevPts,tempPts];
tempMat = cell2mat(tempPts);
C1 = single(reshape(tempMat',[2,length(tempMat)/2]))';
hold on
plot(C1(:,1),C1(:,2),colorstyle(1))
figure,imagesc(masks_obj1)

%%
figure
imshow(I2_pdic);


tempMat = cell2mat(prevPts);
C1 = single(reshape(tempMat',[2,length(tempMat)/2]))';


nextPts = cv.calcOpticalFlowPyrLK(I1, I2, prevPts);
tempMat = cell2mat(nextPts);
C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
% [f2, vpts2] = extractFeatures(I2, C2);

Idpt1 = imread('/Users/jin/Desktop/results/depth125.png');
Idpt1 = imresize(Idpt1,[480,640],'nearest');
Idpt2 = imread('/Users/jin/Desktop/results/depth530.png');
Idpt2 = imresize(Idpt2,[480,640],'nearest');

vpts2 = C2(C2(:,1)>0&C2(:,1)<640&C2(:,2)>0&C2(:,2)<480,:);
vpts1 = C1(C2(:,1)>0&C2(:,1)<640&C2(:,2)>0&C2(:,2)<480,:);


hold on
plot(vpts1(:,1),vpts1(:,2),'g*')

hold on
plot(vpts2(:,1),vpts2(:,2),'r.')
% plot(C2(:,1),C2(:,2),'r.')


figure
imshow(I2_pdic);


obj_dec2 = [1,480,148,500,4;72,400,1,72,3;59,265,445,575,2;64,145,72,115,1;68,125,7,75,6];
objsize = (obj_dec2(:,2)-obj_dec2(:,1)).*(obj_dec2(:,4)-obj_dec2(:,3));
obj_dec2 = [obj_dec2,objsize];

[~,idx] = sort(objsize); 
obj_dec2 = obj_dec2(idx,:);

masks_obj2 = zeros(480,640);
[w,h] = size(obj_dec2);
for i = 1:w
    mask_tmp = zeros(480,640);
    mask_tmp(obj_dec2(i,1):obj_dec2(i,2),obj_dec2(i,3):obj_dec2(i,4)) = 255;
    masks_obj2(masks_obj2(:,:)==0 & mask_tmp(:,:)>0) = obj_dec2(i,5);
end


for i = 1:length(vpts2)
    hold on
    x = round(vpts2(i,1));
    y = round(vpts2(i,2));
%     if (x>0&&x<640&&y>0&&y<480)
     plot(vpts2(i,1),vpts2(i,2),colorstyle( masks_obj2(y,x)+1 ));
%      plot(vpts2(i,1),vpts2(i,2),'ro')
%     end
end

% figure,imagesc(masks_obj2)


obj_stat = zeros(length(vpts1),1); 
for i = 1:length(vpts1)

    x = round(vpts1(i,1));
    y = round(vpts1(i,2));
%     if (x>0&&x<640&&y>0&&y<480)
    colorcode1 = masks_obj1(y,x)+1;
    
    x = round(vpts2(i,1));
    y = round(vpts2(i,2));
	colorcode2 = masks_obj2(y,x)+1;

    if (colorcode2 == colorcode1)
        obj_stat(i) = colorcode1;
    end

    if (colorcode2 < 2 && colorcode1>1)
        obj_stat(i) = colorcode1;
    end
%      plot(vpts2(i,1),vpts2(i,2),'ro')
%     end
end
%%
figure,showMatchedFeatures(I1_pdic,I2_pdic,vpts1(obj_stat>0,:),vpts2(obj_stat>0,:),'montage','PlotOptions',{'yo','yo','y-'});

I1_pdic = imread('/Users/jin/Desktop/results/predictions125.jpg');
I2_pdic = imread('/Users/jin/Desktop/results/predictions530.jpg');
colorstyle = ["g*","b*","c*","m*","r*","mx","b<","r<"];

% obj_stat(obj_stat(:,1)==0) = 8;
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



vpts1 = double(vpts1);
vpts2 = double(vpts2);
pt3d1 = zeros(length(vpts1),4);
pt3d2 = zeros(length(vpts1),4);
for i = 1:length(vpts1)
    x1 = round(vpts1(i,1));
    y1 = round(vpts1(i,2));
    d1 = double(Idpt1(y1,x1))/1000;

    x2 = round(vpts2(i,1));
    y2 = round(vpts2(i,2));
    d2 = double(Idpt2(y2,x2))/1000;

    if (d1>0&&d1<5&&d2>0&&d2<5 )
        pt3d1(i,1) = d1*(double(vpts1(i,2))- cx)/fx;
        pt3d1(i,2) = d1*(double(vpts1(i,1))- cy)/fy;
        pt3d1(i,3) = d1;
        
        pt3d1(i,4) = 1;
        
        pt3d2(i,1) = d2*(double(vpts2(i,2))- cx)/fx;
        pt3d2(i,2) = d2*(double(vpts2(i,1))- cy)/fy;
        pt3d2(i,3) = d2;
        pt3d2(i,4) = 1;
    else
        error = 1
    end
end

obj_stat = obj_stat(pt3d1(:,3)>0,:);
pt3d2 = pt3d2(pt3d1(:,3)>0,:);
vpts2 = vpts2(pt3d1(:,3)>0,:);
vpts1 = vpts1(pt3d1(:,3)>0,:);
pt3d1 = pt3d1(pt3d1(:,3)>0,:);


%%
clc
% close all


[T,inlierPts,inlierObjs] = RANSACpose3d3d_obj(pt3d1,pt3d2,obj_stat);

% T1= eye(4);
% T1(1:3,1:3) = quat2rotm([0.997657 -0.0170876 0.0636572 0.0183562]);
% T1(1:3,4) =[-0.327008 0.109364 -0.234635]';
% T= T1
% inlierObjs'
rotm2eul(T(1:3,1:3))*180/pi
error_Pts = vecnorm(T*pt3d1(inlierPts==1,:)'-pt3d2(inlierPts==1,:)');

mean(error_Pts)

error_Pts = vecnorm(T*pt3d1(:,:)'-pt3d2(:,:)');

for i = 1:7
    if (~isempty(pt3d2(obj_stat(:)==i,:)))
        disp( colorstyle( i) )
        mean(  vecnorm(T*pt3d1(obj_stat(:)==i,:)'-pt3d2(obj_stat(:)==i,:)') )
    end
end


vpts1_inlier = vpts1(inlierPts==1&obj_stat(:)~=5&error_Pts(:)<0.2,:);
vpts2_inlier = vpts2(inlierPts==1&obj_stat(:)~=5&error_Pts(:)<0.2,:);
obj_stat_inler = obj_stat(inlierPts==1&obj_stat(:)~=5&error_Pts(:)<0.2,:);

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


%%
[M, inlierPts] = cv.estimateAffine3D(pt3d1(obj_stat(:)~=5,1:3), pt3d2(obj_stat(:)~=5,1:3));
M
error_Pts = vecnorm(T*pt3d1(inlierPts==1,:)'-pt3d2(inlierPts==1,:)');
mean(error_Pts)
rotm2eul(M(1:3,1:3))*180/pi

vpts1_inlier = vpts1(inlierPts==1,:);
vpts2_inlier = vpts2(inlierPts==1,:);
obj_stat_inler = obj_stat(inlierPts==1,:);

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

%%
clc
idx= find(obj_stat(:)~=5);
essential_matrix = cv.findEssentialMat(vpts1(idx,:), vpts2(idx,:), 'CameraMatrix',[fx 0 cx; 0 fy cy; 0 0 1]);
[R, t, good] = cv.recoverPose(essential_matrix, vpts1(idx,:), vpts2(idx,:),'CameraMatrix',[fx 0 cx; 0 fy cy; 0 0 1]);
good
ep_d = zeros(length(pt3d1),1);
t_x = [0, t(3), t(2);t(3), 0, -t(1);-t(2), t(1), 0];
   
for i = 1:length(pt3d1)
    y1 = pt3d1(i,1:3)./pt3d1(i,3);
    y2 = pt3d2(i,1:3)./pt3d2(i,3);
    ep_d(i) = y2* t_x * R * y1';
    ep_d(i) = y2* t_x * R * y1';
end 

for i = 1:6
    kinderror = ep_d(obj_stat(:)==i);
    i
    disp( colorstyle( i) )
    mean(kinderror)
end