clear
close all
clc

I1 = imread('/Users/jin/Desktop/results/rgb125.png');
I1 = rgb2gray(I1);
I2 = imread('/Users/jin/Desktop/results/rgb530.png');
I2 = rgb2gray(I2);


% points = detectMinEigenFeatures(I1);
% points = detectORBFeatures(I1);
% [idx,C1] = kmeans(points.Location,500);
% % prevPts = num2cell(C,[2,6])';

figure,
imshow(I1)
% hold on
% plot(C1(:,1),C1(:,2),'r*')


prevPts = cv.goodFeaturesToTrack(I1, 'MaxCorners', 500, 'QualityLevel', 0.01, 'MinDistance', 20);

tempMat = cell2mat(prevPts);
C1 = single(reshape(tempMat',[2,length(tempMat)/2]))';
hold on
plot(C1(:,1),C1(:,2),'g*')



%%
% tempMat = cell2mat(prevPts);
% C = single(reshape(tempMat,[length(tempMat)/2,2]));
[f1, vpts1] = extractFeatures(I1, C1);
% vpts1 = C1;
figure
ax = axes;
imshow(imread('/Users/jin/Desktop/results/predictions125.jpg'))

% 1 bottle, 2 blue chair, 3 gray chair, 4 person, 0 background
obj_vpts1 = zeros(length(vpts1),1);

obj_vpts1(vpts1(:,1)>146 & vpts1(:,1)<480 &vpts1(:,2)>0 & vpts1(:,2)<420) = 4;
obj_vpts1(vpts1(:,1)>0 & vpts1(:,1)<125 &vpts1(:,2)>56 & vpts1(:,2)<280) = 3;
obj_vpts1(vpts1(:,1)>500 & vpts1(:,1)<640 &vpts1(:,2)>30 & vpts1(:,2)<155) = 6;
obj_vpts1(vpts1(:,1)>410 & vpts1(:,1)<525 &vpts1(:,2)>56 & vpts1(:,2)<200) = 2;
obj_vpts1(vpts1(:,1)>90 & vpts1(:,1)<130 &vpts1(:,2)>35 & vpts1(:,2)<115) = 1;

hold on
plot(vpts1(obj_vpts1==1,1),vpts1(obj_vpts1==1,2),'r*')
hold on
plot(vpts1(obj_vpts1==2,1),vpts1(obj_vpts1==2,2),'m*')
hold on
plot(vpts1(obj_vpts1==3,1),vpts1(obj_vpts1==3,2),'b*')
hold on
plot(vpts1(obj_vpts1==4,1),vpts1(obj_vpts1==4,2),'c*')
hold on
plot(vpts1(obj_vpts1==6,1),vpts1(obj_vpts1==6,2),'k*')
hold on
plot(vpts1(obj_vpts1==0,1),vpts1(obj_vpts1==0,2),'g*')
%%

figure
imshow(imread('/Users/jin/Desktop/results/predictions530.jpg'));
hold on
nextPts = cv.goodFeaturesToTrack(I2, 'MaxCorners', 500, 'QualityLevel', 0.01, 'MinDistance', 20);

tempMat = cell2mat(nextPts);
C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';
% hold on
% plot(C2(:,1),C2(:,2),'g*')

[f2, vpts2] = extractFeatures(I2, C2);

%%
obj_vpts2 = zeros(length(vpts2),1);

obj_vpts2(vpts2(:,1)>148 & vpts2(:,1)<500 &vpts2(:,2)>0 & vpts2(:,2)<480) = 4;
obj_vpts2(vpts2(:,1)>0 & vpts2(:,1)<72 &vpts2(:,2)>72 & vpts2(:,2)<400) = 3;
obj_vpts2(vpts2(:,1)>445 & vpts2(:,1)<575 &vpts2(:,2)>59 & vpts2(:,2)<265) = 2;
obj_vpts2(vpts2(:,1)>72 & vpts2(:,1)<115 &vpts2(:,2)>64 & vpts2(:,2)<145) = 1;
obj_vpts2(vpts2(:,1)>7 & vpts2(:,1)<75 &vpts2(:,2)>68 & vpts2(:,2)<125) = 5;
hold on
plot(vpts2(obj_vpts2==1,1),vpts2(obj_vpts2==1,2),'r*')
hold on
plot(vpts2(obj_vpts2==2,1),vpts2(obj_vpts2==2,2),'m*')
hold on
plot(vpts2(obj_vpts2==3,1),vpts2(obj_vpts2==3,2),'b*')
hold on
plot(vpts2(obj_vpts2==4,1),vpts2(obj_vpts2==4,2),'c*')
hold on
plot(vpts2(obj_vpts2==5,1),vpts2(obj_vpts2==5,2),'k*')
hold on
plot(vpts2(obj_vpts2==0,1),vpts2(obj_vpts2==0,2),'g*')
%%
figure, 

for i= 1:5
    objid = i-1;
    subplot(5,1,i)
    testf1 = f1(obj_vpts1==objid,:);
    testf2 = f2(obj_vpts2==objid,:);

    testvpts1 = vpts1(obj_vpts1==objid,:); 
    testvpts2 = vpts2(obj_vpts2==objid,:); 
    
    indexPairs = matchFeatures(testf1, testf2,'MatchThreshold',100,'MaxRatio',0.9) ;
    matchedPoints1 = testvpts1(indexPairs(:, 1),:);
    matchedPoints2 = testvpts2(indexPairs(:, 2),:);
    showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage');

end

%%
figure,
testf1 = f1;
testf2 = f2;

testvpts1 = vpts1; 
testvpts2 = vpts2; 

indexPairs = matchFeatures(testf1, testf2,'MatchThreshold',100,'MaxRatio',0.9) ;
matchedPoints1 = testvpts1(indexPairs(:, 1),:);
matchedPoints2 = testvpts2(indexPairs(:, 2),:);
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage');
