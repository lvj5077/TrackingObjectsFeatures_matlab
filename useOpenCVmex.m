s=1000;
fx = 605.0697021484375 ;
fy = 605.45703125;
skew = 0;
cx = 325.93792724609375;
cy = 247.49156188964844;
c_intrinsic = cameraIntrinsics([605.0697021484375,605.45703125],[325.93792724609375,247.49156188964844],[480,640],'RadialDistortion',[0.14949573576450348, -0.49464672803878784, 0.42560458183288574],'TangentialDistortion',[0.00035660676076076925, 0.0001328527432633564]);

% c_intrinsic = cameraIntrinsics([486.9024,486.9024],[316.8537,219.8356],[480,640]);
% fx = 486.9024;
% fy = 486.9024;
% cx = 316.8537;
% cy = 219.8356;
% K = [fx 0 cx;0 fy cy;0 0 1];
I2 = imread('/Users/jin/Desktop/moveL515/color/1608253704.216312740.png');

I1 = imread('/Users/jin/Desktop/moveL515/color/1608253750.800028766.png');
% I1 = imread('/Users/jin/Desktop/moveL515/color/11608253762.937677485.png');
% I1 = I2;
figure
imshowpair(I1, I2, 'montage'); 
title('Original Images');
%%
% Detect feature points
mask = zeros(480,640);
% mask(1:480,146:480)=255;
mask = 255*ones(480,640)-mask;
tempPts = cv.goodFeaturesToTrack(rgb2gray(I1), 'MaxCorners', 500, 'QualityLevel', 0.005, 'MinDistance',15, 'Mask',mask);
prevPts = tempPts;
% Visualize detected points
% figure
% imshow(I1);
% title('Corners from the First Image');

tempMat = cell2mat(tempPts);
C1 = single(reshape(tempMat',[2,length(tempMat)/2]))';
imagePoints1 = C1;
% hold on
% plot(C1(:,1),C1(:,2),'ro')

nextPts = cv.calcOpticalFlowPyrLK(rgb2gray(I1), rgb2gray(I2), prevPts);
tempMat = cell2mat(nextPts);
C2 = single(reshape(tempMat',[2,length(tempMat)/2]))';

matchedPoints2 = C2(C2(:,1)>0&C2(:,1)<640&C2(:,2)>0&C2(:,2)<480,:);
matchedPoints1 = C1(C2(:,1)>0&C2(:,1)<640&C2(:,2)>0&C2(:,2)<480,:);

% Visualize correspondences
% figure
% showMatchedFeatures(I1, I1, matchedPoints1, matchedPoints1);
% title('Tracked Features');


% Estimate the fundamental matrix
[fMatrix, epipolarInliers] = estimateFundamentalMatrix(...
  matchedPoints1, matchedPoints2, 'Method', 'MSAC', 'NumTrials', 100000,'DistanceType','Sampson','InlierPercentage',50); %'Sampson' (default) | 'Algebraic'

% Find epipolar inliers
inlierPoints1 = matchedPoints1(epipolarInliers, :);
inlierPoints2 = matchedPoints2(epipolarInliers, :);

% Display inlier matches
figure
showMatchedFeatures(I1, I2, inlierPoints1, inlierPoints2);


title('Epipolar Inliers with object');
[R, t] = relativeCameraPose(fMatrix, c_intrinsic, inlierPoints1, inlierPoints2);
euler = rotm2eul(R)


epiLines = epipolarLine(fMatrix',matchedPoints2(epipolarInliers,:));
points = lineToBorderPoints(epiLines,size(I1));
line(points(:,[1,3])',points(:,[2,4])');