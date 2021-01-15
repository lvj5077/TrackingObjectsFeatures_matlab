s=1000;
fx = 605.0697021484375 ;
fy = 605.45703125;
skew = 0;
cx = 325.93792724609375;
cy = 247.49156188964844;
c_intrinsic = cameraIntrinsics([605.0697021484375,605.45703125],[325.93792724609375,247.49156188964844],[480,640],'RadialDistortion',[0.14949573576450348, -0.49464672803878784, 0.42560458183288574],'TangentialDistortion',[0.00035660676076076925, 0.0001328527432633564]);

c_intrinsic = cameraIntrinsics([486.9024,486.9024],[316.8537,219.8356],[480,640]);
fx = 486.9024;
fy = 486.9024;
cx = 316.8537;
cy = 219.8356;
K = [fx 0 cx;0 fy cy;0 0 1];
I2 = imread('/Users/jin/Desktop/moveL515/color/1608253889.200758241.png');

I1 = imread('/Users/jin/Desktop/moveL515/color/1608253894.898903038.png');


I1 = imread('/Users/jin/Desktop/results/rgb125.png');
I2 = imread('/Users/jin/Desktop/results/rgb530.png');
% I1 = I2;
figure
imshowpair(I1, I2, 'montage'); 
title('Original Images');
%%
% Detect feature points
mask = zeros(480,640);
mask(1:480,146:480)=255;
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
%%
matchedPoints1 = vpts1(obj_stat(:,1)~=11,:);
matchedPoints2 = vpts2(obj_stat(:,1)~=11,:);
% Visualize correspondences
figure
showMatchedFeatures(I1, I1, matchedPoints1, matchedPoints1);
title('Tracked Features');


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
%%
% Detect dense feature points
imagePoints1 = detectORBFeatures(rgb2gray(I1));
%%
% Create the point tracker
tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% Initialize the point tracker
imagePoints1 = imagePoints1.Location;
initialize(tracker, imagePoints1, I1);

% Track the points
[imagePoints2, validIdx] = step(tracker, I2);
matchedPoints1 = imagePoints1(validIdx, :);
matchedPoints2 = imagePoints2(validIdx, :);
%%
% Compute the camera matrices for each position of the camera
% The first camera is at the origin looking along the X-axis. Thus, its
% rotation matrix is identity, and its translation vector is 0.
camMatrix1 = cameraMatrix(cameraParams, eye(3), [0 0 0]);
camMatrix2 = cameraMatrix(cameraParams, R', -t*R');

% Compute the 3-D points
points3D = triangulate(matchedPoints1, matchedPoints2, camMatrix1, camMatrix2);

% Get the color of each reconstructed point
numPixels = size(I1, 1) * size(I1, 2);
allColors = reshape(I1, [numPixels, 3]);
colorIdx = sub2ind([size(I1, 1), size(I1, 2)], round(matchedPoints1(:,2)), ...
    round(matchedPoints1(:, 1)));
color = allColors(colorIdx, :);

color = color((points3D(:,3))<300 & (points3D(:,3))>0,:);
points3D = points3D((points3D(:,3))<300 & (points3D(:,3))>0,:);
% Create the point cloud
ptCloud = pointCloud(points3D, 'Color', color);

% Visualize the camera locations and orientations
cameraSize = 0.3;
figure
plotCamera('Size', cameraSize, 'Color', 'r', 'Label', '1', 'Opacity', 0);
hold on
grid on
plotCamera('Location', t, 'Orientation', R, 'Size', cameraSize, ...
    'Color', 'b', 'Label', '2', 'Opacity', 0);

% Visualize the point cloud
pcshow(ptCloud, 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 45);

% Rotate and zoom the plot
camorbit(0, -30);
camzoom(1.5);

% Label the axes
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis')

title('Up to Scale Reconstruction of the Scene');