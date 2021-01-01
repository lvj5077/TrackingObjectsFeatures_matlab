s=1000;
fx = 605.0697021484375 ;
fy = 605.45703125;
skew = 0;
cx = 325.93792724609375;
cy = 247.49156188964844;
c_intrinsic = cameraIntrinsics([605.0697021484375,605.45703125],[325.93792724609375,247.49156188964844],[480,640],'RadialDistortion',[0.14949573576450348, -0.49464672803878784, 0.42560458183288574],'TangentialDistortion',[0.00035660676076076925, 0.0001328527432633564]);

K = [fx 0 cx;0 fy cy;0 0 1];
I2 = imread('/Users/jin/Desktop/moveL515/color/moving/1608253888.528166720.png');
I1 = imread('/Users/jin/Desktop/moveL515/color/moving/1608253890.176004345.png');
figure
imshowpair(I1, I2, 'montage'); 
title('Original Images');

% Detect feature points
imagePoints1 = detectORBFeatures(rgb2gray(I1));

% Visualize detected points
figure
imshow(I1, 'InitialMagnification', 50);
title('150 Strongest Corners from the First Image');
hold on
plot(selectStrongest(imagePoints1, 250));

% Create the point tracker
tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% Initialize the point tracker
imagePoints1 = imagePoints1.Location;
initialize(tracker, imagePoints1, I1);

% Track the points
[imagePoints2, validIdx] = step(tracker, I2);
matchedPoints1 = imagePoints1(validIdx, :);
matchedPoints2 = imagePoints2(validIdx, :);

% Visualize correspondences
figure
showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
title('Tracked Features');

% Estimate the fundamental matrix
[fMatrix, epipolarInliers] = estimateFundamentalMatrix(...
  matchedPoints1, matchedPoints2, 'Method', 'MSAC', 'NumTrials', 100000,'DistanceType','Sampson','InlierPercentage',60); %'Sampson' (default) | 'Algebraic'

% Find epipolar inliers
inlierPoints1 = matchedPoints1(epipolarInliers, :);
inlierPoints2 = matchedPoints2(epipolarInliers, :);

% Display inlier matches
figure
showMatchedFeatures(I1, I2, inlierPoints1, inlierPoints2);
title('Epipolar Inliers');

[R, t] = relativeCameraPose(fMatrix, c_intrinsic, inlierPoints1, inlierPoints2)
euler = rotm2eul(R)
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