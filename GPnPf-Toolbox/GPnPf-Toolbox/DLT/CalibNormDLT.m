% CalibNormDLT - computes camera projection matrix from 3d scene points
% and corresponding 2d image points with normalized DLT(direct linear transformation).
%
% Usage:
%           P = CalibNormDLT(x, X)
%
% Input:
%           x : 2xn, nx2, 3xn or nx3 image points
%           X : 3xn, nx3, 4xn or nx4 object points
%
% Output:
%           P : 3x4 camera projection matrix
%
% cf.:
%           x = P*X
%
% This code follows the algorithm given by
% [1] R. Hartley and A. Zisserman "Multiple View Geometry in Computer Vision,"
%     pp.178-181, 2003.
% [2] E. Trucco and A. Verri, "Introductory Techniques for 3-D
%     Computer Vision," pp. 132-134, 1998.
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Jun. 2008  - Original version.
% Jun. 2011  - Input is modified.


function P = CalibNormDLT(x, X)

%% The Number of Points
xSize  = size(x);
XSize  = size(X);
noPnts = length(x);

%% image points
if (xSize(2) == 2)
    x = x';
elseif (xSize(1) == 3)
    x = x(1:2,:);
elseif (xSize(2) == 3)
    x = x(:,1:2)';
end

%% object points
if (XSize(2) == 3)
    X = X';
elseif (XSize(1) == 4)
    X = X(1:3,:);
elseif (XSize(2) == 4)
    X = X(:,1:3)';
end


%% centroids of the points
centroid1 = mean(x(1:2,:)')';
centroid2 = mean(X(1:3,:)')';

%% Shift the origin of the points to the centroid
x(1,:) = x(1,:) - centroid1(1);
x(2,:) = x(2,:) - centroid1(2);

X(1,:) = X(1,:) - centroid2(1);
X(2,:) = X(2,:) - centroid2(2);
X(3,:) = X(3,:) - centroid2(3);

%% Normalize the points so that the average distance from the origin is equal to sqrt(2).
averagedist1 = mean(sqrt(x(1,:).^2 + x(2,:).^2));
averagedist2 = mean(sqrt(X(1,:).^2 + X(2,:).^2 + X(3,:).^2));

scale1 = sqrt(2)/averagedist1;
scale2 = sqrt(3)/averagedist2;

x(1:2,:) = scale1*x(1:2,:);
X(1:3,:) = scale2*X(1:3,:);

%% similarity transform 1
T1 = [scale1         0  -scale1*centroid1(1)
           0    scale1  -scale1*centroid1(2)
           0         0                     1];

%% similarity transform 2
T2 = [scale2       0      0   -scale2*centroid2(1)
           0  scale2      0   -scale2*centroid2(2)
           0       0  scale2  -scale2*centroid2(3)
           0       0       0                     1];

%% Compute the Camera Projection Matrix
A = [X(1,:)'            X(2,:)'             X(3,:)'             ones(noPnts,1) ...
     zeros(noPnts,1)    zeros(noPnts,1)     zeros(noPnts,1)     zeros(noPnts,1) ...
     -x(1,:)'.*X(1,:)'  -x(1,:)'.*X(2,:)'   -x(1,:)'.*X(3,:)'   -x(1,:)';
     zeros(noPnts,1)    zeros(noPnts,1)     zeros(noPnts,1)     zeros(noPnts,1) ...
     X(1,:)'            X(2,:)'             X(3,:)'             ones(noPnts,1) ...
     -x(2,:)'.*X(1,:)'  -x(2,:)'.*X(2,:)'   -x(2,:)'.*X(3,:)'   -x(2,:)'];
     
[U,D,V] = svd(A);
P = reshape(V(:,12),4,3)';

%% Denormalization
P = inv(T1)*P*T2;