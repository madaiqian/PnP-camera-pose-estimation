% MAIN Illustrates how to use the UPnP algorithm described in:
%
%   A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer
%   Exhaustive Linearization for Robust Camera Pose and Focal Length Estimation.
%   In IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013.
%
% Copyright (C) <2013>  <A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer>
%                       Institut de Robotica i Informatica Industrial (CSIC-UPC)
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.


function [A,point,Rt]=generate_planar_noisy_input_data(n,std_noise,f_value,percentage_outliers,orientation_points,draw_plot)

if nargin<6 draw_plot=''; end

%% True 3d position of the points--------------
xini=-2; xend=2; % bounds 
yini=-2; yend=2;

P12 = orientation_points(2,:) - orientation_points(1,:);
P13 = orientation_points(3,:) - orientation_points(1,:);
plane_coefficients = zeros(4,1);
plane_coefficients(1:3) = cross(P12,P13);

plane_coefficients(4) = sum(plane_coefficients(1:3).*-orientation_points(1,:)'); % ax + by + cz + d = 0


%% Intrinsic camera parameters------------------
fx=f_value; % Focal length in pixels
fy=f_value;
f=0.05; % 50 mm of focal length

m=fx/f;
m=1; f=1;

u0=320; v0=240;
width=640; height=480;
A=[fx/m 0 u0/m 0; 0 fy/m v0/m 0; 0 0 1 0];
std_noise=std_noise*(1/m);

i=1;
num_outliers=ceil(n*percentage_outliers/100);  % num_outliers=4;
while i<=n
    point(i).Xcam(1,1)= random(xini,xend,1);
    point(i).Xcam(2,1)= random(yini,yend,1);
    point(i).Xcam(3,1)= (plane_coefficients(1)*point(i).Xcam(1,1) + plane_coefficients(2)*point(i).Xcam(2,1) + plane_coefficients(4)) / (-plane_coefficients(3)) ; % (a*x+b*y+d)/-c = z
    
    % Add little noise to avoid problems with singularities in methods (NaN's with SVD)
    % We add noise in the Normal's direction
    normals_noise = 0;%0.0001;
    point(i).Xcam(:,1)= point(i).Xcam(:,1) + plane_coefficients(1:3)*random(-normals_noise,normals_noise,1);
    
    % Project points into the image plane
    point(i).Ximg_true=project_3d_2d(A,point(i).Xcam);
    point(i).Ximg_true(3)=f;
    
    % Coordinates in pixels
    if i>n-num_outliers % Generate outliers
        urand=random(1,width,1);
        vrand=random(1,height,1);
        point(i).Ximg_pix_true=[urand,vrand]';
    else
        point(i).Ximg_pix_true=point(i).Ximg_true(1:2)*m;
    end
    
    % Check if the point is in the limits of the image
    uv=point(i).Ximg_pix_true;
    if uv(1)>0 & uv(1)<width & uv(2)>0 & uv(2)<height
        i=i+1;
    end
end


%% Add noise
for i=1:n
    noise=randn(1,2)*std_noise;
    point(i).Ximg(1,1)=point(i).Ximg_true(1)+noise(1);
    point(i).Ximg(2,1)=point(i).Ximg_true(2)+noise(2);
    point(i).Ximg(3,1)=f;
    
    % Noisy coordinates in pixels
    point(i).Ximg_pix=point(i).Ximg(1:2)*m;
end

% The observed data will not be in the camera coordinate system. We put the center of the
% world system on the centroid of the data. We also assume that the data is rotated with
% respect to the camera coordenate system.
centroid=[0 0 0]';
for i=1:n
    centroid=centroid+point(i).Xcam;
end
centroid=centroid/n;

x_rotation=random(0,45,1)*pi/180;
y_rotation=random(0,45,1)*pi/180;
z_rotation=random(0,45,1)*pi/180;
tx=0; ty=0; tz=0;
Rt=return_Rt_matrix(x_rotation,y_rotation,z_rotation,tx,ty,tz);
for i=1:n
    point(i).Xworld=transform_3d(Rt,point(i).Xcam-centroid);
end


%% Plot noisy points in the image plane
if ~strcmp(draw_plot,'donotplot')
    figure; hold on;
    for i=1:n
        plot(point(i).Ximg_pix_true(1),point(i).Ximg_pix_true(2),'.','color',[1 0 0]);
        plot(point(i).Ximg_pix(1),point(i).Ximg_pix(2),'o','color',[0 1 0],'markersize',5);
    end
    title('Noise in image plane','fontsize',20);
    grid on;
    xlim([0,width]);
    ylim([0,height]);

    % draw 3d points
    figure; hold on;
    representation_offset=10;
    plot3(0,0,0,'.','color',[0 0 0],'markersize',20);
    for i=1:n
       plot3(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),'.','color',[1 0 0],'markersize',12);
       txt=sprintf('%d',i);
       text(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),txt);
    end
    grid on;
end


function [x]=random(xini,xend,n)

a=xend-xini;
b=xini;

xu=rand(n,1);
x=xu*a+b;
