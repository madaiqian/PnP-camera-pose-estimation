function [A,point,Rt]=generate_noisy_input_data(n,std_noise,f_value,percentage_outliers,draw_plot)

if nargin<5 draw_plot=''; end

%% Bounds where points will be generated
xini=-2; xend=2; 
yini=-2; yend=2;
zini=4; zend=9;

%% Intrinsic camera parameters
fx=f_value;     % Focal length in pixels
fy=f_value;
fz=1;           % distance between optical center and image plane
u0=320; v0=240;
width=640; height=480;
A=[fx 0 u0 0; 0 fy v0 0; 0 0 1 0];

%% Generate points
i=1;
num_outliers=ceil(n*percentage_outliers/100);  
while i<=n
    point(i).Xcam(1,1)=random(xini,xend,1);
    point(i).Xcam(2,1)=random(yini,yend,1);
    point(i).Xcam(3,1)=random(zini,zend,1);
    
    % Project points into the image plane
    point(i).Ximg_true=project_3d_2d(A,point(i).Xcam);
    point(i).Ximg_true(3)=fz;
    
    % Coordinates in pixels
    if i>n-num_outliers % Generate outliers
        urand=random(1,width,1);
        vrand=random(1,height,1);
        point(i).Ximg_pix_true=[urand,vrand]';
    else
        point(i).Ximg_pix_true=point(i).Ximg_true(1:2);
    end
    
    % Check if the point is into the limits of the image
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
    point(i).Ximg(3,1)=fz;
    
    % Noisy coordinates in pixels
    point(i).Ximg_pix=point(i).Ximg(1:2);
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


% Plot noisy points in the image plane
if ~strcmp(draw_plot,'donotplot')
    figure; hold on;
    for i=1:n
        plot(point(i).Ximg_pix_true(1),point(i).Ximg_pix_true(2),'.','color',[1 0 0]);
        plot(point(i).Ximg_pix(1),point(i).Ximg_pix(2),'o','color',[0 0 1],'markersize',5);
    end
    title('Noise in image plane','fontsize',20);
    grid on;
    xlim([0,width]);
    ylim([0,height]);
    legend('True 2D','Observation');
end



function [x]=random(xini,xend,n)

a=xend-xini;
b=xini;

xu=rand(n,1);
x=xu*a+b;
