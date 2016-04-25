function draw_camera(A,R,T,txt1,scale,h)



%clear all; close all; load draw_camera; nargin=3; 

if nargin==6
    figure(h); hold on;
else
    figure;
end


%principal point
Cw=-inv(R)*T;
plot3(Cw(1),Cw(2),Cw(3),'.','markersize',20);  
%text(Cw(1),Cw(2),Cw(3),txt1,'fontsize',20);


%camera coordinate system
w=2; col=[0 0 1]; %color and width of the lines
xc=[2 0 0]'; yc=[0 2 0]'; zc=[0 0 2]';
xw=inv(R)*(xc-T);
yw=inv(R)*(yc-T);
zw=inv(R)*(zc-T);

line_(Cw,xw,col,w); %text(xw(1),xw(2),xw(3),'X_c');
line_(Cw,yw,col,w); text(yw(1),yw(2),yw(3),txt1,'fontsize',20); %text(yw(1),yw(2),yw(3),'Y_c');
line_(Cw,zw,col,w); %text(zw(1),zw(2),zw(3),'Z_c');


%field of view;
a=0.7*scale; b=1.5*scale;
P1c=[-a,a,b]';
P2c=[a,a,b]';
P3c=[a,-a,b]';
P4c=[-a,-a,b]';
P1w=inv(R)*(P1c-T);
P2w=inv(R)*(P2c-T);
P3w=inv(R)*(P3c-T);
P4w=inv(R)*(P4c-T);

col=[1 0 0]; w=2;
line_(Cw,P1w,col,w);
line_(Cw,P2w,col,w);
line_(Cw,P3w,col,w);
line_(Cw,P4w,col,w);
line_(P1w,P2w,col,w);
line_(P2w,P3w,col,w);
line_(P3w,P4w,col,w);
line_(P4w,P1w,col,w);

grid on; axis equal;



