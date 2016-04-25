function Cw=define_control_points_planar2(Xw)

%clear all; close all; load define_control_points_planar2;

m=mean(Xw);
XtX=Xw'*Xw;
[V,S]=eig(XtX);
v1=V(:,3); v1=v1/norm(v1);
v2=V(:,2); v2=v2/norm(v2);

% [u,s,v]=svd(Xw);
% v1_=v(:,1); v1_=v1_/norm(v1_);
% v2_=v(:,2); v2_=v2_/norm(v2_);

Cw=zeros(3,3);
Cw(1,:)=m;
Cw(2,:)=m+v1';
Cw(3,:)=m+v2';


