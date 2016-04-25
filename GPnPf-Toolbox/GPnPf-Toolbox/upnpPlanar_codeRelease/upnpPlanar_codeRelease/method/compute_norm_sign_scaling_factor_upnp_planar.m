function [Cc,Xc,sc] = compute_norm_sign_scaling_factor_upnp_planar(X1,Cw,Alph,Xw)
 
%clear all; close all; load compute_norm_sign_scaling_factor;


%n=size(Xw,1); %number of data points

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(3,3);
for i=1:3
    Cc_(i,:)=X1(3*i-2:3*i);
end

% %position of reference points in camera coordinates
% Xc_=Alph*Cc_;
% 
% %compute distances in world coordinates w.r.t. the centroid
% centr_w=mean(Xw);
% centroid_w=repmat(centr_w,[n,1]);
% tmp1=Xw-centroid_w;
% dist_w=sqrt(sum(tmp1.^2,2));
% 
% %compute distances in camera coordinates w.r.t. the centroid
% centr_c=mean(Xc_);
% centroid_c=repmat(centr_c,[n,1]);
% tmp2=Xc_-centroid_c;
% dist_c=sqrt(sum(tmp2.^2,2));
% 
% %scale_factors=dist_c./dist_w;
% %sc=median(scale_factors); %sc=0.065;
%  
% %least squares solution for the scale factor
% sc=1/(inv(dist_c'*dist_c)*dist_c'*dist_w);

% t=1;
% for i=1:n
%     for j=i+1:n
%         pw=Xw(i,:);
%         qw=Xw(j,:);
%         pc=Xc_(i,:);
%         qc=Xc_(j,:);
%         
%         Dw(t)=sqrt((pw(1)-qw(1))^2+(pw(2)-qw(2))^2+(pw(3)-qw(3))^2);
%         Dc(t)=sqrt((pc(1)-qc(1))^2+(pc(2)-qc(2))^2+(pc(3)-qc(3))^2);
%       
%         t=t+1;
%     end
% end
% scale_factors1=(Dc./Dw)';
% sc1=mean(scale_factors1);
% 
% Cw=[1 0 0; 0 1 0; 0 0 1; 0 0 0];
 t=1;
for i=1:3
    for j=i+1:3
        cw=Cw(i,:);
        dw=Cw(j,:);
        cc=Cc_(i,:);
        dc=Cc_(j,:);
        
        Dcw(t,1)=sqrt((cw(1)-dw(1))^2+(cw(2)-dw(2))^2+(cw(3)-dw(3))^2);
        Dcc(t,1)=sqrt((cc(1)-dc(1))^2+(cc(2)-dc(2))^2+(cc(3)-dc(3))^2);
        t=t+1;
    end
end
sc=1/(inv(Dcc'*Dcc)*Dcc'*Dcw);

% 
% scale_factors2=(Dcc./Dcw)';
% sc2=mean(scale_factors2);


%compute scale factor using least squares approach -- equals to the mean --
% scale_factors3=dist_c./dist_w;
% sc3=mean(scale_factors3); %sc=0.065;
% sc=sc3; 


%sc4
% num=dist_c.*dist_w;
% den=dist_c.*dist_c;
% sc4=sum(den)/sum(num);
% sc=sc4; sc=0.577;


%scale position of the control points
Cc=Cc_/sc;

%rescaled position of the reference points
Xc=Alph*Cc;

%change the sign if necessary. z negative is no possible in camera
%coordinates
neg_z=find(Xc(:,3)<0);
if size(neg_z,1)>=1
    Xc=Xc*(-1);
end


        
        
