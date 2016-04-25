function [R, T]=getrotT(wpts,cpts)

Xw=wpts;
Xc=cpts;

% This routine solves the exterior orientation problem for a point cloud
%  given in both camera and world coordinates. 
  
% wpts = 3D points in arbitrary reference frame
% cpts = 3D points in camera reference frame
  
n=size(wpts,1);
M=zeros(3);

ccent=mean(cpts);
wcent=mean(wpts);

for i=1:3
  cpts(:,i)=cpts(:,i)-ccent(i)*ones(n,1);
  wpts(:,i)=wpts(:,i)-wcent(i)*ones(n,1);
end

for i=1:n
   M=M+cpts(i,:)'*wpts(i,:);
end

[U S V]=svd(M);

detu_detv=det(U)*det(V);

if detu_detv>0
    E=eye(3);
elseif detu_detv<0
    E=eye(3);
    E(3,3)=-1;
end
R=U*E*V';

T=ccent'-R*wcent';

%check approximation
Xc_app=(R*Xw'+repmat(T,[1,n]))';

%error sol2
tmp=Xc_app-Xc;
tmp2=tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2;
error=sqrt(sum(tmp2)/n);






