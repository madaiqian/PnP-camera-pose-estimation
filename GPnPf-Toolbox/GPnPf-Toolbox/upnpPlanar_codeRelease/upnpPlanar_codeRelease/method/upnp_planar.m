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


function [R,T,Xc,f,best_solution,sol,Beta_opt] = upnp_planar(x3d_h,x2d_h,A)

A(1,1)= 0;
A(2,2)= 0;

Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

% Define control points in a world coordinate system (centered on the 3d points centroid)
Cw = define_control_points_planar2(Xw);

% Compute alphas (linear combination of the control points to represent the 3d points)
Alph = compute_alphas_planar2(Xw,Cw);

% Compute M
M = compute_M_ver3_planar(U,Alph,A);

% Compute kernel M
Kmf = kernel_noise(M,4);

A_UPnP = A;


%% Solve assuming dim(ker(M))=1. with f unknown

Kmf1=Kmf(:,end);

% Control points distance constraint
D = compute_constraint_distance_2param_6eq_2unk_f_unk_planar(Kmf1);
dsq = define_distances_btw_control_points_upnp_planar(Cw);
betas_ = inv(D'*D)*D'*dsq;
beta = sqrt(abs(betas_(1))); % We use the positive value
f = sqrt(abs(betas_(2)))/beta; % f has to be positive, therefore we use the positive value
X1 = beta*Kmf1;

% Xci = beta*Vi1
% Yci = beta*Vi2
% Zci = f*beta*Vi3
X1(3:3:size(Kmf1,1)) = X1(3:3:size(Kmf1,1))*f;

[Cc,Xc,sc] = compute_norm_sign_scaling_factor_upnp_planar(X1,Cw,Alph,Xw);

% Get new A matrix
A_UPnP(1,1) = f;
A_UPnP(2,2) = f;

% Points must be in front of the camera
center = mean(Xc);
if (center(3) < 0)
    Xc = -Xc;
end

[R,T] = getrotT(Xw,Xc);
err(1) = reprojection_error_usingRT(Xw,U,R,T,A_UPnP);

sol(1).Xc = Xc;
sol(1).Cc = Cc;
sol(1).R = R;
sol(1).T = T;
sol(1).F = f;
sol(1).betas = beta;
sol(1).sc = sc;
sol(1).error = err(1);

[min_err,best_solution]=min(err);
Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
T=sol(best_solution).T;


%% Gauss Newton Optimization

Betas=sol(best_solution).betas;
sc=sol(best_solution).sc;
f=sol(best_solution).F;

if best_solution==1
    Betas=[0,0,Betas,f];
elseif best_solution==2
    Betas=[0,Betas,f];
elseif best_solution==3
    Betas=[Betas,f];
end

Kmf1=Kmf(:,end-2);
Kmf2=Kmf(:,end-1);
Kmf3=Kmf(:,end);
Kernel=[Kmf1,Kmf2,Kmf3];

Beta_opt = Betas(1:3)/sc;




%% A&D implementation of Horn's RT
function [R, T] = getrotT(wpts,cpts)

Xw=wpts;
Xc=cpts;

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

Xc_app=(R*Xw'+repmat(T,[1,n]))';

tmp=Xc_app-Xc;
tmp2=tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2;
error=sqrt(sum(tmp2)/n);


%% Euclidean Reprojection Error Calculation Function
function [err,Urep] = reprojection_error_usingRT(Xw,U,R,T,A)

n=size(Xw,1);

P=A*[R,T];
Xw_h=[Xw,ones(n,1)];
Urep_=(P*Xw_h')';

Urep=zeros(n,2);
Urep(:,1)=Urep_(:,1)./Urep_(:,3);
Urep(:,2)=Urep_(:,2)./Urep_(:,3);

err_=sqrt((U(:,1)-Urep(:,1)).^2+(U(:,2)-Urep(:,2)).^2);
err=sum(err_)/n;


