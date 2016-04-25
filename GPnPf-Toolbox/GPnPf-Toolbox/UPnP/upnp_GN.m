% upnp_GN UPnP+Gauss Newton algorithm described in:
%
%   A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer
%   Exhaustive Linearization for Robust Camera Pose and Focal Length Estimation.
%   In IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013.
%
% Copyright (C) <2013>  <A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer>
%                       Institut de Robotica i Inform�tica Industrial (CSIC-UPC)
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
%

function [R,T,Xc,f,best_solution,sol,Beta_opt] = upnp_GN(x3d_h,x2d_h,A,maximum_iterations,N3)

A(1,1)= 0;
A(2,2)= 0;

Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

THRESHOLD_REPROJECTION_ERROR3=20;

% Define control points in a world coordinate system (centered on the 3d points centroid)
Cw=define_control_points();

% Compute alphas (linear combination of the control points to represent the 3d points)
Alph=compute_alphas(Xw,Cw);


%Compute M for fu=fv both unknown
Mf=compute_M_ver3(U,Alph,A);
%Compute kernel M for fu=fv both unknown
Kmf=kernel_noise(Mf,4);

Kmf1=Kmf(:,end-2);
Kmf2=Kmf(:,end-1);
Kmf3=Kmf(:,end);
Kernel=[Kmf1,Kmf2,Kmf3];

A_UPnP = A;

%% 1.-Solve assuming dim(ker(M))=1. with f unknown;------------------------------
Kmf1=Kmf(:,end);

%control points distance constraint
D=compute_constraint_distance_2param_6eq_2unk_f_unk(Kmf1);
dsq=define_distances_btw_control_points();
betas_=(D'*D)\D'*dsq;
beta=sqrt(abs(betas_(1))); % We use the positive value
f=sqrt(abs(betas_(2)))/beta; % f has to be positive, therefore we use the positive value
X1=beta*Kmf1;

% Xci = beta*Vi1
% Yci = beta*Vi2
% Zci = f*beta*Vi3
X1(3:3:size(Kmf1,1))=X1(3:3:size(Kmf1,1))*f;

[Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw);

% Get new A matrix
A_UPnP(1,1) = f;
A_UPnP(2,2) = f;

[R,T]=getrotT(Xw,Xc);  % Solve exterior orientation
err(1)=reprojection_error_usingRT(Xw,U,R,T,A_UPnP);

sol(1).Xc=Xc;
sol(1).Cc=Cc;
sol(1).R=R;
sol(1).T=T;
sol(1).f=f;
sol(1).betas=beta;
sol(1).sc=sc;
sol(1).error=err(1);

% ---- Refine the solution iterating over the betas (N=1) ----
Betas=[0,0,sol(1).betas,f];
Beta0=Betas;
Beta0(1:3)=Betas(1:3)/sc;
[Xc_opt,R_opt,T_opt,err_opt,Beta_opt] = UPnP_optimize_betas_gauss_newton(Kernel,Cw,Beta0,Alph,Xw,U,A);

% Just update R,T,Xc if Gauss Newton improves results (which is almost always)
if (err_opt<err(1)) && (Beta_opt(4) > 0) && (Beta_opt(4) < 1e6)
    sol(1).R=R_opt;
    sol(1).T=T_opt;
    sol(1).Xc=Xc_opt;
    sol(1).betas=Beta_opt(3);
    sol(1).f=Beta_opt(4);
    sol(1).error=err_opt;
    err(1)=err_opt;
end


%% 2.-Solve assuming dim(ker(M))=2. with f unknown
Kmf1=Kmf(:,end-1);
Kmf2=Kmf(:,end);

% Control points distance constraint
D=compute_constraint_distance_3param_6eq_6unk_f_unk(Kmf1,Kmf2);
dsq=define_distances_btw_control_points();
betas_=D\dsq;

beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2));

solutions = generate_all_possible_solutions_for_f_unk(betas_, 2);
num_solutions= size(solutions,1);
for i=1:num_solutions
    aux_beta1= solutions(i,1);
    aux_beta2= solutions(i,2);
    aux_f= solutions(i,3);
    
    % Xci = beta*Vi1
    % Yci = beta*Vi2
    % Zci = f*beta*Vi3
    X2=aux_beta1*Kmf1+aux_beta2*Kmf2;
    X2_aux=X2;
    X2_aux(3:3:size(Kmf1,1))=X2(3:3:size(Kmf1,1))*aux_f;
    
    [Cc,Xc]=compute_norm_sign_scaling_factor(X2_aux,Cw,Alph,Xw);
    
    % get new A matrix
    A_UPnP(1,1) = aux_f;
    A_UPnP(2,2) = aux_f;
    
    [R,T]=getrotT(Xw,Xc);  %solve exterior orientation
    error_aux=reprojection_error_usingRT(Xw,U,R,T,A_UPnP);
    
    if ((i==1) || (error_aux < smallest_error))
        beta1= aux_beta1;
        beta2= aux_beta2;
        f= aux_f;
        smallest_error=error_aux;
    end
end

X2=beta1*Kmf1+beta2*Kmf2;
X2(3:3:size(Kmf1,1))=X2(3:3:size(Kmf1,1))*f;
[Cc,Xc,sc]=compute_norm_sign_scaling_factor(X2,Cw,Alph,Xw);

% Get new A matrix
A_UPnP(1,1) = f;
A_UPnP(2,2) = f;

[R,T]=getrotT(Xw,Xc);  % Solve exterior orientation
err(2)=reprojection_error_usingRT(Xw,U,R,T,A_UPnP);

sol(2).Xc=Xc;
sol(2).Cc=Cc;
sol(2).R=R;
sol(2).T=T;
sol(2).f=f;
sol(2).betas=[beta1,beta2];
sol(2).sc=sc;
sol(2).error=err(2);

% ---- Refine the solution iterating over the betas (N=2) ----
Betas=[0,sol(2).betas,f];
Beta0=Betas;
Beta0(1:3)=Betas(1:3)/sc;
[Xc_opt,R_opt,T_opt,err_opt,Beta_opt] = UPnP_optimize_betas_gauss_newton(Kernel,Cw,Beta0,Alph,Xw,U,A);

% Just update R,T,Xc if Gauss Newton improves results (which is almost always)
if (err_opt<err(2)) && (Beta_opt(4) > 0) && (Beta_opt(4) < 1e6)
    sol(2).R=R_opt;
    sol(2).T=T_opt;
    sol(2).Xc=Xc_opt;
    sol(2).betas=Beta_opt(2:3);
    sol(2).f=Beta_opt(4);
    sol(2).error=err_opt;
    err(2)=err_opt;
end


%% 3.-Solve assuming dim(ker(M))=3. with f unknown;------------------------------------------
if ((N3) && (min(err)>THRESHOLD_REPROJECTION_ERROR3)) % Just compute if we do not have a good solution in the previous cases
    
    Kmf1=Kmf(:,end-2);
    Kmf2=Kmf(:,end-1);
    Kmf3=Kmf(:,end);
    
    D=compute_constraint_distance_orthog_4param_9eq_12unk_f_unk(Kmf1,Kmf2,Kmf3);
    dsq=define_distances_btw_control_points();
    
    lastcolumn=[-dsq',0,0,0]';
    D_=[D,lastcolumn];
    Kd=null(D_);
    
    P=compute_permutation_constraint3_f_unk(Kd);
    lambdas_=kernel_noise(P,1);
    
    % -----------   EFFICIENT EXPLORATION OF MINIMAL EQUATION SETS  ----------------------------
    load_order_data= 1;
    n_samples= 100;
    ordered_equations= generate_ordered_list_of_equations(n_samples, load_order_data);
    
    [beta1,beta2,beta3,f] = equation_exploration(ordered_equations,Cw,Alph,Xw,U,A_UPnP,lambdas_,Kd,Kmf1,Kmf2,Kmf3,maximum_iterations);
    
    
    % Xci = beta*Vi1
    % Yci = beta*Vi2
    % Zci = f*beta*Vi3
    X3=beta1*Kmf1+beta2*Kmf2+beta3*Kmf3;
    X3_aux=X3;
    X3_aux(3:3:size(Kmf1,1))=X3(3:3:size(Kmf1,1))*f;
    
    [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X3_aux,Cw,Alph,Xw);
    
    % get new A matrix
    A_UPnP(1,1) = f;
    A_UPnP(2,2) = f;
    
    [R,T]=getrotT(Xw,Xc);  % Solve exterior orientation
    err(3)=reprojection_error_usingRT(Xw,U,R,T,A_UPnP);
    
    sol(3).Xc=Xc;
    sol(3).Cc=Cc;
    sol(3).R=R;
    sol(3).T=T;
    sol(3).f=f;
    sol(3).betas=[beta1,beta2,beta3];
    sol(3).sc=sc;
    sol(3).error=err(3);
    % -----------   END EFFICIENT EXPLORATION OF MINIMAL EQUATION SETS ---------------------
    
    % ---- Refine the solution iterating over the betas (N=3) ----
    Betas=[sol(3).betas,f];
    Beta0=Betas;
    Beta0(1:3)=Betas(1:3)/sc;
    [Xc_opt,R_opt,T_opt,err_opt,Beta_opt] = UPnP_optimize_betas_gauss_newton(Kernel,Cw,Beta0,Alph,Xw,U,A);
    
    % Just update R,T,Xc if Gauss Newton improves results (which is almost always)
    if (err_opt<err(3)) && (Beta_opt(4) > 0) && (Beta_opt(4) < 1e6)
        sol(3).R=R_opt;
        sol(3).T=T_opt;
        sol(3).Xc=Xc_opt;
        sol(2).betas=Beta_opt(1:3);
        sol(3).f=Beta_opt(4);
        sol(3).error=err_opt;
        err(3)=err_opt;
    end
end


%% 4.- We keep the one with less reprojection error

[min_err,best_solution]=min(err);

Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
T=sol(best_solution).T;
f=sol(best_solution).f;




%% Rotation Calculation Funcion
function [R, T]=getrotT(wpts,cpts)

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
R=U*V';
if det(R)<0
    R=-R;
end
T=ccent'-R*wcent';
%

%% Euclidean Reprojection Error Calculation Function
function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)

n=size(Xw,1);

P=A*[R,T];
Xw_h=[Xw,ones(n,1)];
Urep_=(P*Xw_h')';

%project reference points into the image plane
Urep=zeros(n,2);
Urep(:,1)=Urep_(:,1)./Urep_(:,3);
Urep(:,2)=Urep_(:,2)./Urep_(:,3);

%reprojection error
err_=sqrt((U(:,1)-Urep(:,1)).^2+(U(:,2)-Urep(:,2)).^2);
err=sum(err_)/n;
%

