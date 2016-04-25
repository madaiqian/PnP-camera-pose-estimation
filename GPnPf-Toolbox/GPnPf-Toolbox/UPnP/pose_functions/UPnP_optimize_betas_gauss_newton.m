% COMPUTE_BETAS_GAUSS_NEWTON  
%
%       Km: vector of the kernel
%       Cw: position of the control point in world coordinates
%       Beta0: initial guess of the betas and f (Beta0(4))
%
function [Xc_opt,R_opt,T_opt,err_opt,Beta_opt]=UPnP_optimize_betas_gauss_newton(Km,Cw,Beta0,Alph,Xw,U,A)

A(1,1) = 0;
A(2,2) = 0;

n=size(Beta0,2);

[Beta_opt,err,iter]=UPnP_gauss_newton(Km,Cw,Beta0);

% Extract control point camera coordinates from Betas and Kernel
X=zeros(12,1);
for i=1:(n-1)
   X=X+Beta_opt(i)*Km(:,i); 
end
f = Beta_opt(4);

X(3:3:12)=X(3:3:12)*f;

Cc=zeros(4,3);
for i=1:4
    Cc(i,:)=X(3*i-2:3*i);
end


% Check sign of the determinant (keep orientation of the control points)
s_Cw=sign_determinant(Cw);
s_Cc=sign_determinant(Cc);
Cc=Cc*(s_Cw/s_Cc);

Xc_opt=Alph*Cc;
[R_opt,T_opt]=getrotT(Xw,Xc_opt);  % Solve exterior orientation

A(1,1) = f;
A(2,2) = f;
[err_opt,Urep_opt]=reprojection_error_usingRT(Xw,U,R_opt,T_opt,A);

