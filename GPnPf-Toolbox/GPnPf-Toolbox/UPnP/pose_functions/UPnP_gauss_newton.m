% GAUSS_NEWTON  
%
%       K 12x4:  Matrix of kernel vectors
%       Cw 4x3:  Coordinates of the control points in world reference
%       Beta0: initial estimation of Beta's
%
function [Beta_opt,err,iter]=UPnP_gauss_newton(K,Cw,Beta0)

L = compute_L6_12(K);
rho = UPnP_compute_rho(Cw);

current_betas = Beta0;
  
n_iterations=50; % max number of iterations. Usually 4-5 is enough.
for k=1:n_iterations
    [A,b] = UPnP_compute_A_and_b_Gauss_Newton(current_betas,current_betas(4),rho,L);
    dbeta = ((A'*A)\A')*b;
    current_betas = current_betas + dbeta';
    error(k) = b'*b;
    iter(k).betas = current_betas;
    iter(k).error = error(k);
end

[Y,index] = min([iter([1:end]).error]);

Beta_opt = iter(index).betas;
err = iter(index).error;

% The focal lenght must be positive
Beta_opt(4) = abs(Beta_opt(4));


