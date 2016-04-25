% COMPUTE_A_AND_B_GAUSS_NEWTON  
%
%       A: Jacobian
%       b: error vector
%
function [A,b]=UPnP_compute_A_and_b_Gauss_Newton(current_betas,f,rho,L)


A=zeros(6,4);
b=zeros(6,1);

cb=current_betas;

B=[cb(1)*cb(1),
   cb(1)*cb(2),
   cb(1)*cb(3),
   cb(2)*cb(2),
   cb(2)*cb(3),
   cb(3)*cb(3),
   f*f*cb(1)*cb(1),
   f*f*cb(1)*cb(2),
   f*f*cb(1)*cb(3),
   f*f*cb(2)*cb(2),
   f*f*cb(2)*cb(3),
   f*f*cb(3)*cb(3)];

for i=1:6
    A(i,1) = 2*cb(1)*L(i,1) + cb(2)*L(i,2) + cb(3)*L(i,3) + f*f*( 2*cb(1)*L(i,7) + cb(2)*L(i,8) + cb(3)*L(i,9) );               % Beta 1
    A(i,2) = cb(1)*L(i,2) + 2*cb(2)*L(i,4) + cb(3)*L(i,5) + f*f*( cb(1)*L(i,8) + 2*cb(2)*L(i,10) + cb(3)*L(i,11) );             % Beta 2
    A(i,3) = cb(1)*L(i,3) + cb(2)*L(i,5) + 2*cb(3)*L(i,6) + f*f*( cb(1)*L(i,9) + cb(2)*L(i,11) + 2*cb(3)*L(i,12) );             % Beta 3
    A(i,4) = 2*f*( cb(1)*cb(1)*L(i,7) + cb(1)*cb(2)*L(i,8) + cb(1)*cb(3)*L(i,9) + cb(2)*cb(2)*L(i,10) + cb(2)*cb(3)*L(i,11) + cb(3)*cb(3)*L(i,12) );    % f
        
    b(i) = (rho(i)-L(i,:)*B);
end

%% Finished
