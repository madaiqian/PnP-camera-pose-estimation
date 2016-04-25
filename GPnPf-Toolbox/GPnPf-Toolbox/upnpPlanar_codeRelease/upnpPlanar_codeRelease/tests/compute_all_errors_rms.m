function [error_rec,error_repr,error_rot,error_trans]=compute_all_errors_rms(Xw,U,Xc_true,R_true,T_true,Xc,R,T,A)

n=size(Xc_true,1);

% Reconstruction error
N=3;
for i=1:n
   err_rms(i)=sqrt((1/N)*dot(Xc_true(i,:)-Xc(i,:),Xc_true(i,:)-Xc(i,:)));
end
error_rec=mean(err_rms);

% Reprojection error
error_repr=reprojection_error_usingRT(Xw,U,R,T,A);

% Rotation error 
q_true=rot_to_quaternion(R_true);
q=rot_to_quaternion(R);
N=4;
error_rot1=sqrt((1/N)*dot(q_true-q,q_true-q));
error_rot2=sqrt((1/N)*dot(q_true+q,q_true+q));
error_rot=min(error_rot1,error_rot2);

% Translation error
N=3;
error_trans=sqrt((1/N)*dot(T_true-T,T_true-T));

