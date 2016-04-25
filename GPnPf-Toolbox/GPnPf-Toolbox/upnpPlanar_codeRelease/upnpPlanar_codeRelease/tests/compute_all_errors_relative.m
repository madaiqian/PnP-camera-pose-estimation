function [error_rec,error_rot,error_trans]=compute_all_errors_relative(Xc_true,R_true,T_true,Xc,R,T)

n=size(Xc_true,1);

% Reconstruction error
for i=1:n
   err_rel(i)=norm(Xc_true(i,:)-Xc(i,:))/norm(Xc_true(i,:));
end
error_rec=mean(err_rel);

% Rotation error 
q_true=rot_to_quaternion(R_true);
q=rot_to_quaternion(R);
error_rot1=norm(q_true-q)/norm(q_true);
error_rot2=norm(q_true+q)/norm(q_true);
error_rot=min(error_rot1,error_rot2);

% Translation error
error_trans=norm(T_true-T)/norm(T_true);