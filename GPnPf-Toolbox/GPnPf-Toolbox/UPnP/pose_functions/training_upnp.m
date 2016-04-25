function [R,T,Xc,best_solution,sol] = training_upnp(x3d_h,x2d_h,A,X3d_cam,R_true,T_true)

Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

%define control points in a world coordinate system (centered on the 3d
%points centroid)
Cw=define_control_points();

%compute alphas (linear combination of the control points to represent the 3d
%points)
Alph=compute_alphas(Xw,Cw);


Mf=compute_M_ver3(U,Alph,A);
Kmf=kernel_noise(Mf,4);


%% 1.-Solve assuming dim(ker(M))=1. with f unknown;------------------------------
Kmf1=Kmf(:,end);

D=compute_constraint_distance_2param_6eq_2unk_f_unk(Kmf1);
dsq=define_distances_btw_control_points();
betas_=(D'*D)\D'*dsq;
beta=sqrt(abs(betas_(1)));
f=sqrt(abs(betas_(2)))/beta;
X1=beta*Kmf1;

% Xci = beta*Vi1
% Yci = beta*Vi2
% Zci = f*beta*Vi3
X1(3:3:size(Kmf1,1))=X1(3:3:size(Kmf1,1))*f;

[Cc,Xc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw);

[R,T]=getrotT(Xw,Xc);  %solve exterior orientation
err(1)=reprojection_error_usingRT(Xw,U,R,T,A);
[error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,R,T);

sol(1).Xc=Xc;
sol(1).Cc=Cc;
sol(1).R=R;
sol(1).T=T;
sol(1).error_Rec= error_Rec;
sol(1).error_Rot= error_Rot;
sol(1).error_Tras= error_Tras;
sol(1).error_Rep= err(1);



%% 2.-Solve assuming dim(ker(M))=2. with f unknown;------------------------------
Kmf1=Kmf(:,end-1);
Kmf2=Kmf(:,end);

% control points distance constraint
D=compute_constraint_distance_3param_6eq_6unk_f_unk(Kmf1,Kmf2);
dsq=define_distances_btw_control_points();
betas_=D\dsq;

beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2));

solutions= generate_all_possible_solutions_for_f_unk(betas_, 2);
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
    
    [~,Xc]=compute_norm_sign_scaling_factor(X2_aux,Cw,Alph,Xw);
    
    [R,T]=getrotT(Xw,Xc);  %solve exterior orientation
    error_aux=reprojection_error_usingRT(Xw,U,R,T,A);
    [error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,R,T);
    
    solutions_error(i)= error_Rec;
    
    if ((i==1) || (error_aux < smallest_error))
        beta1= aux_beta1;
        beta2= aux_beta2;
        f= aux_f;
        smallest_error=error_aux;
    end
end

% Xci = beta*Vi1
% Yci = beta*Vi2
% Zci = f*beta*Vi3
X2=beta1*Kmf1+beta2*Kmf2;
X2(3:3:size(Kmf1,1))=X2(3:3:size(Kmf1,1))*f;
[Cc,Xc]=compute_norm_sign_scaling_factor(X2,Cw,Alph,Xw);

[R,T]=getrotT(Xw,Xc);  %solve exterior orientation
err(2)=reprojection_error_usingRT(Xw,U,R,T,A);
[error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,R,T);

sol(2).Xc=Xc;
sol(2).Cc=Cc;
sol(2).R=R;
sol(2).T=T;
sol(2).solutions_error=solutions_error;
sol(2).error_Rec= error_Rec;
sol(2).error_Rot= error_Rot;
sol(2).error_Tras= error_Tras;
sol(2).error_Rep= err(2);


%% 3.-Solve assuming dim(ker(M))=3------------------------------------------
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
lambda_kernel_size= 5;

solution_counter= 1;
%-----------   All solutions   ----------------------------
[solutions_to_lambda, lambda_equation_combination]= generate_all_possible_solutions(lambdas_, lambda_kernel_size);
num_solutions_to_lambda= size(solutions_to_lambda,1);
for j=1:num_solutions_to_lambda
    lambda(1)= solutions_to_lambda(j,1);
    lambda(2)= solutions_to_lambda(j,2);
    lambda(3)= solutions_to_lambda(j,3);
    lambda(4)= solutions_to_lambda(j,4);
    lambda(5)= solutions_to_lambda(j,5);
    
    betas_=lambda(1)*Kd(:,1)+lambda(2)*Kd(:,2)+lambda(3)*Kd(:,3)+lambda(4)*Kd(:,4)+lambda(5)*Kd(:,5);
    
    [solutions_to_beta, beta_equation_combination]= generate_all_possible_solutions_for_f_unk(betas_, 3); % B1, B2, B3
    num_solutions_to_beta= size(solutions_to_beta,1);
    
    if(j==1)
        case_error_Rep= zeros(num_solutions_to_lambda*num_solutions_to_beta, 1);
        case_error_Rec= zeros(num_solutions_to_lambda*num_solutions_to_beta, 1);
        case_error_Rot= zeros(num_solutions_to_lambda*num_solutions_to_beta, 1);
        case_error_Trans= zeros(num_solutions_to_lambda*num_solutions_to_beta, 1);
        lambdas= zeros(num_solutions_to_lambda*num_solutions_to_beta, 5);
        lambda_eq=  zeros(num_solutions_to_lambda*num_solutions_to_beta, 5);
        betas_and_f=  zeros(num_solutions_to_lambda*num_solutions_to_beta, 4);
        beta_eq=  zeros(num_solutions_to_lambda*num_solutions_to_beta, 4);
    end
    
    for i=1:num_solutions_to_beta
        aux_beta1= solutions_to_beta(i,1);
        aux_beta2= solutions_to_beta(i,2);
        aux_beta3= solutions_to_beta(i,3);
        aux_f= solutions_to_beta(i,4);
        
        % Xci = beta*Vi1
        % Yci = beta*Vi2
        % Zci = f*beta*Vi3
        X3=aux_beta1*Kmf1+aux_beta2*Kmf2+aux_beta3*Kmf3;
        X3_aux=X3;
        X3_aux(3:3:size(Kmf1,1))=X3(3:3:size(Kmf1,1))*aux_f;
        
        [~,Xc]=compute_norm_sign_scaling_factor(X3_aux,Cw,Alph,Xw);
        
        [R,T]=getrotT(Xw,Xc);  %solve exterior orientation
        error_aux=reprojection_error_usingRT(Xw,U,R,T,A);
        
        [error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,R,T);
        
        if (((j==1)&&(i==1)) || (error_aux < smallest_error))
            beta1= aux_beta1;
            beta2= aux_beta2;
            beta3= aux_beta3;
            f= aux_f;
            smallest_error=error_aux;
        end
        
        case_error_Rep(solution_counter,:)= error_aux;
        case_error_Rec(solution_counter,:)= error_Rec;
        case_error_Rot(solution_counter,:)= error_Rot;
        case_error_Trans(solution_counter,:)= error_Tras;
        lambdas(solution_counter,:)= lambda;
        lambda_eq(solution_counter,:)= lambda_equation_combination(j,:);
        betas_and_f(solution_counter,:)= [aux_beta1, aux_beta2, aux_beta3, aux_f];
        beta_eq(solution_counter,:)= beta_equation_combination(i,:);
        solution_counter= solution_counter+1;
        
    end
end
%----------   END All solutions   -----------------------


% Xci = beta*Vi1
% Yci = beta*Vi2
% Zci = f*beta*Vi3
X3=beta1*Kmf1+beta2*Kmf2+beta3*Kmf3;
X3_aux=X3;
X3_aux(3:3:size(Kmf1,1))=X3(3:3:size(Kmf1,1))*f;

[Cc,Xc]=compute_norm_sign_scaling_factor(X3_aux,Cw,Alph,Xw);

[R,T]=getrotT(Xw,Xc);  % Solve exterior orientation
err(3)=reprojection_error_usingRT(Xw,U,R,T,A);

[error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,R,T);

sol(3).Xc=Xc;
sol(3).Cc=Cc;
sol(3).R=R;
sol(3).T=T;
sol(3).error_Rec= error_Rec;
sol(3).error_Rot= error_Rot;
sol(3).error_Tras= error_Tras;
sol(3).error_Rep= err(3);
sol(3).case_error_Rep= case_error_Rep;
sol(3).case_error_Rec= case_error_Rec;
sol(3).case_error_Rot= case_error_Rot;
sol(3).case_error_Trans= case_error_Trans;
sol(3).lambdas= lambdas;
sol(3).lambda_eq= lambda_eq;
sol(3).betas_and_f= betas_and_f;
sol(3).beta_eq= beta_eq;
end

[~,best_solution]=min(err);
Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
T=sol(best_solution).T;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)
        
%clear all; close all; load reprojection_error_usingRT;
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
