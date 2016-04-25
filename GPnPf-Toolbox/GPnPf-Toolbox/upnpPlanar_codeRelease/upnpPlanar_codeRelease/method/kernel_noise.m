function K=kernel_noise(M,dimker)
%clear all; close all; load kernel_noise;

%tic;
%[U1,S1,V1]=svd(M);

%compute dimension of the kernel=ncols(M)-rank(M)
% [nrows_M,ncols_M]=size(M);
% rank_M=rank(M);
% dim_kernel=ncols_M-rank_M;
% K=V(:,end-dim_kernel+1:end);


%if dim_kernel==0
%    K1=V1(:,end-dimker+1:end);
%end

%t=toc;

%by doing this way, we do not have to compute SVD and it is way faster
MtM=M'*M;
[V,S]=eig(MtM);

K=V(:,dimker:-1:1);