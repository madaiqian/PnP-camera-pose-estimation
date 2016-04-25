addpath('GPnPf-Toolbox/GPnPf-Toolbox');
addpath(genpath('GPnPf-Toolbox\GPnPf-Toolbox\rpnp1.0'));
addpath('p4pf');


m=4000;
tol=0;
for i=1:m
    i
    npt=4;
    f= rand(1)*1800+200; %random focal length in [200,2000]
    %f= 800;
    
    % generate 3d coordinates in camera space
    Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
    t= mean(Xc,2);
    R= rodrigues(randn(3,1));
    XXw= inv(R)*(Xc-repmat(t,1,npt));

    % projection
    xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
    
    M=XXw;m=xx;
    [min_p35p]=p3_5p(M,m,f);
    [min_Bujnak]=p4pf_func(M,m,f);
    [f_R, f_A, f_AT ,f_ATc] = RPnP_f_Test(M,m);
    tol=tol+1;     %×éºÅ
    ans(tol,1)=min_p35p;
    ans(tol,2)=min_Bujnak;  
    if length(f_R)>0 && length(f_A)>0 && length(f_AT)>0 && length(f_ATc)>0
         ans(tol,3) = min(abs(f_R-f)/f);
         ans(tol,4) = min(abs(f_A-f)/f);
         ans(tol,5) = min(abs(f_AT-f)/f);
         ans(tol,6) = min(abs(f_ATc-f)/f);
    end
end