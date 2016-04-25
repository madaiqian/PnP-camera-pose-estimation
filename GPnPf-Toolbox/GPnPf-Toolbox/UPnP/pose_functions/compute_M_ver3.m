function M=compute_M_ver3(U,Alph,A)

n=size(Alph,1); %number of 3d points

% fu=A(1,1);
% fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M=2*n;
ncols_M=12;
M=zeros(nrows_M,ncols_M);

for i=1:n
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    a4=Alph(i,4);
    
    ui=U(i,1);
    vi=U(i,2);
    
    %generate submatrix M
    M_=[a1, 0, a1*(u0-ui), a2, 0, a2*(u0-ui), a3, 0, a3*(u0-ui), a4, 0, a4*(u0-ui);
        0, a1, a1*(v0-vi), 0, a2, a2*(v0-vi), 0, a3, a3*(v0-vi), 0, a4, a4*(v0-vi)];
    
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=M_;
     
end

