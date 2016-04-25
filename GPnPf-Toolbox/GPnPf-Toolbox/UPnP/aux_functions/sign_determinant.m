function s=sign_determinant(C)

c0=C(4,:)';
c1=C(1,:)';
c2=C(2,:)';
c3=C(3,:)';

v1=c1-c0;
v2=c2-c0;
v3=c3-c0;

M=[v1,v2,v3];
detM=det(M);
s=sign(detM);
