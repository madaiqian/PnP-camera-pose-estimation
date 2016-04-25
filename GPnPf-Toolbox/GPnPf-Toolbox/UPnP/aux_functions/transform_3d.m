function X2=transform_3d(R,X1)

X1_h=[X1;1];

X2_h=R*X1_h;

X2(1,1)=X2_h(1);
X2(2,1)=X2_h(2);
X2(3,1)=X2_h(3);
