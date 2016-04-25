function  Ximg=project_3d_2d(A,Xcam)

Xcam_h=[Xcam;1];

Ximg_h=A*Xcam_h;

Ximg(1,1)=Ximg_h(1)/Ximg_h(3);
Ximg(2,1)=Ximg_h(2)/Ximg_h(3);
