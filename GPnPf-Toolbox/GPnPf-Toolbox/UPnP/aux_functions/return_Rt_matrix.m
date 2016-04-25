function R=return_Rt_matrix(alpha,beta,gamma,tx,ty,tz)

R(1,1)=cos(alpha)*cos(gamma)-cos(beta)*sin(alpha)*sin(gamma);
R(2,1)=cos(gamma)*sin(alpha)+cos(alpha)*cos(beta)*sin(gamma);
R(3,1)=sin(beta)*sin(gamma);
R(4,1)=0;

R(1,2)=-cos(beta)*cos(gamma)*sin(alpha)-cos(alpha)*sin(gamma);
R(2,2)=cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma);
R(3,2)=cos(gamma)*sin(beta);
R(4,2)=0;

R(1,3)=sin(alpha)*sin(beta);
R(2,3)=-cos(alpha)*sin(beta);
R(3,3)=cos(beta);
R(4,3)=0;

R(:,4)=[tx,ty,tz,1]';

