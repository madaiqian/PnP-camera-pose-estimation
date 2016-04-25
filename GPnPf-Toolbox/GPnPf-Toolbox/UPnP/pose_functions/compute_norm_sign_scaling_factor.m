function [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw)

n=size(Xw,1); %number of data points

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(4,3);
for i=1:4
    Cc_(i,:)=X1(3*i-2:3*i);
end

%position of reference points in camera coordinates
Xc_=Alph*Cc_;

%compute distances in world coordinates w.r.t. the centroid
centr_w=mean(Xw);
centroid_w=repmat(centr_w,[n,1]);
tmp1=Xw-centroid_w;
dist_w=sqrt(sum(tmp1.^2,2));

%compute distances in camera coordinates w.r.t. the centroid
centr_c=mean(Xc_);
centroid_c=repmat(centr_c,[n,1]);
tmp2=Xc_-centroid_c;
dist_c=sqrt(sum(tmp2.^2,2));
 
%least squares solution for the scale factor
sc=1/(inv(dist_c'*dist_c)*dist_c'*dist_w);

%scale position of the control points
Cc=Cc_/sc;

%rescaled position of the reference points
Xc=Alph*Cc;

%change the sign if necessary. z negative is not possible in camera coordinates
neg_z=find(Xc(:,3)<0);
if size(neg_z,1)>=1
    sc=-sc;
    Xc=Xc*(-1);
end




        
        
