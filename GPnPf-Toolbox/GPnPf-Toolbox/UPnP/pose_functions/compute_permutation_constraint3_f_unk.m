function K=compute_permutation_constraint3_f_unk(V)

%[B11,B12,...,B33,Bf11,Bf12,...,Bf33]=lambda1*v1+lambda2*v2+lambda3*v3

N=5;size(V,2); %dimension of the kernel
n=3; %dimension of Bij
idx=[1 2 3; 2 4 5; 3 5 6];

%1.-Generation of the first set of equations Bii.Bjj=Bij.Bii  (n(n-1)/2 eqs).

% Size as defined in Ansar & Daniilidis paper, in our case we have to bear
% in mind that the total will be 4*nrowsK, because we have to produce the
% combinations with the f; these will be [iijk], [fiijk], [iifjk] & [fiifjk]
nrowsK=n*n*(n-1)/2;%disp(nrowsK);
ncolsK=N*(N+1)/2;%disp(ncolsK);
K=zeros(nrowsK*4,ncolsK);

t=0;
for i=1:n
    for j=1:n
        for k=1:n
            if (i~=j)&&(i~=k)
                offset=1;
                t=t+1;
                for a=1:N
                    for b=a:N
                        if a==b
                            % [iijk]->[ijik] case
                            K(t,offset)= V(idx(i,i),a)*V(idx(j,k),a)-V(idx(i,j),a)*V(idx(i,k),a);
                            
                            % We add 6 because the Bf's start at that 
                            % position in our case because we have 12 unknowns
                            % [fiijk]->[fijik] case; 
                            K(t+nrowsK,offset)= V(idx(i,i)+6,a)*V(idx(j,k),a)-V(idx(i,j)+6,a)*V(idx(i,k),a);
                            
                            % [fiijk]->[fijik] case; 
                            K(t+nrowsK*2,offset)= V(idx(i,i),a)*V(idx(j,k)+6,a)-V(idx(i,j),a)*V(idx(i,k)+6,a);
                            
                            % [fiifjk]->[fijfik] case; 
                            K(t+nrowsK*3,offset)= V(idx(i,i)+6,a)*V(idx(j,k)+6,a)-V(idx(i,j)+6,a)*V(idx(i,k)+6,a);
                        else
                            % [iijk]->[ijik] case
                            K(t,offset)= V(idx(i,i),a)*V(idx(j,k),a)-V(idx(i,j),a)*V(idx(i,k),a)+...
                                         V(idx(i,i),b)*V(idx(j,k),b)-V(idx(i,j),b)*V(idx(i,k),b);
                            
                            % [fiijk]->[fijik] case; 
                            K(t+nrowsK,offset)= V(idx(i,i)+6,a)*V(idx(j,k),a)-V(idx(i,j)+6,a)*V(idx(i,k),a)+...
                                                V(idx(i,i)+6,b)*V(idx(j,k),b)-V(idx(i,j)+6,b)*V(idx(i,k),b);
                            
                            % [iifjk]->[ijfik] case; 
                            K(t+nrowsK*2,offset)= V(idx(i,i),a)*V(idx(j,k)+6,a)-V(idx(i,j),a)*V(idx(i,k)+6,a)+...
                                                  V(idx(i,i),b)*V(idx(j,k)+6,b)-V(idx(i,j),b)*V(idx(i,k)+6,b);
                            
                            % [fiifjk]->[fijfik] case; 
                            K(t+nrowsK*3,offset)= V(idx(i,i)+6,a)*V(idx(j,k)+6,a)-V(idx(i,j)+6,a)*V(idx(i,k)+6,a)+...
                                                  V(idx(i,i)+6,b)*V(idx(j,k)+6,b)-V(idx(i,j)+6,b)*V(idx(i,k)+6,b);
                        end
                        offset=offset+1;
                    end
                end
                        
%                 fprintf('[%d%d%d%d]->[%d%d%d%d]\n',i,i,j,k,i,j,i,k);
%                 fprintf('[f%d%d%d%d]->[%d%df%d%d]\n',i,i,j,k,i,j,i,k);
%                 fprintf('[%d%df%d%d]->[f%d%d%d%d]\n',i,i,j,k,i,j,i,k);
%                 fprintf('[f%d%df%d%d]->[f%d%df%d%d]\n',i,i,j,k,i,j,i,k);
%                 fprintf('t:%d\t offset:%d\n',t,offset);
            end
        end
    end
end

