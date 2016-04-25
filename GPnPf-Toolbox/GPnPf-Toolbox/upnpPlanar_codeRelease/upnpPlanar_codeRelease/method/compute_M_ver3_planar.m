function M=compute_M_ver3_planar(U,Alph,A)

% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

n=size(Alph,1); %number of 3d points

% fu=A(1,1);
% fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M=2*n;
ncols_M=9;
M=zeros(nrows_M,ncols_M);



for i=1:n
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    
    ui=U(i,1);
    vi=U(i,2);
    
    %generate submatrix M
    M_=[a1, 0, a1*(u0-ui), a2, 0, a2*(u0-ui), a3, 0, a3*(u0-ui);
        0, a1, a1*(v0-vi), 0, a2, a2*(v0-vi), 0, a3, a3*(v0-vi)];
    
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=M_;
     
end

