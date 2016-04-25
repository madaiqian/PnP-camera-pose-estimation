% MAIN Illustrates how to use the UPnP algorithm described in:
%
%   A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer
%   Exhaustive Linearization for Robust Camera Pose and Focal Length Estimation.
%   In IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013.
%
% Copyright (C) <2013>  <A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer>
%                       Institut de Robotica i Inform�tica Industrial (CSIC-UPC)
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

clear all; close all; warning off;

addpath data;
addpath pose_functions;
addpath aux_functions;



Number_eq_sets = 500;   % Number of equations sets to explore
Explore_N3 = true;      % Set this to 'false' to avoid exploration of Kernel dimensionality N=3

%% Generate simulated input data
load_points=0;
if ~load_points
    n = 20;                 % Number of points
    std_noise = 10;         % Noise in the measurements (in pixels)
    focal_lenght = 1000;    % Focal lenght
    outliers = 0;           % Outlier percentage
    [A,point,Rt]=generate_noisy_input_data(n,std_noise,focal_lenght,outliers);
    
    save('input_data_noise.mat','A','point','Rt');
else
    load('input_data_noise.mat','A','point','Rt');
    n=size(point,2);
    draw_noisy_input_data(point);
end


%% Compute ground truth R,T
x3d=zeros(n,4);
x2d=zeros(n,3);
A=A(:,1:3);
for i=1:n
    x3d_h(i,:)=[point(i).Xworld',1];
    x2d_h(i,:)=[point(i).Ximg(1:2)',1];
    x2d_h_without_noise(i,:)=[point(i).Ximg_true(1:2)',1];
    
    % world and camera coordinates
    X3d_world(i,:)=point(i).Xworld';
    X3d_cam(i,:)=point(i).Xcam';
end
[R_true,T_true]=getrotT(X3d_world,X3d_cam);


%% UPnP + Gauss Newton
fprintf('\n---------UPnP-GN--------------\n');

Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

A_UPnP = A;
A_UPnP(1,1) = 0;
A_UPnP(2,2) = 0;

[Rp,Tp,Xc,f,best_solution,sol,opt]=upnp_GN(x3d_h,x2d_h,A_UPnP,Number_eq_sets,Explore_N3);

A_UPnP(1,1) = f;
A_UPnP(2,2) = f;

[error_Rec,error_Rep,error_Linf,error_Rot,error_Tras]=compute_all_errors_rms(Xw,U,X3d_cam,R_true,T_true,Xc,Rp,Tp,A_UPnP);
[error_Rec,error_Rot,error_Tras]=compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,Rp,Tp);

fprintf('Reprojection error UPnP-GN: %.3f\n',error_Rep);
fprintf('Rotation error UPnP-GN: %.3f\n',error_Rot);
fprintf('Translation error UPnP-GN: %.3f\n',error_Tras);
fprintf('Estimated Focal Lenght: %.1f (True %.1f)\n',A_UPnP(1,1),A(1,1));

%draw Results
for i=1:n
    point(i).Xcam_est=Xc(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP Gauss Newton',A,R_true,T_true,h);

