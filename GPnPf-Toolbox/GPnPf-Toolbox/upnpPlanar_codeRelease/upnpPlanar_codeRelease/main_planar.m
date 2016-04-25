% MAIN Illustrates how to use the UPnP algorithm described in:
%
%   A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer
%   Exhaustive Linearization for Robust Camera Pose and Focal Length Estimation.
%   In IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013.
%
% Copyright (C) <2013>  <A.Penate-Sanchez, J.Andrade-Cetto, F.Moreno-Noguer>
%                       Institut de Robotica i Informatica Industrial (CSIC-UPC)
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


clear all; close all;
addpath data; addpath tests; addpath method;


%% 1.- Generate simulated input data

n = 100; % Number of points
std_noise = 5; % Noise in the measurements (in pixels)
focal_lenght = 1500; % Focal lenght
outliers = 0; % Outlier percentage

inclination = 2; z_offset = 8;
[A,point,Rt] = generate_planar_noisy_input_data(n,std_noise,focal_lenght,outliers,[1 -1 z_offset+inclination; -1 -1 z_offset+inclination; -1 1 z_offset-inclination], 'donotplot');


%% 2.- Input format

x3d = zeros(n,4);
x2d = zeros(n,3); 
A = A(:,1:3);
for i=1:n
    x3d_h(i,:) = [point(i).Xworld',1]; 
    x2d_h(i,:) = [point(i).Ximg(1:2)',1];
    x2d_h_without_noise(i,:) = [point(i).Ximg_true(1:2)',1];

    % World and camera coordinates
    X3d_world(i,:) = point(i).Xworld';
    X3d_cam(i,:) = point(i).Xcam';
end
[R_true,T_true] = getrotT(X3d_world,X3d_cam);


%% 3.- UPnP planar

fprintf('\n--------- UPnP planar --------------\n');

Xw = x3d_h(:,1:3);
U = x2d_h(:,1:2);

[Rp,Tp,Xc,f,best_solution,sol,opt] = upnp_planar(x3d_h,x2d_h,A);

% Compute error
UPnP_A = A; UPnP_A(1,1) = f; UPnP_A(2,2) = f;

[error_Rec,error_Rep,error_Rot,error_Tras] = compute_all_errors_rms(Xw,U,X3d_cam,R_true,T_true,Xc,Rp,Tp,UPnP_A);
[error_Rec,error_Rot,error_Tras] = compute_all_errors_relative(X3d_cam,R_true,T_true,Xc,Rp,Tp);

fprintf('Error Rot UPnP-GN: %.3f\n',error_Rot);
fprintf('Error Tras UPnP-GN: %.3f\n',error_Tras);
fprintf('Error Rec UPnP-GN: %.3f\n',error_Rec);
fprintf('F: %.3f\n',UPnP_A(1,1));



