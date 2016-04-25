function P=compute_constraint_distance_2param_6eq_2unk_f_unk_planar(m1)

% Copyright (C) <2010>  <Adrián Peñate-Sánchez, Francesc Moreno-Noguer>
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

%redefine variables name, for compatibility with maple
m_1=m1(1); 
m_2=m1(2); 
m_3=m1(3); 
m_4=m1(4); 
m_5=m1(5); 
m_6=m1(6);
m_7=m1(7); 
m_8=m1(8); 
m_9=m1(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = (m_4 ^ 2);
t4 = (m_1 ^ 2);
t5 = (m_5 ^ 2);
t8 = (m_2 ^ 2);
t10 = (m_6 ^ 2);
t13 = (m_3 ^ 2);
t15 = (m_7 ^ 2);
t18 = (m_8 ^ 2);
t22 = (m_9 ^ 2);
P(1,1) = t1 - 2 * m_4 * m_1 + t4 + t5 - 2 * m_5 * m_2 + t8;
P(1,2) = t10 - 2 * m_6 * m_3 + t13;
P(2,1) = t15 - 2 * m_7 * m_1 + t4 + t18 - 2 * m_8 * m_2 + t8;
P(2,2) = t22 - 2 * m_9 * m_3 + t13;
P(3,1) = t15 - 2 * m_7 * m_4 + t1 + t18 - 2 * m_8 * m_5 + t5;
P(3,2) = t22 - 2 * m_9 * m_6 + t10;



