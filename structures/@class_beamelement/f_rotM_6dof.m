%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_rotM_6dof:  class file is part of nlFEM, class beamelement
%               rotation matrix for 6DOF stiffness matrix
%     Author:   Klaus Seywald
%               klaus.seywald@mytum.de
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = f_rotM_6dof(obj,a,b,c)

Lx=[1       0       0
    0   cos(a)  sin(a)
    0   -sin(a) cos(a)];

Ly=[cos(b) 0 -sin(b)
    0      1    0
    sin(b)  0   cos(b)];

Lz=[cos(c) sin(c)   0
    -sin(c) cos(c)  0
    0           0   1];

L=Lz*Ly*Lx;

Z=zeros(3);

T=[L Z Z Z
   Z L Z Z
   Z Z L Z
   Z Z Z L];

end
