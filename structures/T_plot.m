%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function T_plot = T_plot(a,b,c)

Lx=[1       0       0
     0   cos(a)  sin(a)
     0   -sin(a) cos(a)];
 
Ly=[cos(b) 0 -sin(b)
    0      1    0
    sin(b)  0   cos(b)];

Lz=[cos(c) sin(c)   0
    -sin(c) cos(c)  0
     0        0     1];
 
T_plot=Lz*Ly*Lx;

end
