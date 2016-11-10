%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [X_nxt,Euler_nxt,V_K_nxt,omega_K_nxt,V_kin,V_e,alpha_k,beta_k,omega_K_dot,structural_acc]=sixDoF_rigid_body(mass,JJ,JJdot,X,Euler,V_K,omega_K,F,M,dt)
% experiment with time integration schemes

V_K_dot_struct=[1/mass*F(1) 1/mass*F(2) 1/mass*F(3)];

V_K_dot=([1/mass*F(1)  1/mass*F(2) 1/mass*F(3)]-cross(omega_K,V_K));

omega_K_dot=(JJ^(-1)*([M(1) M(2) M(3)]-cross(omega_K,JJ*omega_K')-(JJdot*omega_K')')')';

% time integration step
V_K_nxt=V_K+V_K_dot*dt;


% time integration step
omega_K_nxt=omega_K+omega_K_dot*dt;

V_kin=norm(V_K);
alpha_k=atan(V_K(3)/V_K(1));
beta_k=atan(-V_K(2)/V_K(1));

%% attitude differential equation
T=[1 0 0
    0 1 0
    0 0 1];

Euler_dot=(([1    sin(Euler(1))*tan(Euler(2))     cos(Euler(1))*tan(Euler(2));
    0    cos(Euler(1))                   -sin(Euler(1));
    0   sin(Euler(1))/cos(Euler(2))      cos(Euler(1))/cos(Euler(2));]*omega_K'))';

% time integration step
Euler_nxt=Euler+Euler_dot*dt;



M_BO=[  cos(Euler(3))*cos(Euler(2))                                                 sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
    cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
    cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];


V_e=(M_BO'*V_K_nxt')';
X_nxt=X+(M_BO'*V_K')'*dt;

structural_acc=(T*V_K_dot_struct');
end
