%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [M_fa,K_fa,F_fa] =compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,M_Body,K_body,modeshapes,n_modes,xdot_now)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    % compute nodal coordinates in inertial system
    modeshapes_rom=real(modeshapes)*[eye(n_modes); zeros(size(modeshapes,1)-n_modes,n_modes)];

    Mm=modeshapes_rom(1:end,7:n_modes)'*M_Body(1:end,1:end)*modeshapes_rom(1:end,7:n_modes);
    Km=modeshapes_rom(1:end,7:n_modes)'*K_body(1:end,1:end)*modeshapes_rom(1:end,7:n_modes);

    M_inert=A_Bar*M_Body*A_Bar';
    Mass=I_Hat'*M_inert*I_Hat;
    Inertia=b_Hat_Skew'*M_inert*b_Hat_Skew;

    zw1=I_Hat'*A_Bar*M_Body;
    zw2=b_Hat_Skew'*A_Bar*M_Body;
    zw3=modeshapes_rom(1:end,7:n_modes)'*M_Body*A_Bar'*I_Hat;
    zw4=modeshapes_rom(1:end,7:n_modes)'*M_Body*A_Bar'*b_Hat_Skew;

    M_fa=[Mass                           zeros(3,3)                             zeros(3,length(Mm));
        zeros(3,3)                       Inertia                                zeros(3,length(Mm));
        zw3                             zw4                                    Mm];

    K_fa=M_fa*0;
       
    zw5=I_Hat'*A_Bar*K_body*modeshapes_rom(1:end,7:n_modes)*0;
    zw6=b_Hat_Skew'*A_Bar*K_body*modeshapes_rom(1:end,7:n_modes)*0;
    zw7=Km;

    K_fa(:,7:end)=[zw5;zw6;zw7];
    
    f_r=I_Hat'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*A_Bar*modeshapes_rom(1:end,7:n_modes)*xdot_now(7:end));
    f_om=b_Hat_Skew'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*A_Bar*modeshapes_rom(1:end,7:n_modes)*xdot_now(7:end));
    f_de=modeshapes_rom(1:end,7:end)'*A_Bar'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*A_Bar*modeshapes_rom(1:end,7:n_modes)*xdot_now(7:end));
    
    F_fa=[-f_r;-f_om;-f_de];
end

