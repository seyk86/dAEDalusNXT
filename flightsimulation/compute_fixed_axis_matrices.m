%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [M_fa,K_fa,F_fa] =compute_fixed_axis_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,M_Body,K_body,xdot_now, fixed_frame_bcs )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    % compute nodal coordinates in inertial system  
    M_inert=A_Bar*M_Body*A_Bar';
    
    K_inert=A_Bar*K_body*A_Bar';
    Mass=I_Hat'*M_inert*I_Hat;
    Inertia=b_Hat_Skew'*M_inert*b_Hat_Skew;

    zw1=I_Hat'*A_Bar*M_Body;
    zw2=b_Hat_Skew'*A_Bar*M_Body;
    zw3=M_Body*A_Bar'*I_Hat;
    zw4=M_Body*A_Bar'*b_Hat_Skew;
    

    M_fa=[Mass                     I_Hat'*M_inert*b_Hat_Skew          zw1(:,fixed_frame_bcs);
        b_Hat_Skew'*M_inert*I_Hat  Inertia                            zw2(:,fixed_frame_bcs);
        zw3(fixed_frame_bcs,:)     zw4(fixed_frame_bcs,:)             M_Body(fixed_frame_bcs,fixed_frame_bcs)];
    
%     dzw1=I_Hat'*dM_inert*A_Bar;
%     dzw2=-b_Hat_Skew'*dM_inert*A_Bar;
%     dzw3=A_Bar'*dM_inert*I_Hat;
%     dzw4=-A_Bar'*dM_inert*b_Hat_Skew;
%     dM_tot=[Mass*0                     -I_Hat'*dM_inert*b_Hat_Skew          dzw1(:,fixed_frame_bcs);
%         -b_Hat_Skew'*dM_inert*I_Hat  b_Hat_Skew'*dM_inert*b_Hat_Skew                            dzw2(:,fixed_frame_bcs);
%         dzw3(fixed_frame_bcs,:)                dzw4(fixed_frame_bcs,:)    M_Body(fixed_frame_bcs,fixed_frame_bcs)*0]*0;
    
    K_fa=M_fa*0;

    zw5=I_Hat'*A_Bar*K_body*0;
    zw6=b_Hat_Skew'*A_Bar*K_body*0;
    zw7=K_body;
    
    K_fa(:,7:end)=[zw5(:,fixed_frame_bcs);zw6(:,fixed_frame_bcs);zw7(fixed_frame_bcs,fixed_frame_bcs)];

    zwxdot_now=xdot_now(7:end);
    xdot_now_full=zeros(length(K_body),1);
    for nn=1:length(fixed_frame_bcs)
        xdot_now_full(fixed_frame_bcs(nn))=zwxdot_now(nn);
    end

    f_r=I_Hat'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*xdot_now_full);
    f_om=b_Hat_Skew'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*xdot_now_full);
    f_de=A_Bar'*M_inert*(omega_Hat_Skew*omega_Hat_Skew*b+2*omega_Hat_Skew*A_Bar*xdot_now_full);
    F_fa=[-f_r;-f_om;-f_de(fixed_frame_bcs)];

end

