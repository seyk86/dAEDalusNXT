%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices(Euler,pqr,node_coords,nodal_deflections,fixed_node)
    n_dof=length(nodal_deflections);
    I_Hat=[];
    for i=1:n_dof/6
        I_Hat=[I_Hat; eye(3,3);zeros(3,3)];
    end
    
    M_BI=[  cos(Euler(3))*cos(Euler(2))                                            sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
        cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
        cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
    
    k=1;
    for j=1:1:size((node_coords),1)
        s_inertial(k:k+5)=[M_BI'*node_coords(j,1:3)' [0 0 0]'];
        k=k+6;
    end
    s_inertial=(s_inertial'-I_Hat*M_BI'*node_coords(fixed_node,:)')';
    
    k=1;
    for j=1:1:size((node_coords),1)
        delta_inertial(k:k+2)=M_BI'*nodal_deflections((j-1)*6+1:(j-1)*6+3);
        delta_inertial(k+3:k+5)=M_BI'*nodal_deflections((j-1)*6+4:j*6);
        k=k+6;
    end
    b=s_inertial'+delta_inertial';
   
    b_Hat_Skew=zeros(n_dof,3);
    for bi=1:6:length(b)-1
        b_Hat_Skew(bi:bi+2,1:3)=[0 b(bi+2) -b(bi+1); -b(bi+2) 0 b(bi); b(bi+1) -b(bi) 0];
        % b_Hat_Skew(bi:bi+2,1:3)=[0 -b(bi+2) b(bi+1); b(bi+2) 0 -b(bi); -b(bi+1) b(bi) 0];
    end
    
    for bi=4:6:length(b)-1
        b_Hat_Skew(bi:bi+2,1:3)=[1 b(bi+2)*0 -b(bi+1)*0; -b(bi+2)*0 1 b(bi)*0; b(bi+1)*0 -b(bi)*0 1];
       % b_Hat_Skew(bi:bi+2,1:3)=[1 -b(bi+2) b(bi+1); b(bi+2) 1 -b(bi); -b(bi+1) b(bi) 1];
    end

    omega_Hat_Skew=zeros(n_dof,n_dof);
    for ai=1:6:n_dof
        omega_Hat_Skew(ai:ai+2,ai:ai+2)=[0 -pqr(3) pqr(2); pqr(3) 0 -pqr(1); -pqr(2) pqr(1) 0];
      % omega_Hat_Skew(ai:ai+2,ai:ai+2)=[0 pqr(3) -pqr(2); -pqr(3) 0 pqr(1); pqr(2) -pqr(1) 0];
    end
    for ai=4:6:n_dof
        omega_Hat_Skew(ai:ai+2,ai:ai+2)=1/2*[0 -pqr(3) pqr(2); pqr(3) 0 -pqr(1); -pqr(2) pqr(1) 0];
      % omega_Hat_Skew(ai:ai+2,ai:ai+2)=[0 pqr(3) -pqr(2); -pqr(3) 0 pqr(1); pqr(2) -pqr(1) 0];
    end

    A_Bar=zeros(n_dof,n_dof);
    for ai=1:3:n_dof
        A_Bar(ai:ai+2,ai:ai+2)=M_BI';
    end
end
