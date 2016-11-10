%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Q,Q_0] = compute_GAF_matrix(k,nmod,aircraft,aircraft_structure,state,name,UVLM_settings,exag,varargin)
%check if model is constrained
if aircraft_structure.modefrequencies(7)/aircraft_structure.modefrequencies(6)>100
    free_modes=1;
else
    free_modes=0;
end
% check if trim shape is given
if nargin==9
    shape=varargin{1};
else
    shape=0;
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if shape~=0
    aircraft_structure.nodal_deflections=shape;
    aircraft_structure=aircraft_structure.f_postprocess();
    def=aircraft_structure.f_get_deflections;
    aircraft=aircraft.compute_deflected_grid(def);
    %aircraft.grid=aircraft.grid_deflected;
end

if free_modes==1
    aircraft_structure_back=aircraft_structure;
    aircraft_structure.modeshapes=aircraft_structure.modeshapes(:,7:end);
end

n_cmod=length(aircraft.control_surfaces);
cmod=0;
n_gmod=0;
Q=zeros(6+nmod+n_gmod+n_cmod,6+nmod+n_gmod+n_cmod,length(k));
Q_0=zeros(6+nmod+n_gmod+n_cmod,6+nmod+n_gmod+n_cmod,length(k));
warning off
wingaero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);

for kk=1:length(k)
%     if free_modes
        amplitude=0.02*aircraft.reference.c_ref;
        wingaero=wingaero.solve_unsteady_foraft(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape);
        Q(1,1,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/amplitude;
        Q(2,1,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/amplitude;
        Q(3,1,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/amplitude;
        Q(4,1,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/amplitude;
        Q(5,1,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/amplitude;
        Q(6,1,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/amplitude;
        Q_0(1,1,kk)=wingaero.CX_complex(3);
        Q_0(2,1,kk)=wingaero.CY_complex(3);
        Q_0(3,1,kk)=wingaero.CZ_complex(3);
        Q_0(4,1,kk)=wingaero.CL_complex(3);
        Q_0(5,1,kk)=wingaero.CM_complex(3);
        Q_0(6,1,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,1,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/amplitude;
            Q_0(6+mod,1,kk)=wingaero.C_MODES_complex(mod,3);
        end
        amplitude=0.02*aircraft.reference.c_ref;
        wingaero=wingaero.solve_unsteady_sideheave(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape);
        Q(1,2,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/amplitude;
        Q(2,2,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/amplitude;
        Q(3,2,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/amplitude;
        Q(4,2,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/amplitude;
        Q(5,2,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/amplitude;
        Q(6,2,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/amplitude;
        Q_0(1,2,kk)=wingaero.CX_complex(3);
        Q_0(2,2,kk)=wingaero.CY_complex(3);
        Q_0(3,2,kk)=wingaero.CZ_complex(3);
        Q_0(4,2,kk)=wingaero.CL_complex(3);
        Q_0(5,2,kk)=wingaero.CM_complex(3);
        Q_0(6,2,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,2,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/amplitude;
            Q_0(6+mod,2,kk)=wingaero.C_MODES_complex(mod,3);
        end
        amplitude=0.02*aircraft.reference.c_ref;
        wingaero=wingaero.solve_unsteady_heave(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape);
        Q(1,3,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/amplitude;
        Q(2,3,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/amplitude;
        Q(3,3,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/amplitude;
        Q(4,3,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/amplitude;
        Q(5,3,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/amplitude;
        Q(6,3,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/amplitude;
        Q_0(1,3,kk)=wingaero.CX_complex(3);
        Q_0(2,3,kk)=wingaero.CY_complex(3);
        Q_0(3,3,kk)=wingaero.CZ_complex(3);
        Q_0(4,3,kk)=wingaero.CL_complex(3);
        Q_0(5,3,kk)=wingaero.CM_complex(3);
        Q_0(6,3,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,3,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/amplitude;
            Q_0(6+mod,3,kk)=wingaero.C_MODES_complex(mod,3);
        end
        amplitude=0.1;
        wingaero=wingaero.solve_unsteady_roll(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape);
        Q(1,4,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/(pi/180*amplitude);
        Q(2,4,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/(pi/180*amplitude);
        Q(3,4,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/(pi/180*amplitude);
        Q(4,4,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/(pi/180*amplitude);
        Q(5,4,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/(pi/180*amplitude);
        Q(6,4,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/(pi/180*amplitude);
        Q_0(1,4,kk)=wingaero.CX_complex(3);
        Q_0(2,4,kk)=wingaero.CY_complex(3);
        Q_0(3,4,kk)=wingaero.CZ_complex(3);
        Q_0(4,4,kk)=wingaero.CL_complex(3);
        Q_0(5,4,kk)=wingaero.CM_complex(3);
        Q_0(6,4,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,4,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/(pi/180*amplitude);
            Q_0(6+mod,4,kk)=wingaero.C_MODES_complex(mod,3);
        end
        amplitude=0.1;
        wingaero=wingaero.solve_unsteady_pitch(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape*0);
        Q(1,5,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/(pi/180*amplitude);
        Q(2,5,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/(pi/180*amplitude);
        Q(3,5,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/(pi/180*amplitude);
        Q(4,5,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/(pi/180*amplitude);
        Q(5,5,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/(pi/180*amplitude);
        Q(6,5,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/(pi/180*amplitude);
        Q_0(1,5,kk)=wingaero.CX_complex(3);
        Q_0(2,5,kk)=wingaero.CY_complex(3);
        Q_0(3,5,kk)=wingaero.CZ_complex(3);
        Q_0(4,5,kk)=wingaero.CL_complex(3);
        Q_0(5,5,kk)=wingaero.CM_complex(3);
        Q_0(6,5,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,5,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/(pi/180*amplitude);
            Q_0(6+mod,5,kk)=wingaero.C_MODES_complex(mod,3);
        end
        amplitude=0.1;
        wingaero=wingaero.solve_unsteady_yaw(amplitude,k(kk),nmod,aircraft,aircraft_structure,shape);
        Q(1,6,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/(pi/180*amplitude);
        Q(2,6,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/(pi/180*amplitude);
        Q(3,6,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/(pi/180*amplitude);
        Q(4,6,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/(pi/180*amplitude);
        Q(5,6,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/(pi/180*amplitude);
        Q(6,6,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/(pi/180*amplitude);
        Q_0(1,6,kk)=wingaero.CX_complex(3);
        Q_0(2,6,kk)=wingaero.CY_complex(3);
        Q_0(3,6,kk)=wingaero.CZ_complex(3);
        Q_0(4,6,kk)=wingaero.CL_complex(3);
        Q_0(5,6,kk)=wingaero.CM_complex(3);
        Q_0(6,6,kk)=wingaero.CN_complex(3);
        for mod=1:nmod
            Q(6+mod,6,kk)=(wingaero.C_MODES_complex(mod,1)+1i*wingaero.C_MODES_complex(mod,2))/(pi/180*amplitude);
            Q_0(6+mod,6,kk)=wingaero.C_MODES_complex(mod,3);
        end
%     else
%        Q(1:6+nmod,1:6,kk)=0;
%       Q_0(1:6+nmod,1:6,kk)=0;
%     end
    for mod=1:nmod
        amplitude=exag(mod);
        wingaero=wingaero.solve_unsteady_mode(amplitude,k(kk),nmod,aircraft,aircraft_structure,mod,shape);
        Q(1,6+mod,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/amplitude;
        Q(2,6+mod,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/amplitude;
        Q(3,6+mod,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/amplitude;
        Q(4,6+mod,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/amplitude;
        Q(5,6+mod,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/amplitude;
        Q(6,6+mod,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/amplitude;
        Q_0(1,6+mod,kk)=wingaero.CX_complex(3);
        Q_0(2,6+mod,kk)=wingaero.CY_complex(3);
        Q_0(3,6+mod,kk)=wingaero.CZ_complex(3);
        Q_0(4,6+mod,kk)=wingaero.CL_complex(3);
        Q_0(5,6+mod,kk)=wingaero.CM_complex(3);
        Q_0(6,6+mod,kk)=wingaero.CN_complex(3);
        for i_mod=1:nmod
            Q(6+i_mod,6+mod,kk)=(wingaero.C_MODES_complex(i_mod,1)+1i*wingaero.C_MODES_complex(i_mod,2))/amplitude;
            Q_0(6+i_mod,6+mod,kk)=wingaero.C_MODES_complex(i_mod,3);
        end
    end

    for cmod=1:n_cmod
        for set0=1:n_cmod
           aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{set0},aircraft.control_deflections{cmod});
        end
        amplitude=0.5;
        wingaero=wingaero.solve_unsteady_controlmode(amplitude,k(kk),nmod,aircraft,aircraft_structure,aircraft.control_surfaces{cmod},aircraft.control_deflections{cmod},shape);
        Q(1,6+nmod+cmod,kk)=(  wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/(pi/180*amplitude);
        Q(2,6+nmod+cmod,kk)=(  wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/(pi/180*amplitude);
        Q(3,6+nmod+cmod,kk)=(  wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/(pi/180*amplitude);
        Q(4,6+nmod+cmod,kk)=(  wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/(pi/180*amplitude);
        Q(5,6+nmod+cmod,kk)=(  wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/(pi/180*amplitude);
        Q(6,6+nmod+cmod,kk)=(  wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/(pi/180*amplitude);
        Q_0(1,6+nmod+cmod,kk)=wingaero.CX_complex(3);
        Q_0(2,6+nmod+cmod,kk)=wingaero.CY_complex(3);
        Q_0(3,6+nmod+cmod,kk)=wingaero.CZ_complex(3);
        Q_0(4,6+nmod+cmod,kk)=wingaero.CL_complex(3);
        Q_0(5,6+nmod+cmod,kk)=wingaero.CM_complex(3);
        Q_0(6,6+nmod+cmod,kk)=wingaero.CN_complex(3);
        for i_mod=1:nmod
            Q(6+i_mod,6+nmod+cmod,kk)=( wingaero.C_MODES_complex(i_mod,1)+1i*wingaero.C_MODES_complex(i_mod,2))/(pi/180*amplitude);
            Q_0(6+i_mod,6+nmod+cmod,kk)=wingaero.C_MODES_complex(i_mod,3);
        end
    end
   
%     for gmod=1:n_gmod
%         amplitude=5;
%         wingaero=wingaero.solve_unsteady_gust_mode(amplitude,k(kk),nmod,aircraft,fullstructure);
%         Q(1,6+nmod+n_cmod+gmod,kk)=( wingaero.CX_complex(1)+1i*wingaero.CX_complex(2))/amplitude;
%         Q(2,6+nmod+n_cmod+gmod,kk)=( wingaero.CY_complex(1)+1i*wingaero.CY_complex(2))/amplitude;
%         Q(2,6+nmod+n_cmod+gmod,kk)=( wingaero.CZ_complex(1)+1i*wingaero.CZ_complex(2))/amplitude;
%         Q(4,6+nmod+n_cmod+gmod,kk)=( wingaero.CL_complex(1)+1i*wingaero.CL_complex(2))/amplitude;
%         Q(3,6+nmod+n_cmod+gmod,kk)=( wingaero.CM_complex(1)+1i*wingaero.CM_complex(2))/amplitude;
%         Q(6,6+nmod+n_cmod+gmod,kk)=( wingaero.CN_complex(1)+1i*wingaero.CN_complex(2))/amplitude;
%         Q_0(1,6+mod+cmod+gmod,kk)=wingaero.CX_complex(3);
%         Q_0(2,6+mod+cmod+gmod,kk)=wingaero.CY_complex(3);
%         Q_0(3,6+mod+cmod+gmod,kk)=wingaero.CZ_complex(3);
%         Q_0(4,6+mod+cmod+gmod,kk)=wingaero.CL_complex(3);
%         Q_0(5,6+mod+cmod+gmod,kk)=wingaero.CM_complex(3);
%         Q_0(6,6+mod+cmod+gmod,kk)=wingaero.CN_complex(3);
%         for i_mod=1:nmod
%             Q(6+i_mod,6+nmod+n_cmod+gmod,kk)=( wingaero.C_MODES_complex(i_mod,1)+1i*wingaero.C_MODES_complex(i_mod,2))/amplitude;
%             Q_0(6+i_mod,6+nmod+n_cmod+gmod,kk)=wingaero.C_MODES_complex(i_mod,3);
%         end
%     end
end


mkdir(aircraft.name)
if free_modes==1
    aircraft_structure=aircraft_structure_back;
end

save([aircraft.name '/' name],'Q','Q_0','k','aircraft_structure','aircraft','-v7.3');



end

