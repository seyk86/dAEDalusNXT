%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,aircraft_structure,steady_trim_aero,state_VLM,state_UVLM] = trim_elastic_aircraft(aircraft,aircraft_structure,trim_state,aeroelastic_solver_settings,varargin)

if nargin==5
    uvlm_flag=varargin{1};
else
    uvlm_flag=1;
end

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
aircraft_structure=aircraft_structure.f_postprocess();

total_mass = aircraft_structure.f_compute_totalMass;

mean_axis_origin = aircraft_structure.f_compute_CG;

[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,[0 0 0]');
[M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff_lumped,aircraft_structure.Kff,aircraft_structure.modeshapes,length(aircraft_structure.Kff),zeros(length(aircraft_structure.Kff),1));
    
%aeroelastic_solver_settings=class_aeroelastic_solver_settings;

%aircraft.reference.p_ref=mean_axis_origin';

% set current weight
aircraft.weights.W=aircraft_structure.f_compute_totalMass;
trim_state.aircraft_state.weight=aircraft.weights.W;
trim_state.aircraft_state.CG_ref=mean_axis_origin';
trim_state.aircraft_state.I_xyz=M_tot_mean(4:6,4:6);
[aircraft,aircraft_structure,steady_trim_aero,state_VLM] = flight_state_loop(aircraft,aircraft_structure,trim_state,aeroelastic_solver_settings);

if uvlm_flag==1
    state_UVLM=0;
    [aircraft,state_UVLM,wingaero] = trim_aircraft_unsteady(aircraft,trim_state,aircraft_structure,'elastic');
else
    state_UVLM=[];
end

Mass=M_tot_mean(1,1);
Inertia=M_tot_mean(4:6,4:6);

mkdir(aircraft.name)
save([aircraft.name '/MassInertia'],'Mass','Mass','Inertia','Inertia');


end

