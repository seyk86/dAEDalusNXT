%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,aircraft_structure,wingaero,trimmed_AESSM_state,ASE] = trim_elastic_aircraft_AESSM(aircraft,aircraft_structure,AeroelasticSSM,trim_state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
aircraft_structure=aircraft_structure.f_postprocess();

total_mass=0;

for i=1:length(aircraft_structure.beam)
    for j=1:length(aircraft_structure.beam(i).beamelement)
           total_mass=total_mass+aircraft_structure.beam(i).beamelement(j).le*aircraft_structure.beam(i).beamelement(j).m;
           if isa(aircraft_structure.beam(i),'class_wing')
               total_mass=total_mass+aircraft_structure.beam(i).beamelement(j).el_m_fuel;
           end
    end
end

mean_axis_origin=[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/total_mass;
                      sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/total_mass;
                      sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/total_mass];
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,[0 0 0]');
[M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff_lumped,aircraft_structure.Kff,aircraft_structure.modeshapes,length(aircraft_structure.Kff),zeros(length(aircraft_structure.Kff),1));
    
aeroelastic_solver_settings=class_aeroelastic_solver_settings;

%aircraft.reference.p_ref=mean_axis_origin';

% set current weight
aircraft.weights.W=M_tot_mean(1,1);
trim_state.aircraft_state.weight=M_tot_mean(1,1);
trim_state.aircraft_state.CG_ref=mean_axis_origin';
trim_state.aircraft_state.I_xyz=M_tot_mean(4:6,4:6);

[aircraft,trimmed_AESSM_state,wingaero,ASE]=trim_aircraft_AESSM(aircraft,trim_state,AeroelasticSSM);
end

