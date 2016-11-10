%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus aeroelasticity
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aircraft,aircraft_structure,wingaero,state ] = flight_state_loop(aircraft,aircraft_structure,state,aeroelastic_solver_settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('     performing aeroelastic flight state loop')
%% set state
aircraft_structure=aircraft_structure.f_set_state(state);

total_mass = aircraft_structure.f_compute_totalMass;

%% calculate initial state
%run aerodynamics module
[aircraft,state,wingaero]=trim_aircraft(aircraft,state);
%pass on aerodynamic forces to structures module
aircraft=aircraft.compute_beam_forces(wingaero.F_body,wingaero.M_body,aircraft_structure);

if aeroelastic_solver_settings.CG_accelerations==1
    aircraft=aircraft.compute_acceleration(wingaero,aircraft_structure,state);
    %aircraft.acc=[0 0 12.26 0 0 0];
    % aircraft.acc=[0 0 1.151780E+01 0 0 0];
    aircraft_structure=aircraft_structure.f_set_acceleration(-aircraft.acc,wingaero.reference.p_ref);
end

for i=1:length(aircraft_structure.beam)
    if  isa(aircraft_structure.beam(i),'class_wing')
        aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
    end
end
%run structures module
for i=1:length(aircraft_structure.beam)
    aircraft_structure.beam(i).update_M=1;
    aircraft_structure.beam(i).update_K=1;
    aircraft_structure.beam(i).update_Q=1;
end

if aeroelastic_solver_settings.unrestrained
    aircraft_structure=aircraft_structure.f_solve_unrestrained();
else
    aircraft_structure=aircraft_structure.f_solve();
end

%% iteration to reach static aeroelastic equilibrium
converged=0;
j=1;
df=1;
prev_error=1E12;
error=100;
damping=0;
def=0;
prev_def=aircraft_structure.nodal_deflections*0;
prv_def=0;
prv_prv_def=0;
prv_prv_prv_def=0;
%ds_lim=1;
disp('     entering aeroelastic loop');


while~converged
    
    j=j+1;
    % update aerodynamic geometry with structural deflections
    if damping==1
        disp(['     CONVERGENCE PROBLEMS: load stepping enabled:']);
        if j<=4
            prv_prv_prv_def=prv_def;
        else
            prv_prv_prv_def=prv_prv_def;
        end
        prv_prv_def=prv_def;
        prv_def=def;
        def=aircraft_structure.f_get_deflections;
        for i=1:length(aircraft_structure.beam)
            if  isa(aircraft_structure.beam(i),'class_wing')
                def(i).def=(prv_prv_prv_def(i).def+prv_prv_def(i).def+2*(1-error/100)*prv_def(i).def+2*(error/100)*def(i).def)/4;
            end
        end
        aircraft=aircraft.compute_deflected_grid(def);
        [aircraft,state,wingaero]=trim_aircraft(aircraft,state,aircraft_structure,'elastic',-5,def);
    else
        prv_prv_prv_def=prv_prv_def;
        prv_prv_def=prv_def;
        prv_def=def;
        def=aircraft_structure.f_get_deflections;
        aircraft=aircraft.compute_deflected_grid(def);
        [aircraft,state,wingaero]=trim_aircraft(aircraft,state,aircraft_structure,'elastic');
    end
    
    
    aircraft=aircraft.compute_beam_forces(wingaero.F_body,wingaero.M_body,aircraft_structure);
    
    if aeroelastic_solver_settings.CG_accelerations==1
        aircraft=aircraft.compute_acceleration(wingaero,aircraft_structure,state);
        %   aircraft.acc=[0 0 12.26 0 0 0];
        %          aircraft.acc=[0 0 1.151780E+01 0 0 0];
        aircraft_structure=aircraft_structure.f_set_acceleration(-aircraft.acc,wingaero.reference.p_ref);
    end
    
    for i=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(i),'class_wing')
            aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
        end
    end
    % solve beam to get new deflections
    for i=1:length(aircraft_structure.beam)
        aircraft_structure.beam(i).update_M=1;
        aircraft_structure.beam(i).update_K=1;
        aircraft_structure.beam(i).update_Q=1;
    end
    if aeroelastic_solver_settings.unrestrained
        aircraft_structure=aircraft_structure.f_solve_unrestrained();
    else
        aircraft_structure=aircraft_structure.f_solve();
    end
    
    if aeroelastic_solver_settings.unrestrained
		mean_axis_origin = aircraft_structure.f_compute_CG;
    	state.aircraft_state.CG_ref=mean_axis_origin';
	end
    
    prev_error=error;
    
    error=sum(sqrt((aircraft_structure.nodal_deflections-prev_def).^2))/sum(sqrt((aircraft_structure.nodal_deflections.^2)))*100;
    prev_def=aircraft_structure.nodal_deflections;
    %error=100*abs((aircraft_structure.beam(1).nodal_deflections(end-3)-prev_def))/aircraft_structure.beam(1).wing_frontview_length;
    if j>2
        if error>prev_error
            damping=1;
        end
    end
    
    disp(['     aeroelastic loop iteration: ' num2str(j) '  error: ' num2str(error) '%']);
    
    if((error<=aeroelastic_solver_settings.convergence_tol)&&(error>=0))
        if (j<=aeroelastic_solver_settings.max_it) && (df==1)
            converged=1;
        elseif df==length(0.0:0.1:1)+1
            converged=1;
        end
    end
    
    if j == aeroelastic_solver_settings.max_it
        disp('not converged, possible wing divergence');
        return
    end
    wingaero.plot_cp;
end

def=aircraft_structure.f_get_deflections;
aircraft=aircraft.compute_deflected_grid(def);
wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state.aerodynamic_state,aircraft.reference);
wingaero=wingaero.f_solve_std();
%      figure;
%      wingaero.plot_cp;
end

