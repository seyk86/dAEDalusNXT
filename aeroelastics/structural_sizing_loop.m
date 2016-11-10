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

function  [aircraft,aircraft_structure,wingaero ]  = structural_sizing_loop(aircraft,aircraft_structure,state,weights,aeroelastic_solver_settings,varargin)

if nargin==6
    overwrite=varargin{1};
else
    overwrite=1;
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
disp('     performing aeroelastic flight state loop')
%% set state
%set design speed to design dive speed
state.V=[state.VD 0 0];

%aircraft=aircraft.set_state(state);
aircraft_structure=aircraft_structure.f_set_state(state);
%% calculate initial state
%run aerodynamics module
[aircraft,state,wingaero]=trim_aircraft(aircraft,state);
%pass on aerodynamic forces to structures module
aircraft=aircraft.compute_beam_forces(wingaero.F_body,wingaero.F_body*0,aircraft_structure);

if aeroelastic_solver_settings.CG_accelerations==1
    aircraft=aircraft.compute_acceleration(wingaero,aircraft_structure,state);
    aircraft_structure=aircraft_structure.f_set_acceleration(-aircraft.acc,wingaero.reference.p_ref);
end

for i=1:length(aircraft_structure.beam)
    if  isa(aircraft_structure.beam(i),'class_wing')
        aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
        
    end
end
%run structures module
if aeroelastic_solver_settings.unrestrained
    aircraft_structure=aircraft_structure.f_load_based_self_design_unrestrained(weights,overwrite);
else
   aircraft_structure=aircraft_structure.f_load_based_self_design(weights,overwrite); 
end

%% iteration to reach static aeroelastic equilibrium
converged=0;
j=1;
df=1;

%ds_lim=1;
disp('     entering aeroelastic loop');


prev_def=aircraft_structure.nodal_deflections*0;
prev_error=0;
error=100;
damping=0;
prv_def=0;
prv_prv_def=0;
prv_prv_prv_def=0;
def=0;
n=2;

step=0:0.05:1;
n_step=length(step);
step_ctr=0;

while~converged
    
    j=j+1;
    figure;
    wingaero.plot_cp;
    
    % update aerodynamic geometry with structural deflections
    
    if damping==1
        disp(['     CONVERGENCE PROBLEMS: load stepping enabled:']);
        step_ctr=step_ctr+1;
%         if j<=4
%             prv_prv_prv_def=prv_def;
%         else
%             prv_prv_prv_def=prv_prv_def;
%         end
%         prv_prv_def=prv_def;
         prv_def=def;

        def=aircraft_structure.f_get_deflections;
        for i=1:length(aircraft_structure.beam)
            if  isa(aircraft_structure.beam(i),'class_wing')
              %  def(i).def=(prv_prv_prv_def(i).def+prv_prv_def(i).def+2*(1-error/100)*prv_def(i).def+2*(error/100)*def(i).def)/4;
              def(i).def=def(i).def*step(step_ctr);%(prv_prv_prv_def(i).def+prv_prv_def(i).def+prv_def(i).def+def(i).def)/4;
            end
        end
        
        if(step_ctr==n_step)
            damping=0;
            step_ctr=0;
        end
        aircraft=aircraft.compute_deflected_grid(def);
        [aircraft,state,wingaero]=trim_aircraft(aircraft,state,aircraft_structure,'elastic',1,def);
    else
        prv_prv_prv_def=prv_prv_def;
        prv_prv_def=prv_def;
        prv_def=def;
        def=aircraft_structure.f_get_deflections;
        aircraft=aircraft.compute_deflected_grid(def);
        [aircraft,state,wingaero]=trim_aircraft(aircraft,state,aircraft_structure,'elastic');
    end
    
%     figure
%     hold on
%     plot3(aircraft.wings(1).wingbox_c4(1,:)+def(1).def(end/2-2:6:end)',aircraft.wings(1).wingbox_c4(2,:)+def(1).def(end/2-1:6:end)',aircraft.wings(1).wingbox_c4(3,:)+def(1).def(end/2:6:end)','-xb')
%     plot3(aircraft.wings(1).wingbox_c4(1,:),aircraft.wings(1).wingbox_c4(2,:),aircraft.wings(1).wingbox_c4(3,:),'-xr')
%      plot3(aircraft.wings(1).wingbox_c4(1,:)+def(1).def(end/2-2:-6:1)',-aircraft.wings(1).wingbox_c4(2,:)+def(1).def(end/2-1:-6:1)',aircraft.wings(1).wingbox_c4(3,:)+def(1).def(end/2:-6:1)','-xb')
%      plot3(aircraft.wings(1).wingbox_c4(1,:),-aircraft.wings(1).wingbox_c4(2,:),aircraft.wings(1).wingbox_c4(3,:),'-xr')
%      aircraft.plot_grid_deflected
     
     aircraft=aircraft.compute_beam_forces(wingaero.F_aero,wingaero.F_aero*0,aircraft_structure);
    
    if aeroelastic_solver_settings.CG_accelerations==1
        aircraft=aircraft.compute_acceleration(wingaero,aircraft_structure,state);
        aircraft_structure=aircraft_structure.f_set_acceleration(-aircraft.acc,wingaero.reference.p_ref);
    end
    
    for i=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(i),'class_wing')
            aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
        end
    end

    % solve beam to get new deflections
    if aeroelastic_solver_settings.unrestrained
        for i=1:length(aircraft_structure.beam)
            aircraft_structure.beam(i).update_M=1;
            aircraft_structure.beam(i).update_K=1;
            aircraft_structure.beam(i).update_Q=1;
        end
        aircraft_structure=aircraft_structure.f_load_based_self_design_unrestrained(weights,overwrite);
    else
        aircraft_structure=aircraft_structure.f_load_based_self_design(weights,overwrite);
    end
    
    total_mass = aircraft_structure.f_compute_totalMass;
    
    mean_axis_origin = aircraft_structure.f_compute_CG;
             
    state.aircraft_state.CG_ref=mean_axis_origin';
    
    if aeroelastic_solver_settings.unrestrained
        state.aircraft_state.weight=total_mass;
    end
        
    prev_error=error;

    error=sum(sqrt((aircraft_structure.nodal_deflections-prev_def).^2))/sum(sqrt((aircraft_structure.nodal_deflections.^2)))*100;
    prev_def=aircraft_structure.nodal_deflections;
    
    if j>2
        if error>prev_error
            damping=1;
        end
    end
    
    disp(['     aeroelastic loop iteration: ' num2str(j) '  error: ' num2str(error) '%']);
    
    if((error<=aeroelastic_solver_settings.convergence_tol)&&(error>=0))
        if (j<=aeroelastic_solver_settings.max_it) && (step_ctr==0)
            converged=1;
        elseif step_ctr==n_step
            converged=1;
        end
    end
    
    if j == aeroelastic_solver_settings.max_it
        disp('not converged, possible wing divergence');
        return
    end
end

