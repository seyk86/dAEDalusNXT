%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function aircraft_structure=f_compute_loads_rigid(aircraft,aircraft_structure,state,aeroelastic_solver_settings)

for state_ctr=1:length(state)
    %% sets loadfactor to structure
    aircraft_structure=aircraft_structure.f_set_state(state(state_ctr));
    %% compute rigid aeroloads
    [aircraft,~,wingaero]=trim_aircraft(aircraft,state(state_ctr));
    % transform forces from aero to structural grid
    aircraft=aircraft.compute_beam_forces(wingaero.F_body,aircraft_structure);
    
    if aeroelastic_solver_settings.CG_accelerations==1
        aircraft=aircraft.compute_acceleration(wingaero,aircraft_structure,state(state_ctr));
        %next line is wrong
        aircraft.acc=[0 0 0 0 0 0];
        aircraft_structure=aircraft_structure.f_set_acceleration(-aircraft.acc,wingaero.reference.p_ref);
    end
    
    for i=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(i),'class_wing')
            aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
        end
    end
    
    nonlin=0;
    loadstep=0;
    
    aircraft_structure=aircraft_structure.f_assemble(loadstep,nonlin);
    
    loadcase(state_ctr,:)=aircraft_structure.Fglob;
    
end

%% condense loads
for fvec_ctr=1:length(loadcase(1,:))
    sizing_loads(fvec_ctr)=max(abs(loadcase(:,fvec_ctr)));
end

aircraft_structure.Fglob=sizing_loads;

for i=1:length(aircraft_structure.beam)
    aircraft_structure.beam(i).update_Q=0;
end

end
