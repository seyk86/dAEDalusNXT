%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%clear all
close all

%% generate required paths for program execution
addpath(genpath('../stdlib'));
addpath(genpath('../geometry'));
addpath(genpath('../input'));
addpath(genpath('../aerodynamics'));
addpath(genpath('../structures'));
addpath(genpath('../aircraft'));
addpath(genpath('../aeroelastics'));
addpath(genpath('../criticaldesign'));
addpath(genpath('../xml_toolbox'));

structure_solver_settings=class_wingstructure_solver_settings;
structure_solver_settings.gravity=0;
aero_solver_settings=class_aero_solver_settings;
aeroelastic_solver_settings=class_aeroelastic_solver_settings;

fuselage_structure_solver_settings=class_fuselagestructure_solver_settings; 
fuselage_structure_solver_settings.nonlinear=0;

%% Define Weights
%clear all
Uinf=100;
alpha=10;
beta=0;
rho_air=1.225;
Ma=0;

state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);

%aircraft=class_aircraft('FLEX_OP.xml',1);
aircraft=class_aircraft('DISS_CWING.xml',1);
%aircraft=class_aircraft('A320dAEDalus_fuselage.xml',1);

result_path=['results/' aircraft.name '/'];

%set  structural grid size
aircraft.grid_settings.dy_max_struct_grid=aircraft.reference.b_ref/40;
aircraft.grid_settings.aerodynamic_fuselage=0; %deactivate aerodynamic fuselage model
% set aerodynamic panel sizes
aircraft.grid_settings.x_max_grid_size=aircraft.reference.c_ref/8;
aircraft.grid_settings.y_max_grid_size=aircraft.reference.b_ref/40;

% compute loftline (compute first)
aircraft=aircraft.compute_shell_coords();
aircraft=aircraft.compute_grid();


% addtional fuselage aero grid
%aircraft=aircraft.compute_fuselage_grid();

%aircraft.plot_grid_vol();

aircraft.write_tecplot_volgrid([result_path aircraft.name]);

[aircraft,aircraft_structure]=create_structural_model(aircraft);
%aircraft_structure.coupling_condition(4).node_idx(1)=19;
aircraft_structure=aircraft_structure.f_set_solver_settings(structure_solver_settings);

aircraft_structure.write_tecplot([result_path aircraft.name '_undeformed']);

% define reference state
ref_state=critical_ref_state(aircraft,0.4,0); % variables:  Mach and altitude equivalent to max dynamic pressure(knee of placard diagram = max dyn pressure)ref_state=critical_ref_state(aircraft,0.76,7620);

aircraft.write_tecplot_wingbox([result_path 'wingbox_geometry']);
%% define critical loadcases
% 2.5 g maneuver case
critical_state(1)=critical_g_maneuver_state(ref_state,2.5);
critical_state(1).loadcase_index=1;
% critical_state(2)=critical_g_maneuver_state(ref_state,-1.0);
% critical_state(2).loadcase_index=2;
critical_state(2)=critical_sideslip_maneuver_state(ref_state,10);
critical_state(2).loadcase_index=3;
% critical_state(4)=critical_sideslip_maneuver_state(ref_state,18);
% critical_state(3).loadcase_index=3;
 critical_state(3)=critical_aileron_roll_load_state(ref_state,2.5*2/3,'aileron',-25);
 critical_state(3).loadcase_index=3;
% critical_state(3)=critical_aileron_roll_load_state(ref_state,2.5*2/3,'aileron',25);
% critical_state(3).loadcase_index=3;

%critical_state(1).V=[critical_state(1).VD 0 0];


% manually set control surface defelctions for tests
 elev=0;
 aircraft=aircraft.f_set_control_surface('elevator',elev);
 aircraft=aircraft.f_set_control_surface('rudder',elev);
 aircraft=aircraft.f_set_control_surface('aileron_left',elev);
 aircraft=aircraft.f_set_control_surface('aileron_right',elev);
 % compute aerodynamic grid of aircraft wings
aircraft=aircraft.compute_grid();
 % compute aerodynamic grid of fuselage
%aircraft=aircraft.compute_fuselage_grid();
%  flight_state=critical_state(1);
%  flight_state.aerodynamic_state.p_ref=flight_state.aircraft_state.CG_ref;
%  wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,aircraft.reference);
%  wingaero.Uinf=204.18*[cosd(1) 0 sind(1)];
 %wingaero=wingaero.f_solve_std();

% aircraft.trim_surfaces{1}='elevon_left';
% aircraft.trim_surfaces{2}='elevon_right';

%aircraft=aircraft.compute_grid();
%aircraft=aircraft.compute_fuselage_grid();
%% INITIAL SIZING WITH 2.5G WITH RIGID LOADS
% force transformation matrix between aero grid and structure grid
aircraft=aircraft.compute_force_interpolation_matrix(aircraft_structure);
close all
% trim aircraft in critical_state(1)
[aircraft,trimmed_state,wingaero_trim]=trim_aircraft(aircraft,critical_state(1));
% 
%aircraft=aircraft.compute_forces(wingaero.F_body,wingaero.panels,wingaero.grid);
% transforms aeroloads to structure
aircraft=aircraft.compute_beam_forces(wingaero_trim.F_body,aircraft_structure);
for i=1:length(aircraft_structure.beam)
    if  isa(aircraft_structure.beam(i),'class_wing')
        aircraft_structure.beam(i)=aircraft_structure.beam(i).f_set_aeroloads(aircraft.wings(i));
        
    end
end
%run structures module
aircraft_structure.settings.engine=1; %boolean flag to consider engine inertia relief or not
aircraft_structure.settings.fuel_mass=1; %boolean flag to set fuel tank full or empty
% a320structure.settings.landing_gear=1;
% a320structure.plot_geometry

%aircraft_structure=aircraft_structure.f_set_state(critical_state(1));

% initial structural sizing with clamped aircraft (restrained at point defined in xml file)
aircraft_structure=aircraft_structure.f_load_based_self_design(aircraft.weights);
def=aircraft_structure.f_get_deflections;
aircraft=aircraft.compute_deflected_grid(def);
% structural sizing free flying
aircraft_structure=aircraft_structure.f_load_based_self_design_unrestrained(aircraft.weights);


%% ITERATIVE AEROELASTIC SIZING
% aircraft_structure.plot_externalforces(0,0,0,0)

total_mass=0;
for i=1:length(aircraft_structure.beam)
    for j=1:length(aircraft_structure.beam(i).beamelement)
        total_mass=total_mass+aircraft_structure.beam(i).beamelement(j).le*aircraft_structure.beam(i).beamelement(j).m;
        if isa(aircraft_structure.beam(i),'class_wing')
            total_mass=total_mass+aircraft_structure.beam(i).beamelement(j).el_m_fuel;
        end
    end
end

mean_axis_origin=[sum(aircraft_structure.Mglob_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords_full(:,1)))/total_mass;
                  sum(aircraft_structure.Mglob_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords_full(:,2)))/total_mass;
                  sum(aircraft_structure.Mglob_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords_full(:,3)))/total_mass];
             
critical_state(1).aircraft_state.CG_ref=mean_axis_origin';              
              
aeroelastic_solver_settings.CG_accelerations=1; %inertial relief on/off
aeroelastic_solver_settings.convergence_tol=0.5; % convergence tolerance for aeroelastic loop (deflection changes squared and normed
aeroelastic_solver_settings.unrestrained=1;
 [aircraft,aircraft_structure,wingaero]= structural_sizing_loop(aircraft,aircraft_structure,critical_state(1),aircraft.weights,aeroelastic_solver_settings);
% [aircraft,aircraft_structure,wingaero]= critical_case_layout(aircraft,aircraft_structure,critical_state([1 3]),aircraft.weights,aeroelastic_solver_settings);
%[aircraft,aircraft_structure]  = structural_sizing_loop_diss(aircraft,aircraft_structure,critical_state,aircraft.weights,aeroelastic_solver_settings);


aircraft.grid_settings.aerodynamic_fuselage=0;

aircraft_structure.write_tecplot([result_path aircraft.name '_deformed'])
aircraft_structure.write_tecplot_modes(result_path,20,30);
aircraft.write_modes_in_tecplot(aircraft_structure,20,40)
wingaero.write_tecplot([result_path aircraft.name '_aero_deformed'],1)
wingaero_trim.write_tecplot([result_path aircraft.name '_aero_undeformed'],1)
aircraft_structure.write_tecplot_mass([result_path aircraft.name '_mass_mtow']);
aircraft_structure.write_fuel_tecplot([result_path aircraft.name '_fuel']);
% 
% owe=1;
% aircraft_structure_owe=aircraft_structure.f_copy();
% for i=1:length(aircraft_structure_owe.beam)
%     if isa(aircraft_structure_owe.beam(i),'class_wing')
%         for j=1:length(aircraft_structure_owe.beam(i).beamelement)
%             if owe==1
%                 aircraft_structure_owe.beam(i).beamelement(j).is_fueled=0;
%             else
%                 aircraft_structure_owe.beam(i).beamelement(j).el_m_fuel=aircraft_structure_owe.beam(i).beamelement(j).el_fuel_vol*aircraft_structure_owe.beam(i).beamelement(j).is_fueled*aircraft_structure_owe.beam(i).fuel_density;
%             end
%         end
%     elseif isa(aircraft_structure_owe.beam(i),'class_fuselage')
% %         for j=1:length( aircraft_structure.beam(i).beamelement)  
% %             if owe==1
% %               %  aircraft_structure.beam(i).beamelement(j).el_m_s=0;
% %                 aircraft.weights.FuselageNonStructuralEstimate=aircraft.weights.FuselageNonStructuralEstimate*0;
% %             else
% %               %  aircraft_structure.beam(i).beamelement(j).el_m_s=0;
% %               aircraft.weights.FuselageNonStructuralEstimate=FuselageNonStructuralEstimate;
% %             end
% %         end
%     end
%     aircraft_structure_owe=aircraft_structure_owe.f_calc_mass(aircraft.weights);
% end
% 
% for i=1:length(aircraft_structure_owe.beam)
%     aircraft_structure_owe.beam(i).update_M=1;
% end
% 
% aircraft_structure_owe=aircraft_structure_owe.f_assemble_free(1,0);
% aircraft_structure_owe=aircraft_structure_owe.f_calc_mass(aircraft.weights);
% aircraft_structure_owe.write_tecplot_mass([result_path aircraft.name '_mass_owe']);
% aircraft_structure_owe.write_structure_tecplot([result_path aircraft.name '_structure']);

% aircraft.control.control_allocation_matrix=[0     1      0   0; ...
%                                             0     1       0   0; ...
%                                             -1    0       0   -1; ...
%                                             1     0       0   1; ...
%                                             0     0       1   0;];
%                                 
% %aircraft.control.trim_surfaces={'elevon_left','elevon_right'};
% aircraft.control.trim_surface_idx=[3,4];
% aircraft.control.lon_ctrl_idx=[1,2,5,6]; 
% aircraft.control.control_allocation_matrix_lon=ControlAllocationMatrix(aircraft.control.lon_ctrl_idx(3:end)-2,:);
% aircraft.control.lon_names={'Thr_L','Thr_R','Elevon_Left','Elevon_Right'};
% aircraft.control.lat_ctrl_idx=[3,4,7];
% aircraft.control.control_allocation_matrix_lat=ControlAllocationMatrix(aircraft.control.lat_ctrl_idx(1:end)-2,:);
% aircraft.control.lat_names={'Delta_AilLeft','Delta_AilRight','Delta_Rudder'};
