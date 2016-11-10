%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%% dAEDalus time domain flight simulation
% before this script, at least generate_aircraft_model.m has to be executed

%% calculation of starting trimmed state
trim_state=critical_ref_state(aircraft,0.4,4000);
aeroelastic_solver_settings=class_aeroelastic_solver_settings;
[aircraft,aircraft_structure,steady_trim_aero,trimmed_state,trimmed_stateUVLM] = trim_elastic_aircraft(aircraft,aircraft_structure,trim_state,aeroelastic_solver_settings);

%% different control input signals might be simulated. some examples:
%                                   generate_input_signal(aircraft,cs_name,                          name,              type,   step_length,    tstart,     tend,   dt,     amplitude)
[ ctrl_input,t_vec,maneuver_name] = generate_input_signal(aircraft,{'aileron_left','aileron_right'},'aileron_doublet','doublet',  1,            1,          8,     0.005,   [5 5]);

%% FREE FLYING RIGID
simulate_free_flying_rigid
simulate_free_flying_rigid_unsteady

%% FREE FLYING FIXED AXIS Simulation
dt=t_vec(2)-t_vec(1);
fixed_node=1;
movie_on=1;

FixedAxisEL=simulate_free_flying_fixed_axis(aircraft,aircraft_structure,fixed_node,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
FixedAxisEL.Name='fixed-axis-el'; 

%% FREE FLYING FIXED AXIS Simulation QUASISTEDY
FixedAxisELQS=simulate_free_flying_fixed_axis_qs(aircraft,aircraft_structure,fixed_node,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
FixedAxisELQS.Name='fixed-axis-el-qs'; 

%% FREE FLYING MEAN AXIS Simulation
MeanAxisEL=simulate_free_flying_mean_axis(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisEL.Name='mean-axis-el'; 

%% FREE FLYING MEAN AXIS Simulation QUASISTEADY 
MeanAxisELQS=simulate_free_flying_mean_axis_qs(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELQS.Name='mean-axis-el-qs'; 

%% FREE FLYING MEAN AXIS Simulation Inertialy Decoupled
MeanAxisELID=simulate_free_flying_mean_axis_inertially_decoupled(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELID.Name='mean-axis-el-id';

%% FREE FLYING MEAN AXIS Simulation Modal
MeanAxisELMOD=simulate_free_flying_mean_axis_modal(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELMOD.Name='mean-axis-el-mod';

%% FREE FLYING MEAN AXIS Simulation Modal Inertialy Decoupled
MeanAxisELMODID=simulate_free_flying_mean_axis_modal_idcpl(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELMODID.Name='mean-axis-el-mod-id';

%% FREE FLYING MEAN AXIS Simulation Modal Quasisteady
MeanAxisELMODQS=simulate_free_flying_mean_axis_modal_qs(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELMODQS.Name='mean-axis-el-mod-qs';

%% FREE FLYING MEAN AXIS Simulation Modal Quasisteady Inertialy Decoupled
MeanAxisELMODIDQS=simulate_free_flying_mean_axis_modal_qs_id(aircraft,aircraft_structure,trimmed_stateUVLM,ctrl_input,t_vec,dt,maneuver_name,movie_on);
MeanAxisELMODIDQS.Name='mean-axis-el-mod-id-qs';

%% FREE FLYING MEAN AXIS Simulation Partitioned 
simulate_free_flying_mean_axis_partitioned

%% FREE FLYING MEAN AXIS Simulation Partitioned w Coefficients
simulate_free_flying_mean_axis_partitioned_w_coefficients

%% Hybrid simulations
% input AeroelasticSSM
%% evtl Roger
simulate_free_flying_hybrid

%% Roger + uvlm
simulate_free_flying_hybrid_unsteady

%% minimum state input AeroelasticSSM
simulate_free_flying_hybrid_ms

%% minimum state mit uvlm
simulate_free_flying_hybrid_ms_unsteady

%% executes simulink model with ctrl input as defined before. here simulink model is needed 
simulate_free_flying_hybrid_simulink	

%% to check what this is:
simulate_free_flying_stat





%% GUST RESPONSE
UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=0;
Uinf=trim_state;
alpha=0;
beta=0;
rho_air=1.225;
Ma=0.15;
state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);

aircraft_ow.grid_settings.wake=1;
aircraft_ow=aircraft_ow.compute_grid();

t_start=0.5;
wingaero_gx=class_UVLM_solver(aircraft_ow.name,aircraft_ow.grid,aircraft_ow.is_te,aircraft_ow.panels,state,aircraft_ow.grid_wake,aircraft_ow.panels_wake,aircraft_ow.reference,UVLM_settings);

H=0.5; % gust length
Uds=10; %gust speed
s1=0:(2*H/(20)):2*H;
t_d1=2*H/norm(wingaero_gx.Uinf);
wingaero_gx.settings.movie=0;
wingaero_gx=wingaero_gx.solve_unsteady_aeroelastic_gust_response(aircraft_ow,aircraft_structure_ow,Uds,H,t_start,6,aircraft_structure_ow.nodal_deflections*0);


