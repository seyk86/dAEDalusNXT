%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [trimstate,trimreport,mean_axis_origin]=trim_free_flying_hybrid_simulink(aircraft,aircraft_structure,Aircraft,AeroelasticSSM,trim_state)

%% generate required paths for program execution
aircraft.weights.W= sum(aircraft_structure.Mff(1:6:end));

mean_axis_origin=-aircraft_structure.f_compute_CG;
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections*0,-mean_axis_origin,[0 0 0]');
[M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,aircraft_structure.modeshapes,AeroelasticSSM.nE+6,zeros(AeroelasticSSM.nE+6,1));
 
Trim_State=Flight_State_SIM;
Trim_State.AircraftState.CG=mean_axis_origin';

[~,~,~,rho_itr]=atmosisa(0);
% height=0;
% while (trim_state.aerodynamic_state.rho_air-rho_itr)<0.0001;
%     [~,~,~,rho_itr]=atmosisa(height); 
%     height=height+0.5;
% end
height=trim_state.h;

Trim_State.Xe=[0 0 -height];
Trim_State.AircraftState.m=M_tot_mean(1,1);
Trim_State.AircraftState.Ixyz=Aircraft.Inertia;
Trim_State.VCAS=trim_state.aerodynamic_state.V_A;
Trim_State.Alpha=trim_state.aerodynamic_state.alpha;
Trim_State.Euler=[0 trim_state.aerodynamic_state.alpha*pi/180 0];
Trim_State.rho=trim_state.aerodynamic_state.rho_air;
Trim_State.AircraftState.Trim_Deflection=0;
Trim_State.AircraftState.DeltaThr=0;
Trim_State.AircraftState.Conf=0;
Trim_State.AeroelasticState=0;
Trim_State.x_stat_ms=AeroelasticSSM.msASee^-1*1/2*trim_state.aerodynamic_state.rho_air*trim_state.aerodynamic_state.V_A^2*aircraft.reference.S_ref*[zeros(AeroelasticSSM.nE,1);AeroelasticSSM.Q_0(7:6+AeroelasticSSM.nE,1);zeros(size(AeroelasticSSM.msASre(:,:),1)-2*AeroelasticSSM.nE,1)];
Trim_State.AeroelasticityActive=1;
Environment.ACTIVATE_GUST=0;
Environment.ATMOSPHERE_g=9.81;
SimMode=0;
% control_in.time=ctrl_input.time;
% control_in.signals.values=zeros(length(control_in.time),ctrl_input.signals.dimensions+4+7);
% control_in.signals.dimensions=ctrl_input.signals.dimensions+4+7;
% for n_ctrl=1:length(control_surface_states)
%     control_in.signals.values(:,2+n_ctrl)=ctrl_input.signals.values(:,n_ctrl)+trimmed_state.aircraft_state.control_deflections{n_ctrl};
% end 


assignin('base','Environment',Environment);
%assignin('base','control_in',control_in);
assignin('base','SimMode',SimMode);
assignin('base','Aircraft',Aircraft);



[trimstate,trimreport]=f_trim_aircraft_model(Aircraft,Trim_State,'horizontal',trim_state.aerodynamic_state.V_A,0);
%     
%     [Init_Flight_State,Aircraft] = f_flightstate_init(Aircraft,trimstate,Trim_State.AircraftState.Conf,0,Trim_State.AircraftState.m,CG_percentage,0,0,0,Trim_State.AeroelasticityActive);
%     control_in.time=[0 5];
%     control_in.signals.values(1:2,1:2)=trimstate.Inputs(1).u(1);
%     control_in.signals.values(1:2,3)=trimstate.Inputs(1).u(3);
%     control_in.signals.values(1:2,4)=trimstate.Inputs(1).u(4);
%     CutDynamics=0;
%     SimMode=0;
%  [Init_Flight_State] = f_flightstate_init(Aircraft,trimstate(i),conf,0,60000,trimstate_cgs(i),0,0,0);



