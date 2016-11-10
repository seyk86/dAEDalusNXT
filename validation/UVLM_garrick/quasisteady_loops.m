%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
clear all
%close all
% generate required paths for program execution
addpath(genpath('../../aerodynamics'));
addpath(genpath('../../geometry'));
addpath(genpath('../../aircraft'));
addpath(genpath('../../stdlib'));
addpath(genpath('../../input'));
M=16;
UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.wakelength_factor=0.9;
UVLM_settings.n_osc=8;
%minimum steps per period
UVLM_settings.spp=16;
aero_solver_settings=class_aero_solver_settings;

%% define aircraft geometry
%aircraft=class_aircraft('AGARD445.6Wing.xml',1);

%aircraft=class_aircraft('optimwing.xml',1)

aircraft=class_aircraft('AR40wing.xml',1);
% aircraft.plot
% aircraft.weights=weights;
aircraft.grid_settings.x_max_grid_size=1.05*aircraft.reference.c_ref/M;
aircraft.grid_settings.y_max_grid_size=20;
% aircraft.grid_settings.x_max_grid_size=0.18;
% aircraft.grid_settings.y_max_grid_size=0.25;
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();
  
Uinf=10;
alpha=0;
beta=0;
Ma=0.0;

rho_air=0.397;
state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);


k_UVLM=[0.001];
 hh=[10]; % h hat -> h hat = h*b (h @ palacios paper)
 ah=1*pi/180;
 aa=-0.5;
 wingaeros=[];
 Q_UVLM=zeros(3,3,length( k_UVLM));
 kk=1
 wingaero.t_vec=0:0.01:1000;
 wingaero.reference=aircraft.reference;
 generate_garrick_data;
 %% plot garrick data
 figure
  hold on
  plot(alpha_eff*180/pi,Cd_heave,'--r','LineWidth',1.5)

%% quasisteady aeff sim
 t_step=aircraft_reference/(Uinf*M) ;
 
 cdi=[];
 for j=1:length(t_step)
    qs_state=state;
    qs_state.alpha=max(alpha_eff)*180/pi;
    testaero_settings=UVLM_settings;
    testaero_settings.wakelength_factor=0.45;
    testaero=class_UVLM_solver(aircraft.name,aircraft.grid,aircraft.is_te,aircraft.panels,qs_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,testaero_settings);
    testaero=testaero.initialize_time_domain_solution(t_step);
    testaero.CX/2
 end

%% steady VLM
    wingaeroVLM=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,qs_state,aircraft.reference); 
    wingaeroVLM=wingaeroVLM.f_solve_std();
    wingaeroVLM.CX/2
