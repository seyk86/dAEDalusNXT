%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function FixedAxis=simulate_constrained(aircraft,aircraft_structure,fixed_node,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

if ctrl_input.signals.dimensions~=length(aircraft.control_surfaces)
    disp('Control Input Signal Incorrect')
end

% set newmark beta scheme parameters
beta_newmark=1/4;
gamma_newmark=1/2;
% initialize flight state variables

FixedAxis.Alpha=zeros(length(t_vec),1);
FixedAxis.V=zeros(length(t_vec),1);
FixedAxis.qinf=zeros(length(t_vec),1);
FixedAxis.Ma=zeros(length(t_vec),1);
% for testing
gravity_on=1;

    aircraft_structure=aircraft_structure.f_solve();
% initialize state vectors to zero
x_now=zeros(length(aircraft_structure.Kff),1);
x_nxt=x_now;
xdot_now=x_now;
xdotdot_now=x_now;
xdotdot_nxt=x_now;

aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();

aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);

x_now_init= aircraft_structure.nodal_deflections;

x_now=x_now_init;

aeroelastic_solver_settings=class_aeroelastic_solver_settings;

UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=movie_on;
aircraft.reference.p_ref=trimmed_state.aerodynamic_state.p_ref;
aircraft_aero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
aircraft_aero=aircraft_aero.initialize_time_domain_solution(dt);


FixedAxis.Time=t_vec;
x_prv=x_now;
xdot_prv=xdot_now;
xdotdot_prv=xdotdot_now;

control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;
rho_air=trimmed_state.aerodynamic_state.rho_air;


    
    M_tot_mean = aircraft_structure.Mff;
    K_tot_mean = aircraft_structure.Kff;
        C_mean=K_tot_mean*0;
    C_mean(7:end,7:end)=K_tot_mean(7:end,7:end)*0.002+M_tot_mean(7:end,7:end)*0.002;
    
        gravity=[0 0 -9.81]'*gravity_on;
    for bb=1:length(aircraft_structure.beam)
        for be=1:length(aircraft_structure.beam(bb).beamelement)
            aircraft_structure.beam(bb).beamelement(be).ax=gravity(1);
            aircraft_structure.beam(bb).beamelement(be).ay=gravity(2);
            aircraft_structure.beam(bb).beamelement(be).az=gravity(3);
            aircraft_structure.beam(bb).update_Q=1;
        end
    end
    
% Defining fixed flight parameters for the constrained plane

beta=1/4;
gamma=1/2;
% figure()
% plot(aircraft_structure.nodal_deflections,'x')
% hold on
for i = 1:length(t_vec)
    
    
    aircraft_aero=aircraft_aero.solve_time_domain_aerodynamics(aircraft,[0 0 0 0 trimmed_state.aerodynamic_state.alpha 0],trimmed_state.aerodynamic_state.V_A,trimmed_state.aerodynamic_state.alpha,0,...
        [0 0 0] ,rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/constrained_flutter']);
    
    aircraft_aero.F_body(1,:)=-sind(trimmed_state.aerodynamic_state.alpha).*aircraft_aero.F_body(3,:);
    
    aircraft=aircraft.compute_beam_forces(aircraft_aero.F_body,aircraft_structure);
    for bb=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(bb),'class_wing')
            aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
        end
    end
   
    aircraft_structure=aircraft_structure.f_solve();
    
    F_tot_mean = aircraft_structure.Ftest;
    M_tot_mean = aircraft_structure.Mff;
    K_tot_mean = aircraft_structure.Kff;
    C_mean=K_tot_mean*0;
    C_mean(7:end,7:end)=K_tot_mean(7:end,7:end)*0.002+M_tot_mean(7:end,7:end)*0.002;
    
    
%     % ORIGINAL
%     x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
%         +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));    
    x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
            +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));
    xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
    xdot_nxt=xdot_now+dt*((1-gamma_newmark)*xdotdot_now+gamma_newmark*xdotdot_nxt);
    
    % FROM class_UVLM_solver.time_domain_flutter_analysis
%     x_nxt=linsolve(aircraft_structure.Mff+beta*dt^2*aircraft_structure.Kff,beta*dt^2*aircraft_structure.Ftest+aircraft_structure.Mff*x_now+dt*aircraft_structure.Mff*xdot_now+dt^2*aircraft_structure.Mff*(1/2-beta)*xdotdot_now);
%     xdotdot_nxt=1/(beta*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta)*xdotdot_now);         %Match with original
%     xdot_nxt=xdot_now+dt*((1-gamma)*xdotdot_now+gamma*xdotdot_nxt);                          %Match with original
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    aircraft_structure.nodal_deflections=x_nxt;
%     plot(aircraft_structure.nodal_deflections)
    aircraft_structure=aircraft_structure.f_postprocess();
    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    x_now=x_nxt;
    xdot_now=xdot_nxt;
    xdotdot_now=xdotdot_nxt; 
    
%     if max(x_nxt)>1
%         I=i;
%         
%         plot(x_nxt)
%         drawnow
%         pause
%     end
%     
end
