%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [FixedAxis,  R_tip, velocity_t]=simulate_constrained_flutter(aircraft,aircraft_structure,aero_state,t_vec,maneuver_name,movie_on, velIncr)
dt=t_vec(2)-t_vec(1);
mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');
modalAnalysis=1;
numberModes=30;
breakevery=50;
flutCrit=2e-5;
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
    if modalAnalysis
        aircraft_structure=aircraft_structure.f_solve_modes();
        modalBasis=aircraft_structure.modeshapes(:,1:numberModes);
    end
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();
% initialize state vectors to zero




aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
if modalAnalysis
    x_now=zeros(numberModes,1);
    x_nxt=x_now;
    xdot_nxt=x_now;
    xdot_now=x_now;
    xdotdot_now=x_now;
    xdotdot_nxt=x_now;
    init_cond=aircraft_structure.modeshapes^-1*aircraft_structure.nodal_deflections;
    x_now_init=init_cond(1:numberModes);
    x_now=x_now_init;
else
    x_now=zeros(length(aircraft_structure.Kff),1);
    x_nxt=x_now;
    xdot_nxt=x_now;
    xdot_now=x_now;
    xdotdot_now=x_now;
    xdotdot_nxt=x_now;
    x_now_init= aircraft_structure.nodal_deflections;
    x_now=x_now_init;
end

aeroelastic_solver_settings=class_aeroelastic_solver_settings;

UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=movie_on;
UVLM_settings.wakelength_factor=2;
aircraft.reference.p_ref=aero_state.p_ref;
aircraft_aero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,aero_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
aircraft_aero=aircraft_aero.initialize_time_domain_solution(dt);


FixedAxis.Time=t_vec;
x_prv=x_now;
xdot_prv=xdot_now;
xdotdot_prv=xdotdot_now;


rho_air=aero_state.rho_air;
C_mod=diag(0.02*(2*pi*aircraft_structure.modefrequencies).^2);
if modalAnalysis
    M_tot_mean = modalBasis'*aircraft_structure.Mff*modalBasis;
    K_tot_mean = modalBasis'*aircraft_structure.Kff*modalBasis;
    C_mean=C_mod(1:numberModes,1:numberModes);
else
    M_tot_mean = aircraft_structure.Mff;
    K_tot_mean = aircraft_structure.Kff;
    C_mean=aircraft_structure.modeshapes'^-1*C_mod*aircraft_structure.modeshapes^-1;
    %     C_mean=K_tot_mean*0;
%     C_mean(1:end,1:end)=K_tot_mean(1:end,1:end)*0.0035+M_tot_mean(1:end,1:end)*0.0035;

end
    
        gravity=[0 0 0*-9.81]'*gravity_on;
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

%white noise toggle
whiteNoise=0;
whiteNoiseAmplitude=0.01;    %amount of freestream velocity
alphaNoise=0;
betaNoise=0;
velocity=aero_state.V_A;
velocity_t=[];
R_tip=[];
X_tip=[];

for i = 1:length(t_vec)
    
    disp(t_vec(i))
    if whiteNoise && mod(i,5)==0 && i>5
        speedNoise=aero_state.V_A*whiteNoiseAmplitude*[rand(1)-0.5 rand(1)-0.5 rand(1)-0.5];
        alphaNoise=atand(speedNoise(3)/aero_state.V_A);
        betaNoise=atand(speedNoise(2)/aero_state.V_A);
    end
    if i>51 && mean(abs(R_tip(end-50:end)))<flutCrit && mod(i,100)==0
        velocity=velocity+velIncr;
    end
    aircraft_aero=aircraft_aero.solve_time_domain_aerodynamics(aircraft,[0 0 0 0 deg2rad(aero_state.alpha) 0],velocity,aero_state.alpha+alphaNoise,0+betaNoise,...
        [0 0 0] ,rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/constrained_flutter']);
    
   % aircraft_aero.F_body(1,:)=-sind(aero_state.alpha).*aircraft_aero.F_body(3,:);
    
    aircraft=aircraft.compute_beam_forces(aircraft_aero.F_body,aircraft_structure);
    for bb=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(bb),'class_wing')
            aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
        end
    end
   
    %aircraft_structure=aircraft_structure.f_solve();
    aircraft_structure=aircraft_structure.f_set_acceleration([0,0,-9.81,0,0,0],[0 0 0]);
    aircraft_structure=aircraft_structure.f_assemble(1,0);
    if modalAnalysis
        F_tot_mean = modalBasis'*aircraft_structure.Ftest;
    else
        F_tot_mean = aircraft_structure.Ftest;
    end
    
%     % ORIGINAL
%     x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
%         +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));    
    x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
            +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));
    if i>1
        if i>0
            xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
        else
            xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
        end
        xdot_nxt=xdot_now+dt*((1-gamma_newmark)*xdotdot_now+gamma_newmark*xdotdot_nxt);
    end
    
    % FROM class_UVLM_solver.time_domain_flutter_analysis
%     x_nxt=linsolve(aircraft_structure.Mff+beta*dt^2*aircraft_structure.Kff,beta*dt^2*aircraft_structure.Ftest+aircraft_structure.Mff*x_now+dt*aircraft_structure.Mff*xdot_now+dt^2*aircraft_structure.Mff*(1/2-beta)*xdotdot_now);
%     xdotdot_nxt=1/(beta*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta)*xdotdot_now);         %Match with original
%     xdot_nxt=xdot_now+dt*((1-gamma)*xdotdot_now+gamma*xdotdot_nxt);                          %Match with original
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if modalAnalysis
        aircraft_structure.nodal_deflections=modalBasis*x_nxt;
    else
        aircraft_structure.nodal_deflections=x_nxt;
    end
    aircraft_structure=aircraft_structure.f_postprocess();
    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
    
    %compute difference to previous step
    if modalAnalysis
        R=x_nxt-x_now;
        X_tip=[X_tip; x_nxt(1)];
        R_tip=[R_tip; R(1)];
        velocity_t=[velocity_t; velocity];
        if mod(i,breakevery)==0
            figure
            hold on
            subplot(3,1,1)
            plot(X_tip)
            ylabel('Mode 1')
            subplot(3,1,2)
            plot(R_tip)
            ylabel('delta Mode 1')
            subplot(3,1,3)
            plot(velocity_t)
            ylabel('velocity [m/s]')
        end
    else
        R=x_nxt-x_now;
        X_tip=[X_tip; x_nxt(6*0+3)];
        R_tip=[R_tip; R(6*0+3)];
        velocity_t=[velocity_t; velocity];
        if mod(i,breakevery)==0
            figure
            hold on
            subplot(3,1,1)
            plot(X_tip)
            ylabel('tip deflection')
            subplot(3,1,2)
            hold on
            plot(R_tip)
            plot(repmat(flutCrit,i),'r')
            plot(repmat(-flutCrit,i),'r')
            ylabel('delta tip deflection')
            subplot(3,1,3)
            plot(velocity_t)
            ylabel('velocity [m/s]')
        end
    end
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
end
