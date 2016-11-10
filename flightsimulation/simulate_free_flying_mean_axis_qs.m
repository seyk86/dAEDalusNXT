%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxis=simulate_free_flying_mean_axis_qs(aircraft,aircraft_structure,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

if ctrl_input.signals.dimensions~=length(aircraft.control_surfaces)
    disp('Control Input Signal Incorrect')
end

%newmark beta parameters
beta_newmark=1/4;
gamma_newmark=1/2;

% initialize flight state variables
MeanAxis.Xe=zeros(length(t_vec),3);
MeanAxis.Euler=zeros(length(t_vec),3);
MeanAxis.Vb=zeros(length(t_vec),3);
MeanAxis.Ve=zeros(length(t_vec),3);
MeanAxis.pqr=zeros(length(t_vec),3);
MeanAxis.pqr_dot=zeros(length(t_vec),3);
MeanAxis.Alpha=zeros(length(t_vec),1);
MeanAxis.Beta=zeros(length(t_vec),1);
MeanAxis.V=zeros(length(t_vec),1);
MeanAxis.qinf=zeros(length(t_vec),1);
MeanAxis.Ma=zeros(length(t_vec),1);

% for testing
gravity_on=1;

% initialize state vectors to zero
x_now=zeros(length(aircraft_structure.Kff)+6,1);
x_nxt=x_now;
xdot_now=x_now;
xdotdot_now=x_now;
xdotdot_nxt=x_now;

%aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
%aircraft_structure=aircraft_structure.f_postprocess();
%
%
%[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections,[0 0 0]',[0 0 0]');
%[M_tot_mean,K_tot_mean,F_acc_mean]=compute_mean_axis_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,xdot_now);
%mean_axis_origin=-(I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat)^-1*I_Hat'*(A_Bar*aircraft_structure.Mff*A_Bar'*s_inertial');
%[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(x_now(4:6),xdot_now(4:6),aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,[0 0 0]');
%[M_tot_mean,K_tot_mean,F_acc_mean]=compute_mean_axis_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,xdot_now);
 
% set current weight
%aircraft.weights.W=M_tot_mean(1,1);
%trimmed_state.aircraft_state.weight=M_tot_mean(1,1);
%trimmed_state.aircraft_state.CG_ref=-mean_axis_origin';
% trim & initialize aerodynamics

aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();

aeroelastic_solver_settings=class_aeroelastic_solver_settings;

% UVLM_settings=class_UVLM_computation_settings();
% UVLM_settings.debug=0;
% UVLM_settings.movie=movie_on;
% aircraft_aero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
% aircraft_aero=aircraft_aero.initialize_time_domain_solution(dt);
aircraft_aero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.reference);
aircraft_aero=aircraft_aero.f_solve_std();


x_now(7:end)=aircraft_structure.nodal_deflections;
x_now(1:6)=[0 0 0 0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
xdot_now(1:3)=-[norm(trimmed_state.aerodynamic_state.V_inf) 0 0];

[~,~,~,~,~,~,A_Bar]=compute_eqm_matrices_mean(x_now(4:6),[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections,[0 0 0]',[0 0 0]');
MeanAxis.Time=t_vec;

control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;
rho_air=trimmed_state.aerodynamic_state.rho_air;
disp('Simulating Maneouver: ')
for i=1:length(t_vec)
    fprintf('\b\b\b\b %03d',round(t_vec(i)/t_vec(end)*100));
    % compute important states
    MeanAxis.Xe(i,:)=[-1 0 0;0 1 0; 0 0 -1]*x_now(1:3);
    MeanAxis.Euler(i,:)=[-1 0 0;0 1 0; 0 0 -1]*x_now(4:6);
    MeanAxis.Vb(i,:)=[-1 0 0;0 1 0; 0 0 -1]*A_Bar(1:3,1:3)'*xdot_now(1:3);
    MeanAxis.Ve(i,:)=[-1 0 0;0 1 0; 0 0 -1]*xdot_now(1:3);
    MeanAxis.pqr(i,:)=([1   sin(MeanAxis.Euler(i,1))*tan(MeanAxis.Euler(i,2))       cos(MeanAxis.Euler(i,1))*tan(MeanAxis.Euler(i,2));
                        0   cos(MeanAxis.Euler(i,1))                                -sin(MeanAxis.Euler(i,1));
                        0   sin(MeanAxis.Euler(i,1))/cos(MeanAxis.Euler(i,2))       cos(MeanAxis.Euler(i,1))/cos(MeanAxis.Euler(i,2));]^(-1)*[-1 0 0;0 1 0; 0 0 -1]*xdot_now(4:6));
                       
    MeanAxis.Alpha(i)=atan(MeanAxis.Vb(i,3)/MeanAxis.Vb(i,1));
    MeanAxis.Beta(i)=atan(MeanAxis.Vb(i,2)/(sqrt(MeanAxis.Vb(i,1)^2+MeanAxis.Vb(i,3)^2)));
    MeanAxis.V(i)=norm(MeanAxis.Vb(i,:));
    
    mean_axis_origin=[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/aircraft.weights.W];
    mean_axis_orientation=[0 0 0]'*pi/180;
    [s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(x_now(4:6),xdot_now(4:6),aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,mean_axis_orientation);
    
    [M_tot_mean,K_tot_mean,F_acc_mean]=compute_mean_axis_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,xdot_now);
    
    % set control surface deflections
    control_surface_states_prv=control_surface_states;
    for n_ctrl=1:length(control_surface_states)
        control_surface_states(n_ctrl)=interp1(ctrl_input.time,ctrl_input.signals.values(:,n_ctrl),t_vec(i));
    end
    % only update if changed to previous step
    if sum(control_surface_states_prv==control_surface_states)~=length(control_surface_states)
        for n_ctrl=1:length(control_surface_states)
            aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{n_ctrl},control_surface_states(n_ctrl)+control_surface_states_0{n_ctrl});   
        end
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    end
    
    aero_state=class_aero_state(MeanAxis.V(i),MeanAxis.Alpha(i)*180/pi,-MeanAxis.Beta(i)*180/pi,0,rho_air);
    aircraft_aero=aircraft_aero.f_set_state(aero_state);
    aircraft_aero=aircraft_aero.set_grid(aircraft.grid,aircraft.panels);
    aircraft_aero=aircraft_aero.f_solve_dynamic(-MeanAxis.pqr(i,1),MeanAxis.pqr(i,2),-MeanAxis.pqr(i,3),aircraft_aero.colloc*0);
%     aircraft_aero=aircraft_aero.solve_time_domain_aerodynamics(aircraft,[([-1 0 0;0 1 0; 0 0 -1]*MeanAxis.Xe(i,:)')' MeanAxis.Euler(i,:)],MeanAxis.V(i),MeanAxis.Alpha(i),MeanAxis.Beta(i),...
%         [-1 0 0;0 1 0; 0 0 -1]*MeanAxis.pqr(i,:)',rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/mean_axis']);

    aircraft=aircraft.compute_beam_forces(aircraft_aero.F_body,aircraft_structure);  
    
    for bb=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(bb),'class_wing')
            aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
        end
    end
    % graviation
    gravity=A_Bar(1:3,1:3)'*[0 0 -9.81]'*gravity_on;
    for bb=1:length(aircraft_structure.beam)
        for be=1:length(aircraft_structure.beam(bb).beamelement)
            aircraft_structure.beam(bb).beamelement(be).ax=gravity(1);
            aircraft_structure.beam(bb).beamelement(be).ay=gravity(2);
            aircraft_structure.beam(bb).beamelement(be).az=gravity(3);
            aircraft_structure.beam(bb).update_Q=1;
        end
    end  
    % assemble
    aircraft_structure=aircraft_structure.f_assemble_free(1,0);
    % plus inertial forces
    F_tot_mean=F_acc_mean+[I_Hat'*A_Bar*aircraft_structure.Ftest;b_Hat_Skew'*A_Bar*aircraft_structure.Ftest;aircraft_structure.Ftest];
     % no structural damping
    C_mean=K_tot_mean*0;
    C_mean(7:end,7:end)=K_tot_mean(7:end,7:end)*0.001+M_tot_mean(7:end,7:end)*0.001;
    % newmark beta scheme
    x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
        +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));
    xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
    xdot_nxt=xdot_now+dt*((1-gamma_newmark)*xdotdot_now+gamma_newmark*xdotdot_nxt);   
    % set physical deflections in structural mesh
    aircraft_structure.nodal_deflections=x_nxt(7:end);
    aircraft_structure=aircraft_structure.f_postprocess();
    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    % step forward
    x_now=x_nxt;
    xdot_now=xdot_nxt;
    xdotdot_now=xdotdot_nxt;
end

MeanAxis.CM=aircraft_aero.CM;
MeanAxis.CZ=aircraft_aero.CZ;
MeanAxis.CX=aircraft_aero.CX;
save(['results/' aircraft.name '/' maneuver_name '/Simulation_Mean_AxisQS'],'MeanAxis','MeanAxis');
end
