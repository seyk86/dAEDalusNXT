%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxisModal=simulate_free_flying_mean_axis_modal_idcpl(aircraft,aircraft_structure,n_modes,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

if ctrl_input.signals.dimensions~=length(aircraft.control_surfaces)
    disp('Control Input Signal Incorrect')
end

%newmark beta parameters
beta_newmark=1/4;
gamma_newmark=1/2;

% initialize flight state variables
MeanAxisModal.Xe=zeros(length(t_vec),3);
MeanAxisModal.Euler=zeros(length(t_vec),3);
MeanAxisModal.Vb=zeros(length(t_vec),3);
MeanAxisModal.Ve=zeros(length(t_vec),3);
MeanAxisModal.pqr=zeros(length(t_vec),3);
MeanAxisModal.pqr_dot=zeros(length(t_vec),3);
MeanAxisModal.Alpha=zeros(length(t_vec),1);
MeanAxisModal.Beta=zeros(length(t_vec),1);
MeanAxisModal.V=zeros(length(t_vec),1);
MeanAxisModal.qinf=zeros(length(t_vec),1);
MeanAxisModal.Ma=zeros(length(t_vec),1);

% for testing
gravity_on=1;

% initialize state vectors to zero
x_now=zeros(n_modes,1);
x_nxt=x_now;
xdot_now=x_now;
xdotdot_now=x_now;
xdotdot_nxt=x_now;

% aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
% aircraft_structure=aircraft_structure.f_postprocess();
% Euler=[0 0 0]*pi/180;
% 
% aircraft_structure=aircraft_structure.f_solve_free_modes(1,0)
% aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
% [s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections,[0 0 0]',[0 0 0]');
% [M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,aircraft_structure.modeshapes,n_modes,xdot_now(1:end));
%     
% mean_axis_origin=-(I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat)^-1*I_Hat'*(A_Bar*aircraft_structure.Mff*A_Bar'*s_inertial');
% 
% trim & initialize aerodynamics

aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();

aeroelastic_solver_settings=class_aeroelastic_solver_settings;

UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=movie_on;
XHALE_aero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
XHALE_aero=XHALE_aero.initialize_time_domain_solution(dt);

modeshapes_rom=real(aircraft_structure.modeshapes)*[eye(n_modes); zeros(size(aircraft_structure.modeshapes,1)-n_modes,n_modes)];
mdef=aircraft_structure.modeshapes^(-1)*aircraft_structure.nodal_deflections;
x_now(7:end)=mdef(7:n_modes);
x_now(1:6)=[0 0 0 0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
xdot_now(1:3)=-[norm(trimmed_state.aerodynamic_state.V_inf) 0 0];



[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(x_now(4:6),[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections*0,[0 0 0]',[0 0 0]');

MeanAxisModal.Time=t_vec;

x_AE=zeros(length(t_vec),n_modes);

control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;
rho_air=trimmed_state.aerodynamic_state.rho_air;
disp('Simulating Maneouver: ')
for i=1:length(t_vec)
    fprintf('\b\b\b\b %03d',round(t_vec(i)/t_vec(end)*100));
    MeanAxisModal.Xe(i,:)=[-1 0 0;0 1 0; 0 0 -1]*x_now(1:3);
    MeanAxisModal.Euler(i,:)=[-1 0 0;0 1 0; 0 0 -1]*x_now(4:6);
    MeanAxisModal.Vb(i,:)=[-1 0 0;0 1 0; 0 0 -1]*A_Bar(1:3,1:3)'*xdot_now(1:3);
    MeanAxisModal.Ve(i,:)=[-1 0 0;0 1 0; 0 0 -1]*xdot_now(1:3);
    MeanAxisModal.pqr(i,:)=([1      sin(MeanAxisModal.Euler(i,1))*tan(MeanAxisModal.Euler(i,2))         cos(MeanAxisModal.Euler(i,1))*tan(MeanAxisModal.Euler(i,2));
                             0      cos(MeanAxisModal.Euler(i,1))                                       -sin(MeanAxisModal.Euler(i,1));
                             0      sin(MeanAxisModal.Euler(i,1))/cos(MeanAxisModal.Euler(i,2))         cos(MeanAxisModal.Euler(i,1))/cos(MeanAxisModal.Euler(i,2));]^(-1)*[-1 0 0;0 1 0; 0 0 -1]*xdot_now(4:6));
                       
    MeanAxisModal.Alpha(i)=atan(MeanAxisModal.Vb(i,3)/MeanAxisModal.Vb(i,1));
    MeanAxisModal.Beta(i)=atan(MeanAxisModal.Vb(i,2)/(sqrt(MeanAxisModal.Vb(i,1)^2+MeanAxisModal.Vb(i,3)^2)));
    MeanAxisModal.V(i)=norm(MeanAxisModal.Vb(i,:));                  
                       
    mean_axis_origin=-[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/aircraft.weights.W];
    %mean_axis_origin=-(I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat)^-1*I_Hat'*(A_Bar*aircraft_structure.Mff*A_Bar'*b);
    mean_axis_orientation=[0 0 0]'*pi/180;
    [s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(x_now(4:6),xdot_now(4:6),aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,mean_axis_orientation);
    
    [M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices_dcpl(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,aircraft_structure.modeshapes,n_modes,xdot_now);
    
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
    
    XHALE_aero=XHALE_aero.solve_time_domain_aerodynamics(aircraft,[([-1 0 0;0 1 0; 0 0 -1]*MeanAxisModal.Xe(i,:)')' MeanAxisModal.Euler(i,:)],MeanAxisModal.V(i),MeanAxisModal.Alpha(i),MeanAxisModal.Beta(i),...
        [-1 0 0;0 1 0; 0 0 -1]*MeanAxisModal.pqr(i,:)',rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/mean_axis_modal']);
    
    aircraft=aircraft.compute_beam_forces(XHALE_aero.F_body,aircraft_structure);  
    for bb=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(bb),'class_wing')
        aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
        end
    end
    
    gravity=A_Bar(1:3,1:3)'*[0 0 -9.81]'*gravity_on;
    
    for bb=1:length(aircraft_structure.beam)
        for be=1:length(aircraft_structure.beam(bb).beamelement)
            aircraft_structure.beam(bb).beamelement(be).ax=gravity(1);
            aircraft_structure.beam(bb).beamelement(be).ay=gravity(2);
            aircraft_structure.beam(bb).beamelement(be).az=gravity(3);
            aircraft_structure.beam(bb).update_Q=1;
        end
    end  
 
    aircraft_structure=aircraft_structure.f_assemble_free(1,0);
    
    Fm=modeshapes_rom(1:end,7:end)'*aircraft_structure.Ftest;
    F_tot_mean=F_tot_mean+[I_Hat'*A_Bar*aircraft_structure.Ftest;b_Hat_Skew'*A_Bar*aircraft_structure.Ftest;Fm];

    % no structural damping
    C_mean=K_tot_mean*0.001+M_tot_mean*0.001;
    
    x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
        +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));
    xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
    xdot_nxt=xdot_now+dt*((1-gamma_newmark)*xdotdot_now+gamma_newmark*xdotdot_nxt);   

    aircraft_structure.nodal_deflections=modeshapes_rom(1:end,7:end)*x_nxt(7:end);
    
    aircraft_structure=aircraft_structure.f_postprocess();
    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    
   % step forward
    x_now=x_nxt;
    x_AE(i,:)=x_now;
    xdot_now=xdot_nxt;
    xdotdot_now=xdotdot_nxt;
end
MeanAxisModal.CM=XHALE_aero.CM;
MeanAxisModal.CZ=XHALE_aero.CZ;
MeanAxisModal.x_AE=x_AE;
save('MeanAxisModal','MeanAxisModal','MeanAxisModal');
%save('AeroelasticModalStates','x_AE','x_AE')
save(['results/' aircraft.name '/' maneuver_name '/Simulation_Mean_Axis_Modal_idcpl'],'MeanAxisModal','MeanAxisModal');
%save(['results/' aircraft.name '/' maneuver_name '/AeroelasticModalStates'],'x_AE','x_AE');
end

