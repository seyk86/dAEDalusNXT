%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxisPartitioned=simulate_free_flying_mean_axis_partitioned(aircraft,aircraft_structure,n_modes,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

if ctrl_input.signals.dimensions~=length(aircraft.control_surfaces)
    disp('Control Input Signal Incorrect')
end

%newmark beta parameters
beta_newmark=1/4;
gamma_newmark=1/2;

% initialize flight state variables
MeanAxisPartitioned.Xe=zeros(length(t_vec),3);
MeanAxisPartitioned.Euler=zeros(length(t_vec),3);
MeanAxisPartitioned.Vb=zeros(length(t_vec),3);
MeanAxisPartitioned.Ve=zeros(length(t_vec),3);
MeanAxisPartitioned.pqr=zeros(length(t_vec),3);
MeanAxisPartitioned.pqr_dot=zeros(length(t_vec),3);
MeanAxisPartitioned.Alpha=zeros(length(t_vec),1);
MeanAxisPartitioned.Beta=zeros(length(t_vec),1);
MeanAxisPartitioned.V=zeros(length(t_vec),1);
MeanAxisPartitioned.qinf=zeros(length(t_vec),1);
MeanAxisPartitioned.Ma=zeros(length(t_vec),1);

structure_solver_settings=class_wingstructure_solver_settings;
structure_solver_settings.gravity=0;

% for testing
gravity_on=1;

% initialize state vectors to zero
x_now=zeros(length(aircraft_structure.Kff)-6,1);
x_nxt=x_now;
xdot_now=x_now;
xdotdot_now=x_now;
xdotdot_nxt=x_now; 

aircraft_structure.nodal_deflections=zeros(length(aircraft_structure.Kff),1);
aircraft_structure=aircraft_structure.f_postprocess();

Mdiag=aircraft_structure.modeshapes'*aircraft_structure.Mff_lumped*aircraft_structure.modeshapes;
Kdiag=aircraft_structure.modeshapes'*aircraft_structure.Kff*aircraft_structure.modeshapes;
Qdiag=aircraft_structure.modeshapes'*aircraft_structure.Ftest;

structural_acc=[0 0 0 0 0 0];
omega_K=[0 0 0];
Euler=[0 0 0];
X=[0 0 0];


mean_axis_origin=-[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/aircraft.weights.W];
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(Euler',omega_K,aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,[0 0 0]');

Euler=[0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
V_K=-trimmed_state.aerodynamic_state.V_inf;
V_e=[-trimmed_state.aerodynamic_state.V_A 0 0];
MeanAxisPartitioned.Vb(1,:)=[-1 0 0;0 1 0; 0 0 -1]*V_K';
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();

UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=movie_on;
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();
aircraft.reference.p_ref=-(A_Bar(1:3,1:3)'*mean_axis_origin)';%plus or minus?
XHALE_aero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
XHALE_aero=XHALE_aero.initialize_time_domain_solution(dt);

Mass=I_Hat'*A_Bar*aircraft_structure.Mff_lumped*A_Bar'*I_Hat;
Inertia=b_Hat_Skew'*A_Bar*aircraft_structure.Mff_lumped*A_Bar'*b_Hat_Skew;

shape=zeros(length(aircraft_structure.Kff),1);
aircraft=aircraft.compute_beam_forces(XHALE_aero.F_body,aircraft_structure);

for bb=1:length(aircraft_structure.beam)
    if  isa(aircraft_structure.beam(bb),'class_wing')
        aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
    end
end

gravity=[0 0 -9.81]'*gravity_on;

structural_acc(1:3)=structural_acc(1:3)+gravity';
aircraft_structure=aircraft_structure.f_set_acceleration(structural_acc,XHALE_aero.reference.p_ref);
aircraft_structure=aircraft_structure.f_assemble_free(1,0);

Qdiag=aircraft_structure.modeshapes'*aircraft_structure.Ftest;
x_now=Kdiag(7:end,7:end)^-1*Qdiag(7:end);
for ss=1:length(aircraft_structure.modefrequencies)-6
    shape=shape+aircraft_structure.modeshapes(:,6+ss)*x_now(ss);
end
aircraft_structure.nodal_deflections=shape;

aircraft_structure=aircraft_structure.f_postprocess();
alpha_k=trimmed_state.aerodynamic_state.alpha*pi/180;
beta_k=0;

MeanAxisPartitioned.Time=t_vec;

control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;
rho_air=trimmed_state.aerodynamic_state.rho_air;
%load(['results/' aircraft.name '/' maneuver_name '/Simulation_Fixed_Axis']);
disp('Simulating Maneouver: ')
for i=1:length(t_vec)
    fprintf('\b\b\b\b %03d',round(t_vec(i)/t_vec(end)*100));   
    MeanAxisPartitioned.Xe(i,:)=[-1 0 0;0 1 0; 0 0 -1]*X';
    MeanAxisPartitioned.Euler(i,:)=Euler;
    MeanAxisPartitioned.Vb(i,:)=[-1 0 0;0 1 0; 0 0 -1]*V_K';
    MeanAxisPartitioned.Ve(i,:)=[-1 0 0;0 1 0; 0 0 -1]*V_e';
    MeanAxisPartitioned.pqr(i,:)=[-1 0 0;0 1 0; 0 0 -1]*omega_K'; 
    MeanAxisPartitioned.Alpha(i)=alpha_k;
    MeanAxisPartitioned.Beta(i)=beta_k;
    MeanAxisPartitioned.V(i)=norm(V_e);
    MeanAxisPartitioned.qinf(i)=XHALE_aero.qinf;
    
    M_BI=[  cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
        cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
        cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];

    gravity=M_BI*[0 0 -9.81]'*gravity_on;
    
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
    
    %XHALE_aero=XHALE_aero.update_p_ref(-(M_BI*mean_axis_origin)');
    XHALE_aero=XHALE_aero.solve_time_domain_aerodynamics(aircraft,[([-1 0 0;0 1 0; 0 0 -1]*MeanAxisPartitioned.Xe(i,:)')' MeanAxisPartitioned.Euler(i,:)],MeanAxisPartitioned.V(i),MeanAxisPartitioned.Alpha(i),MeanAxisPartitioned.Beta(i),...
        [-1 0 0;0 1 0; 0 0 -1]*MeanAxisPartitioned.pqr(i,:)',rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/mean_axis_partitioned']);
    aircraft=aircraft.compute_beam_forces(XHALE_aero.F_body,aircraft_structure);  

    for bb=1:length(aircraft_structure.beam)
        if  isa(aircraft_structure.beam(bb),'class_wing')
            aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
        end
    end
    
    structural_acc(1:3)=gravity';
    structural_acc(4:6)=0;
    aircraft_structure=aircraft_structure.f_set_acceleration(structural_acc,XHALE_aero.reference.p_ref);
    aircraft_structure=aircraft_structure.f_assemble_free(1,0);     
    Qdiag=aircraft_structure.modeshapes'*aircraft_structure.Ftest;
    
    K_tot_mean=Kdiag(7:end,7:end);
    M_tot_mean=Mdiag(7:end,7:end);
    F_tot_mean=Qdiag(7:end);
    C_mean=K_tot_mean*0.00+M_tot_mean*0.00;
    x_nxt=linsolve(M_tot_mean+beta_newmark*dt^2*K_tot_mean+gamma_newmark*dt*C_mean,beta_newmark*dt^2*F_tot_mean+M_tot_mean*x_now+dt*M_tot_mean*xdot_now+dt^2*M_tot_mean*(1/2-beta_newmark)*xdotdot_now...
        +C_mean*beta_newmark*dt^2*(gamma_newmark/(beta_newmark*dt)*x_now+(gamma_newmark/beta_newmark-1)*xdot_now+1/2*dt*(gamma_newmark/beta_newmark-2)*xdotdot_now));
    xdotdot_nxt=1/(beta_newmark*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta_newmark)*xdotdot_now);
    xdot_nxt=xdot_now+dt*((1-gamma_newmark)*xdotdot_now+gamma_newmark*xdotdot_nxt);   
    
    shape=shape*0;
    for ss=1:length(aircraft_structure.modefrequencies)-6
        shape=shape+aircraft_structure.modeshapes(:,6+ss)*x_nxt(ss);
    end
    aircraft_structure.nodal_deflections=shape;
    aircraft_structure=aircraft_structure.f_postprocess();     
    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);     
  
    MeanAxisPartitioned.cF(i,:)=[XHALE_aero.qinf*XHALE_aero.CX(i)*XHALE_aero.reference.S_ref+gravity(1)*Mass(1,1);
       XHALE_aero.qinf*XHALE_aero.CY(i)*XHALE_aero.reference.S_ref+gravity(2)*Mass(1,1);
       XHALE_aero.qinf*XHALE_aero.CZ(i)*XHALE_aero.reference.S_ref+gravity(3)*Mass(1,1)];   
    MeanAxisPartitioned.cM(i,:)=[XHALE_aero.qinf*XHALE_aero.CL(i)*XHALE_aero.reference.S_ref*XHALE_aero.reference.b_ref;
       XHALE_aero.qinf*XHALE_aero.CM(i)*XHALE_aero.reference.S_ref*XHALE_aero.reference.c_ref;
       XHALE_aero.qinf*XHALE_aero.CN(i)*XHALE_aero.reference.S_ref*XHALE_aero.reference.b_ref];

    F=[sum(aircraft_structure.Ftest(1:6:end));
       sum(aircraft_structure.Ftest(2:6:end));
       sum(aircraft_structure.Ftest(3:6:end))];
   
   MeanAxisPartitioned.fF(i,:)=F;
  
    M= [ XHALE_aero.qinf*XHALE_aero.CL(i)*XHALE_aero.reference.S_ref*XHALE_aero.reference.b_ref;
       -sum((aircraft_structure.node_coords(:,1)-XHALE_aero.state.p_ref(1)).*aircraft_structure.Ftest(3:6:end))+sum((aircraft_structure.node_coords(:,3)-XHALE_aero.state.p_ref(3)).*aircraft_structure.Ftest(1:6:end))+sum(aircraft_structure.Ftest(5:6:end));
       XHALE_aero.qinf*XHALE_aero.CN(i)*XHALE_aero.reference.S_ref*XHALE_aero.reference.b_ref];
   
   
    MeanAxisPartitioned.fM(i,:)=M;
    [X,Euler,V_K,omega_K,V_kin,V_e,alpha_k,beta_k,omega_K_dot,structural_acc]=sixDoF_rigid_body(Mass(1,1),Inertia,Inertia*0,X,Euler,V_K,omega_K,F,M,dt);
    if i>100
        structural_acc=[structural_acc' -omega_K];
    else
        structural_acc=[structural_acc' -omega_K];
    end
    %sav_acc(:,i)=structural_acc;
    x_now=x_nxt;
    xdot_now=xdot_nxt;
    xdotdot_now=xdotdot_nxt;
end

MeanAxisPartitioned.CM=XHALE_aero.CM(2:end);
MeanAxisPartitioned.CZ=XHALE_aero.CZ(2:end);
MeanAxisPartitioned.CX=XHALE_aero.CX(2:end);

MeanAxisPartitioned.Time=t_vec;
save(['results/' aircraft.name '/' maneuver_name '/Simulation_Mean_Axis_Partitioned'],'MeanAxisPartitioned','MeanAxisPartitioned');
end
