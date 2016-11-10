%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function Rigid=simulate_free_flying_rigid(aircraft,aircraft_structure,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

% initialize flight state variables
Rigid.Xe=zeros(length(t_vec),3);
Rigid.Euler=zeros(length(t_vec),3);
Rigid.Vb=zeros(length(t_vec),3);
Rigid.Ve=zeros(length(t_vec),3);
Rigid.pqr=zeros(length(t_vec),3);
Rigid.pqr_dot=zeros(length(t_vec),3);
Rigid.Alpha=zeros(length(t_vec),1);
Rigid.Beta=zeros(length(t_vec),1);
Rigid.V=zeros(length(t_vec),1);
Rigid.qinf=zeros(length(t_vec),1);
Rigid.Ma=zeros(length(t_vec),1);

% for testing
gravity_on=1;
% initialize state vectors to zero

[~,~,~,rho_itr]=atmosisa(0);
height=0;
while (trimmed_state.aerodynamic_state.rho_air-rho_itr)<0.0001;
    [~,~,~,rho_itr]=atmosisa(height);
    height=height+0.5;
end


structural_acc=[0 0 0 0 0 0];
omega_K=[0 0 0];
Euler=[0 trimmed_state.aerodynamic_state.alpha 0]*pi/180;
X=[0 0 -height];

mean_axis_origin=[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/aircraft.weights.W;
    sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/aircraft.weights.W*0;
    sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/aircraft.weights.W];
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(Euler',omega_K,aircraft_structure.node_coords,aircraft_structure.nodal_deflections*0,mean_axis_origin,[0 0 0]');

%Euler=[0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
V_K=trimmed_state.aerodynamic_state.V_inf;
V_e=[trimmed_state.aerodynamic_state.V_A 0 0];
Rigid.Vb(1,:)=[1 0 0;0 1 0; 0 0 1]*V_K';

aircraft.grid_settings.wake=2;
for n_ctrl=1:length(aircraft.control_deflections)
    aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{n_ctrl},trimmed_state.aircraft_state.control_deflections{n_ctrl});
end

aircraft=aircraft.compute_grid();
aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
aircraft=aircraft.compute_grid();
aircraft.reference.p_ref=mean_axis_origin';
aircraft_aero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.reference);
aircraft_aero=aircraft_aero.f_solve_std();


Mass=I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat;
Inertia=b_Hat_Skew'*A_Bar*aircraft_structure.Mff*A_Bar'*b_Hat_Skew;

shape=zeros(length(aircraft_structure.Kff),1);
aircraft=aircraft.compute_beam_forces(aircraft_aero.F_body,aircraft_structure);

gravity=[0 0 9.81]'*gravity_on;

structural_acc(1:3)=structural_acc(1:3)+gravity';
% aircraft_structure=aircraft_structure.f_set_acceleration(structural_acc,XHALE_aero.reference.p_ref);
% aircraft_structure=aircraft_structure.f_assemble_free(1,0);
%
% aircraft_structure.nodal_deflections=shape*0;
%
% aircraft_structure=aircraft_structure.f_postprocess();
alpha_k=trimmed_state.aerodynamic_state.alpha*pi/180;
beta_k=0;

aircraft=aircraft.compute_grid();
aircraft_aero=aircraft_aero.set_grid(aircraft.grid,aircraft.panels);
elastic_forces=zeros(length(t_vec),6);
rigid_forces=zeros(length(t_vec),6);

control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;

ctrl_state=zeros(1,length(control_surface_states)*3);
ctrl_state_0=ctrl_state;
for n_ctrl=1:length(control_surface_states)
    ctrl_state_0(n_ctrl)=trimmed_state.aircraft_state.control_deflections{n_ctrl};
end
T=[-1 0 0;
    0 1 0;
    0 0 -1];
Tm=[1 0 0;
    0 1 0;
    0 0 1];
%load(['results/' aircraft.name '/' maneuver_name '/Simulation_Fixed_AxisQS'])
rho_air=trimmed_state.aerodynamic_state.rho_air;
for i=1:length(t_vec)
    t_vec(i)
    Rigid.Xe(i,:)=[1 0 0;0 1 0; 0 0 1]*X';
    Rigid.Euler(i,:)=Euler;
    Rigid.Vb(i,:)=[1 0 0;0 1 0; 0 0 1]*V_K';
    Rigid.Ve(i,:)=[1 0 0;0 1 0; 0 0 1]*V_e';
    Rigid.pqr(i,:)=[1 0 0;0 1 0; 0 0 1]*omega_K';
    Rigid.Alpha(i)=alpha_k;
    Rigid.Beta(i)=beta_k;
    Rigid.V(i)=norm(V_K);
    Rigid.qinf(i)=aircraft_aero.qinf;
    [T, a, P, rho_air] = atmosisa(-Rigid.Xe(i,3));
    M_BI=[  cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
        cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
        cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
    
    gravity=M_BI*[0 0 9.81]'*gravity_on;
    
    % set control surface deflections
    control_surface_states_prv=control_surface_states;
    for n_ctrl=1:length(control_surface_states)
        control_surface_states(n_ctrl)=interp1(ctrl_input.time,ctrl_input.signals.values(:,n_ctrl),t_vec(i));
    end
    % only update if changed to previous step
    if sum(control_surface_states_prv==control_surface_states)~=length(control_surface_states)
        for n_ctrl=1:length(control_surface_states)
            aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{n_ctrl},control_surface_states(n_ctrl)+control_surface_states_0{n_ctrl});
            ctrl_state(n_ctrl)=control_surface_states(n_ctrl);
            ctrl_state(length(control_surface_states)+n_ctrl)=(control_surface_states(n_ctrl)-control_surface_states_prv(n_ctrl))/dt;
        end
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    end
    
    aero_state=class_aero_state(Rigid.V(i),Rigid.Alpha(i)*180/pi,-Rigid.Beta(i)*180/pi,trimmed_state.aerodynamic_state.Ma,rho_air);
    aero_state.p_ref=mean_axis_origin';
    aircraft_aero=aircraft_aero.f_set_state(aero_state);
    aircraft_aero=aircraft_aero.set_grid(aircraft.grid,aircraft.panels);
    aircraft_aero=aircraft_aero.f_solve_dynamic(-Rigid.pqr(i,1),Rigid.pqr(i,2),-Rigid.pqr(i,3),aircraft_aero.colloc*0);
    Rigid.CX(i)=-aircraft_aero.Cl*sin(Rigid.Alpha(i));%XHALE_aero.CX;
    Rigid.CZ(i)=aircraft_aero.CZ;
    Rigid.CM(i)=aircraft_aero.CM;
    Rigid.CY(i)=aircraft_aero.CY;
    Rigid.CL(i)=aircraft_aero.CL;
    Rigid.CN(i)=aircraft_aero.CN;
    Rigid.qinf(i)=1/2*rho_air*Rigid.V(i)^2;
    
    F=[aircraft_aero.qinf*aircraft_aero.Cl*sin(Rigid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(1)*Mass(1,1);
        aircraft_aero.qinf*aircraft_aero.CY*aircraft_aero.reference.S_ref+gravity(2)*Mass(1,1);
        -aircraft_aero.qinf*aircraft_aero.Cl*cos(Rigid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(3)*Mass(1,1)-trimmed_state.aircraft_state.engine_thrust*length(aircraft.engines)*sin(atan(aircraft_aero.Uinf(3)/aircraft_aero.Uinf(1)))];
    
    M=[-aircraft_aero.qinf*aircraft_aero.CL*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref;
        aircraft_aero.qinf*aircraft_aero.CM*aircraft_aero.reference.S_ref*aircraft_aero.reference.c_ref-(aircraft.engines(1).cg_pos(3)-trimmed_state.aerodynamic_state.p_ref(3))*trimmed_state.aircraft_state.engine_thrust*length(aircraft.engines);
        -aircraft_aero.qinf*aircraft_aero.CN*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref];

    
    Rigid.cF(i,:)=F;
    Rigid.cM(i,:)=M;
    
    rigid_forces(i,1:3)=F;
    rigid_forces(i,1)=rigid_forces(i,1)-gravity(1)*Mass(1,1);
    rigid_forces(i,2)=rigid_forces(i,2)-gravity(2)*Mass(1,1);
    rigid_forces(i,3)=rigid_forces(i,3)-gravity(3)*Mass(1,1);
    rigid_forces(i,4:6)=M;
    
    [X,Euler,V_K,omega_K,V_kin,V_e,alpha_k,beta_k,omega_K_dot,structural_acc]=sixDoF_rigid_body(Mass(1,1),Inertia,Inertia*0,X,Euler,V_K,omega_K,F,M,dt);
end

Rigid.rigid_forces=rigid_forces;
Rigid.elastic_forces=elastic_forces;
Rigid.Time=t_vec;

Rigid=Rigid;

save(['results/' aircraft.name '/' maneuver_name '/Rigid'],'Rigid','Rigid');
end

