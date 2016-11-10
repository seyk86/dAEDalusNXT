%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxisHybrid=simulate_free_flying_hybrid_ms(aircraft,aircraft_structure,AeroelasticSSM,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

mkdir(['results/' aircraft.name],maneuver_name);
mkdir(['results/' aircraft.name '/' maneuver_name],'movie');

MeanAxisHybrid.Xe=zeros(length(t_vec),3);
MeanAxisHybrid.Euler=zeros(length(t_vec),3);
MeanAxisHybrid.Vb=zeros(length(t_vec),3);
MeanAxisHybrid.Ve=zeros(length(t_vec),3);
MeanAxisHybrid.pqr=zeros(length(t_vec),3);
MeanAxisHybrid.pqr_dot=zeros(length(t_vec),3);
MeanAxisHybrid.Alpha=zeros(length(t_vec),1);
MeanAxisHybrid.Beta=zeros(length(t_vec),1);
MeanAxisHybrid.V=zeros(length(t_vec),1);
MeanAxisHybrid.qinf=zeros(length(t_vec),1);
MeanAxisHybrid.Ma=zeros(length(t_vec),1);
n_elastic_modes=AeroelasticSSM.nE;
% for testing
gravity_on=1;
% initialize state vectors to zero

structural_acc=[0 0 0 0 0 0];
omega_K=[0 0 0];
Euler=[0 0 0];
X=[0 0 -1000];
n_elastic_modes=AeroelasticSSM.nE;
mean_axis_origin=-[sum(aircraft_structure.Mff_lumped(1:6:end,1:6:end)*(aircraft_structure.node_coords(:,1)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(2:6:end,2:6:end)*(aircraft_structure.node_coords(:,2)))/aircraft.weights.W;
                      sum(aircraft_structure.Mff_lumped(3:6:end,3:6:end)*(aircraft_structure.node_coords(:,3)))/aircraft.weights.W];
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(Euler',omega_K,aircraft_structure.node_coords,aircraft_structure.nodal_deflections,mean_axis_origin,[0 0 0]');
[M_tot_mean,K_tot_mean,F_tot_mean]=compute_mean_axis_modal_matrices(A_Bar,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,aircraft_structure.Mff,aircraft_structure.Kff,aircraft_structure.modeshapes,AeroelasticSSM.nE+6,zeros(AeroelasticSSM.nE+6,1));
 

Euler=[0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
V_K=trimmed_state.aerodynamic_state.V_inf;
V_e=[trimmed_state.aerodynamic_state.V_A 0 0];
MeanAxisHybrid.Vb(1,:)=[1 0 0;0 1 0; 0 0 1]*V_K';
aircraft.grid_settings.wake=2;

for n_ctrl=1:length(aircraft.control_deflections)
     aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{n_ctrl},trimmed_state.aircraft_state.control_deflections{n_ctrl});  
end
aircraft=aircraft.compute_grid();
aircraft.reference.p_ref=-mean_axis_origin';
aircraft_aero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.reference);
aircraft_aero=aircraft_aero.f_solve_std();

Mass=I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat;
Inertia=b_Hat_Skew'*A_Bar*aircraft_structure.Mff*A_Bar'*b_Hat_Skew;

shape=zeros(length(aircraft_structure.Kff),1);
aircraft=aircraft.compute_beam_forces(aircraft_aero.F_body,aircraft_structure);

gravity=[0 0 9.81]'*gravity_on;

structural_acc(1:3)=structural_acc(1:3)+gravity';
aircraft_structure=aircraft_structure.f_set_acceleration(structural_acc,aircraft_aero.reference.p_ref);
aircraft_structure=aircraft_structure.f_assemble_free(1,0);

aircraft_structure.nodal_deflections=shape*0;

aircraft_structure=aircraft_structure.f_postprocess();
alpha_k=trimmed_state.aerodynamic_state.alpha*pi/180;
beta_k=0;

aircraft=aircraft.compute_grid();
aircraft_aero=aircraft_aero.set_grid(aircraft.grid,aircraft.panels);
elastic_forces=zeros(length(t_vec),6);
rigid_forces=zeros(length(t_vec),6);
x_AE=zeros(length(t_vec),length(AeroelasticSSM.msASee));
x_AE_dot=zeros(length(AeroelasticSSM.msASee),1);
control_surface_states_0=trimmed_state.aircraft_state.control_deflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);
control_surface_states_prv=control_surface_states;

ctrl_state=zeros(1,length(control_surface_states)*3);
ctrl_state_0=ctrl_state;
for n_ctrl=1:length(control_surface_states)
   ctrl_state_0(n_ctrl)=trimmed_state.aircraft_state.control_deflections{n_ctrl};
end
T=[1 0 0;
    0 0 0;
    0 0 0];
Tm=[1 0 0;
    0 1 0;
    0 0 1];

MeanAxisHybrid.Alpha(1)=alpha_k;


x_AE(1,:)=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 MeanAxisHybrid.Alpha(1)*0  0   ([135 0 0]+(Tm*MeanAxisHybrid.Vb(1,:)')')*0   MeanAxisHybrid.pqr(1,:)*0 zeros(1,6)]'-AeroelasticSSM.msASce*(ctrl_state_0')*pi/180*0);
x_trim=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 MeanAxisHybrid.Alpha(1)  0   zeros(1,12)]'-AeroelasticSSM.msASce*(ctrl_state_0')*pi/180);


disp('Simulating Maneouver: ')

%x_trim=[-14 -0.6 -0.02 0.6];
% trim_shape=zeros(length(aircraft_structure.modefrequencies),1);
% for i=1:AeroelasticSSM.nE
%    trim_shape(:)=trim_shape(:)+real(aircraft_structure.modeshapes(1:size(aircraft_structure.modeshapes,1),i+6)).*x_AE(1,i);
% end
% 
% aircraft_structure.nodal_deflections=trim_shape(:);
% aircraft_structure=aircraft_structure.f_postprocess();
% aircraft=aircraft.compute_grid;
% aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
% aircraft.write_grid_deflected([aircraft.name '/Trim_Shape_MSA'],[0 0 0]',[0 alpha_k  0]'*pi/180,[0 0 0]');
% trim_shape=zeros(length(aircraft_structure.modefrequencies),1);
% aircraft_structure.nodal_deflections=trim_shape(:);
% aircraft_structure=aircraft_structure.f_postprocess();
rho_air=trimmed_state.aerodynamic_state.rho_air;

x_stat_ms=AeroelasticSSM.msASee^-1*1/2*trimmed_state.aerodynamic_state.rho_air*trimmed_state.aerodynamic_state.V_A^2*aircraft.reference.S_ref*[zeros(n_elastic_modes,1);AeroelasticSSM.Q_0(7:6+n_elastic_modes,1);zeros(size(AeroelasticSSM.msASre(:,:),1)-2*n_elastic_modes,1)];

for i=1:length(t_vec)
    t_vec(i)
    MeanAxisHybrid.Xe(i,:)=[1 0 0;0 1 0; 0 0 1]*X';
    MeanAxisHybrid.Euler(i,:)=Euler;
    MeanAxisHybrid.Vb(i,:)=[1 0 0;0 1 0; 0 0 1]*V_K';
    MeanAxisHybrid.Ve(i,:)=[1 0 0;0 1 0; 0 0 1]*V_e';
    MeanAxisHybrid.pqr(i,:)=[1 0 0;0 1 0; 0 0 1]*omega_K'; 
    MeanAxisHybrid.Alpha(i)=alpha_k;
    MeanAxisHybrid.Beta(i)=beta_k;
    MeanAxisHybrid.V(i)=norm(V_e);
    [Tmp, a, P, rho_air] = atmosisa(-MeanAxisHybrid.Xe(i,3));
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
            ctrl_state(length(control_surface_states)+n_ctrl)=0*(control_surface_states(n_ctrl)-control_surface_states_prv(n_ctrl))/dt; 
        end
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    end
    
    aero_state=class_aero_state(MeanAxisHybrid.V(i),MeanAxisHybrid.Alpha(i)*180/pi,-MeanAxisHybrid.Beta(i)*180/pi,trimmed_state.aerodynamic_state.Ma,rho_air);
    aircraft_aero=aircraft_aero.f_set_state(aero_state);
    aircraft_aero=aircraft_aero.set_grid(aircraft.grid,aircraft.panels);
    aircraft_aero=aircraft_aero.f_solve_dynamic(-MeanAxisHybrid.pqr(i,1),MeanAxisHybrid.pqr(i,2),-MeanAxisHybrid.pqr(i,3),aircraft_aero.colloc*0);
    MeanAxisHybrid.CX(i)=-aircraft_aero.qinf*aircraft_aero.Cl*sin(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref;
    MeanAxisHybrid.CZ(i)=aircraft_aero.CZ;
    MeanAxisHybrid.CM(i)=aircraft_aero.CM;
        MeanAxisHybrid.CY(i)=aircraft_aero.CY;
    MeanAxisHybrid.CN(i)=aircraft_aero.CN;
        MeanAxisHybrid.CL(i)=aircraft_aero.CL;

    xAE=x_AE(i,:)';  
    AlphaIn=MeanAxisHybrid.Alpha(i)-MeanAxisHybrid.Alpha(1);
    VIn=-0*((Tm*MeanAxisHybrid.Vb(i,:)')'-(Tm*MeanAxisHybrid.Vb(1,:)')');
    %VIn(3)=VIn(3)+(MeanAxisHybrid.Alpha(i)-sin(MeanAxisHybrid.Alpha(1)))*MeanAxisHybrid.Vb(i,1);
    for its=1:30 % fixed point iteration
        x_AE_dot=AeroelasticSSM.msASee*xAE+AeroelasticSSM.msASce*((ctrl_state+ctrl_state_0*0)')*(pi/180) ...
            +AeroelasticSSM.msASre*[0 0 0 MeanAxisHybrid.Euler(1)*0 AlphaIn  MeanAxisHybrid.Beta(i)*0  VIn  (T*0*MeanAxisHybrid.pqr(i,:)')' 0 0 0 0 0 0  ]';
        %    +AeroelasticSSM.ASre*[0 0 0 0 MeanAxisHybrid.Alpha(i) MeanAxisHybrid.Beta(i)  (Tm*MeanAxisHybrid.Vb(1,:)')'-(Tm*MeanAxisHybrid.Vb(i,:)')'  omega_K 0 0 0 0 0 0  ]';
        xAE=x_AE(i,:)'+x_AE_dot*dt;
    end
    x_AE(i+1,:)=x_AE(i,:)'+x_AE_dot*dt;
    
    F_el=-AeroelasticSSM.msCer*((x_AE(i,1:length(AeroelasticSSM.msCer))')+x_stat_ms(1:length(AeroelasticSSM.msCer))+x_trim(1:length(AeroelasticSSM.msCer)));
    elastic_forces(i,:)=F_el;

    F=[aircraft_aero.qinf*aircraft_aero.Cl*sin(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(1)*Mass(1,1);
       aircraft_aero.qinf*aircraft_aero.CY*aircraft_aero.reference.S_ref+gravity(2)*Mass(1,1);
       -aircraft_aero.qinf*aircraft_aero.CZ*aircraft_aero.reference.S_ref+gravity(3)*Mass(1,1)];
   
    M=[-aircraft_aero.qinf*aircraft_aero.CL*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref;
       aircraft_aero.qinf*aircraft_aero.CM*aircraft_aero.reference.S_ref*aircraft_aero.reference.c_ref;
       -aircraft_aero.qinf*aircraft_aero.CN*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref];
   
   rigid_forces(i,1:3)=F;
   rigid_forces(i,1)=rigid_forces(i,1)-gravity(1)*Mass(1,1);
   rigid_forces(i,2)=rigid_forces(i,2)-gravity(2)*Mass(1,1);
   rigid_forces(i,3)=rigid_forces(i,3)-gravity(3)*Mass(1,1);
   rigid_forces(i,4:6)=M;

   [X,Euler,V_K,omega_K,V_kin,V_e,alpha_k,beta_k,omega_K_dot,structural_acc]=sixDoF_rigid_body(Mass(1,1),Inertia,Inertia*0,X,Euler,V_K,omega_K,F-F_el(1:3),M+F_el(4:6),dt);
end

MeanAxisHybrid.rigid_forces=rigid_forces;
MeanAxisHybrid.elastic_forces=elastic_forces;
MeanAxisHybrid.Time=t_vec;
MeanAxisHybrid.x_AE=x_AE;

save(['results/' aircraft.name '/' maneuver_name '/MeanAxisHybrid'],'MeanAxisHybrid','MeanAxisHybrid');
end

