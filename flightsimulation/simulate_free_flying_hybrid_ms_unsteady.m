%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxisHybrid=simulate_free_flying_hybrid_ms_unsteady(aircraft,aircraft_structure,AeroelasticSSM,trimmed_state,ctrl_input,t_vec,dt,maneuver_name,movie_on)

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
rtf=10;
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
[s_inertial,delta_inertial,I_Hat,b,b_Hat_Skew,omega_Hat_Skew,A_Bar]=compute_eqm_matrices_mean(Euler',omega_K,aircraft_structure.node_coords,aircraft_structure.nodal_deflections*0,mean_axis_origin,[0 0 0]');

Euler=[0 trimmed_state.aerodynamic_state.alpha*pi/180 0];
V_K=trimmed_state.aerodynamic_state.V_inf;
V_e=[trimmed_state.aerodynamic_state.V_A 0 0];
MeanAxisHybrid.Vb(1,:)=[1 0 0;0 1 0; 0 0 1]*V_K';
aircraft.grid_settings.wake=2;
for n_ctrl=1:length(aircraft.control_deflections)
     aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{n_ctrl},trimmed_state.aircraft_state.control_deflections{n_ctrl});  
end
aircraft=aircraft.compute_grid();
aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);

aircraft=aircraft.compute_grid();
aircraft.reference.p_ref=-mean_axis_origin';


Mass=I_Hat'*A_Bar*aircraft_structure.Mff*A_Bar'*I_Hat;
Inertia=b_Hat_Skew'*A_Bar*aircraft_structure.Mff*A_Bar'*b_Hat_Skew;

alpha_k=trimmed_state.aerodynamic_state.alpha*pi/180;
beta_k=0;
aircraft=aircraft.compute_grid();
UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.movie=movie_on;
aircraft.reference.p_ref=-mean_axis_origin';
aircraft_aero=class_UVLM_solver(aircraft.name,aircraft.grid,aircraft.is_te,aircraft.panels,trimmed_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
aircraft_aero=aircraft_aero.initialize_time_domain_solution(dt*rtf);

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

MeanAxisHybrid.Alpha(1)=alpha_k;

% aircraft_structure.nodal_deflections=trim_shape(:);
% aircraft_structure=aircraft_structure.f_postprocess();
% aircraft=aircraft.compute_grid;
% aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
% aircraft.write_grid_deflected([aircraft.name '/Trim_Shape_MSA'],[0 0 0]',[0 alpha_k  0]'*pi/180,[0 0 0]');
% trim_shape=zeros(length(aircraft_structure.modefrequencies),1);
% aircraft_structure.nodal_deflections=trim_shape(:);
% aircraft_structure=aircraft_structure.f_postprocess();
rho_air=trimmed_state.aerodynamic_state.rho_air;




x_AE(1,:)=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 MeanAxisHybrid.Alpha(1)*0  0   ([135 0 0]+(Tm*MeanAxisHybrid.Vb(1,:)')')*0   MeanAxisHybrid.pqr(1,:)*0 zeros(1,6)]'-AeroelasticSSM.msASce*(ctrl_state_0')*pi/180*0);
elastic_forces=zeros(length(t_vec),6);
rigid_forces=zeros(length(t_vec),6);
disp('Simulating Maneouver: ')

% %x_trim=[-14 -0.6 -0.02 0.6];
% trim_shape=zeros(length(aircraft_structure.modefrequencies),1);
% for i=1:AeroelasticSSM.nE
%    trim_shape(:)=trim_shape(:)+real(aircraft_structure.modeshapes(1:size(aircraft_structure.modeshapes,1),i+6)).*x_AE(1,i);
% end

x_stat_ms=AeroelasticSSM.msASee^-1*1/2*trimmed_state.aerodynamic_state.rho_air*trimmed_state.aerodynamic_state.V_A^2*aircraft.reference.S_ref*[zeros(n_elastic_modes,1);AeroelasticSSM.Q_0(7:6+n_elastic_modes,1);zeros(size(AeroelasticSSM.msASre(:,:),1)-2*n_elastic_modes,1)];
x_trim_ms=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 alpha_k  0  zeros(1,12)]'-AeroelasticSSM.msASce*(ctrl_state_0')*pi/180);

for i=1:length(t_vec)
    t_vec(i)
    MeanAxisHybrid.Xe(i,:)=[1 0 0;0 1 0; 0 0 1]*X';
    MeanAxisHybrid.Euler(i,:)=Euler;
    MeanAxisHybrid.Vb(i,:)=[1 0 0;0 1 0; 0 0 1]*V_K';
    MeanAxisHybrid.Ve(i,:)=[1 0 0;0 1 0; 0 0 1]*V_e';
    MeanAxisHybrid.pqr(i,:)=[1 0 0;0 1 0; 0 0 1]*omega_K'; 
    MeanAxisHybrid.Alpha(i)=alpha_k;
    MeanAxisHybrid.Beta(i)=beta_k;
    MeanAxisHybrid.V(i)=norm(V_K);
    MeanAxisHybrid.qinf(i)=aircraft_aero.qinf;
    [T, a, P, rho_air] = atmosisa(-MeanAxisHybrid.Xe(i,3));
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
            ctrl_state(length(control_surface_states)+n_ctrl)=(control_surface_states(n_ctrl)-control_surface_states_prv(n_ctrl))*0/dt; 
        end
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
    end
    
    if mod(i,rtf)==0
        aircraft_aero=aircraft_aero.solve_time_domain_aerodynamics(aircraft,[([-1 0 0;0 1 0; 0 0 -1]*MeanAxisHybrid.Xe(i,:)')' MeanAxisHybrid.Euler(i,:)],MeanAxisHybrid.V(i),MeanAxisHybrid.Alpha(i),-MeanAxisHybrid.Beta(i),...
            [-1 0 0;0 1 0; 0 0 -1]*MeanAxisHybrid.pqr(i,:)',rho_air,i,['results/' aircraft.name '/' maneuver_name '/movie/mean_axis_rigid_unsteady']);
    end
    MeanAxisHybrid.qinf(i)=1/2*rho_air*MeanAxisHybrid.V(i)^2;
    
    xAE=x_AE(i,:)';  
    AlphaIn=MeanAxisHybrid.Alpha(i)-MeanAxisHybrid.Alpha(1);
    VIn=-0*((Tm*MeanAxisHybrid.Vb(i,:)')'-(Tm*MeanAxisHybrid.Vb(1,:)')');
    %VIn(3)=VIn(3)+(MeanAxisHybrid.Alpha(i)-sin(MeanAxisHybrid.Alpha(1)))*MeanAxisHybrid.Vb(i,1);
    for its=1:30 % fixed point iteration
        x_AE_dot=AeroelasticSSM.msASee*xAE+AeroelasticSSM.msASce*((ctrl_state+ctrl_state_0*0)')*(pi/180) ...
            +AeroelasticSSM.msASre*[0 0 0 0 AlphaIn  MeanAxisHybrid.Beta(i)*0  VIn  (T*0.0*MeanAxisHybrid.pqr(i,:)')' 0 0 0 0 0 0  ]';
        %    +AeroelasticSSM.ASre*[0 0 0 0 MeanAxisHybrid.Alpha(i) MeanAxisHybrid.Beta(i)  (Tm*MeanAxisHybrid.Vb(1,:)')'-(Tm*MeanAxisHybrid.Vb(i,:)')'  omega_K 0 0 0 0 0 0  ]';
        xAE=x_AE(i,:)'+x_AE_dot*dt;
    end
    x_AE(i+1,:)=x_AE(i,:)'+x_AE_dot*dt;
    
    F_el=-AeroelasticSSM.msCer*((x_AE(i,1:length(AeroelasticSSM.msCer))')+x_trim_ms(1:length(AeroelasticSSM.msCer))+x_stat_ms(1:length(AeroelasticSSM.msCer)));
    elastic_forces(i,:)=F_el;
    if mod(i,rtf)==0
        F=[aircraft_aero.qinf*aircraft_aero.CZ(i/rtf)*sin(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(1)*Mass(1,1);
            aircraft_aero.qinf*aircraft_aero.CY(i/rtf)*aircraft_aero.reference.S_ref+gravity(2)*Mass(1,1);
            -aircraft_aero.qinf*aircraft_aero.CZ(i/rtf)*cos(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(3)*Mass(1,1)];
        
        M=[-aircraft_aero.qinf*aircraft_aero.CL(i/rtf)*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref;
            aircraft_aero.qinf*aircraft_aero.CM(i/rtf)*aircraft_aero.reference.S_ref*aircraft_aero.reference.c_ref;
            -aircraft_aero.qinf*aircraft_aero.CN(i/rtf)*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref];
    end

    if i==1
        F=[aircraft_aero.qinf*aircraft_aero.CZ(1)*sin(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(1)*Mass(1,1);
            aircraft_aero.qinf*aircraft_aero.CY(1)*aircraft_aero.reference.S_ref+gravity(2)*Mass(1,1);
            -aircraft_aero.qinf*aircraft_aero.CZ(1)*cos(MeanAxisHybrid.Alpha(i))*aircraft_aero.reference.S_ref+gravity(3)*Mass(1,1)];
        
        M=[-aircraft_aero.qinf*aircraft_aero.CL(1)*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref;
            aircraft_aero.qinf*aircraft_aero.CM(1)*aircraft_aero.reference.S_ref*aircraft_aero.reference.c_ref;
            -aircraft_aero.qinf*aircraft_aero.CN(1)*aircraft_aero.reference.S_ref*aircraft_aero.reference.b_ref];  
    end
    
    if isnan(Euler(2))
        break;
    end

   [X,Euler,V_K,omega_K,V_kin,V_e,alpha_k,beta_k,omega_K_dot,structural_acc]=sixDoF_rigid_body(Mass(1,1),Inertia,Inertia*0,X,Euler,V_K,omega_K,F-F_el(1:3),M+F_el(4:6),dt);
end

MeanAxisHybrid.rigid_forces=rigid_forces;
MeanAxisHybrid.elastic_forces=elastic_forces;
MeanAxisHybrid.Time=t_vec;
MeanAxisHybrid.x_AE=x_AE;
    MeanAxisHybrid.CX=-aircraft_aero.CZ*sin(MeanAxisHybrid.Alpha(i));%XHALE_aero.CX;
    MeanAxisHybrid.CZ=aircraft_aero.CZ;
    MeanAxisHybrid.CM=aircraft_aero.CM;
    MeanAxisHybrid.CY=aircraft_aero.CY;
    MeanAxisHybrid.CL=aircraft_aero.CL;
    MeanAxisHybrid.CN=aircraft_aero.CN;
save(['results/' aircraft.name '/' maneuver_name '/MeanAxisHybrid'],'MeanAxisHybrid','MeanAxisHybrid');
end

