%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,flight_state,wingaero,x_AE ] = trim_aircraft_AESSM(aircraft,flight_state,AeroelasticSSM,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
trim_name3=[];
trim_name4=[];
consider_thrust=1;
for i=1:length(flight_state.aircraft_state.control_deflections)
    %if (flight_state.aircraft_state.control_deflections{i}>=1E-4)
        aircraft=aircraft.f_set_control_surface(flight_state.aircraft_state.control_surfaces{i},flight_state.aircraft_state.control_deflections{i});
   % end
end

def=0;
elev=0;
d_elev=2;
no_redeform=0;
trim_name=[];
trim_name2=[];
trim_name3=[];
trim_name4=[];
%aircraft.reference.p_ref=flight_state.aircraft_state.CG_ref;
trim_ref=aircraft.reference;
trim_ref.p_ref=flight_state.aircraft_state.CG_ref;
if nargin==5
    wingstructure=varargin{1};
    def=wingstructure.f_get_deflections_c4;
end
if nargin==7
    wingstructure=varargin{1};
    elev=0;
    def=varargin{4};
end
trim_idx=1;

if isempty(aircraft.trim_surfaces)
    trim_name='elevator';
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
else
    trim_name=aircraft.trim_surfaces{1};
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
end

elev=aircraft.control_deflections{trim_idx};
if length(aircraft.trim_surfaces)==2
    trim_name2=aircraft.trim_surfaces{2};
end

if length(aircraft.trim_surfaces)==4
    trim_name2=aircraft.trim_surfaces{2};
    trim_name3=aircraft.trim_surfaces{3};
    trim_name4=aircraft.trim_surfaces{4};
end

flight_state.aerodynamic_state.p_ref=flight_state.aircraft_state.CG_ref;
aircraft=aircraft.f_set_control_surface(trim_name,elev);
if ~isempty(trim_name2)
   aircraft=aircraft.f_set_control_surface(trim_name2,-elev); 
end
if ~isempty(trim_name3)
   aircraft=aircraft.f_set_control_surface(trim_name3,elev); 
end

if ~isempty(trim_name4)
   aircraft=aircraft.f_set_control_surface(trim_name4,elev); 
end

aircraft=aircraft.compute_grid();

if nargin>3
    aircraft=aircraft.compute_deflected_grid(def);
    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);   
else
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);
end
disp(['     	trimming aircraft for Cl: ' num2str(flight_state.get_Cl(aircraft.reference.S_ref))]);
wingaero=wingaero.f_solve_for_Cl(flight_state.get_Cl(aircraft.reference.S_ref));



aircraft=aircraft.compute_CD_f(flight_state.aerodynamic_state,aircraft.reference.S_ref);
if consider_thrust==1
    wingaero=wingaero.f_solve_full();
    T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)*cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
    delta_CM=0;
    delta_Cl=0;
    for i=1:length(aircraft.engines)
        aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end


x_AE=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 wingaero.state.alpha*pi/180  0   zeros(1,6)  zeros(1,6)]'-AeroelasticSSM.msASce*[cell2mat(aircraft.control_deflections)';cell2mat(aircraft.control_deflections)'*0;cell2mat(aircraft.control_deflections)'*0]*pi/180);
F_el=AeroelasticSSM.msCer*x_AE(1:length(AeroelasticSSM.msCer));
dCl_AE=F_el(3)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
dCm_AE=F_el(5)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);

CM0=wingaero.CM+delta_CM+dCm_AE;


elev=elev+d_elev;
aircraft=aircraft.f_set_control_surface(trim_name,elev);
if ~isempty(trim_name2)
   aircraft=aircraft.f_set_control_surface(trim_name2,-elev); 
end
if ~isempty(trim_name3)
   aircraft=aircraft.f_set_control_surface(trim_name3,elev); 
end
if ~isempty(trim_name4)
   aircraft=aircraft.f_set_control_surface(trim_name4,elev); 
end



aircraft=aircraft.compute_grid();
if nargin>3
        aircraft=aircraft.compute_deflected_grid(def);
        wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);
else
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);
end
wingaero=wingaero.f_solve_for_Cl(flight_state.get_Cl(aircraft.reference.S_ref)-delta_Cl-dCl_AE);
if consider_thrust==1
    wingaero=wingaero.f_solve_full();
    T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)*cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
    delta_CM=0;
    delta_Cl=0;
    for i=1:length(aircraft.engines)
        aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end
x_AE=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 wingaero.state.alpha*pi/180  0   zeros(1,6)  zeros(1,6)]'-AeroelasticSSM.msASce*[cell2mat(aircraft.control_deflections)';cell2mat(aircraft.control_deflections)'*0;cell2mat(aircraft.control_deflections)'*0]*pi/180);
F_el=AeroelasticSSM.msCer*x_AE(1:length(AeroelasticSSM.msCer));
dCl_AE=F_el(3)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
dCm_AE=F_el(5)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);  
    
CM0p=wingaero.CM+delta_CM+dCm_AE;

dCMelev=(CM0p-CM0)/d_elev;
k=1;

while(abs(wingaero.CM+delta_CM+dCm_AE)>1E-5)
    fprintf('CM= %f, elev=%f, alpha=%f\n',wingaero.CM+delta_CM+dCm_AE,elev, wingaero.state.alpha);
    CM_prev=wingaero.CM+delta_CM+dCm_AE;
    d_elev=-CM_prev/dCMelev;
    elev=elev+d_elev;
    aircraft=aircraft.f_set_control_surface(trim_name,elev);
    if ~isempty(trim_name2)
        aircraft=aircraft.f_set_control_surface(trim_name2,-elev);
    end
    if ~isempty(trim_name3)
        aircraft=aircraft.f_set_control_surface(trim_name3,elev);
    end
    if ~isempty(trim_name4)
        aircraft=aircraft.f_set_control_surface(trim_name4,elev);
    end
    
    aircraft=aircraft.compute_grid();
    if nargin>3
        aircraft=aircraft.compute_deflected_grid(def);
        wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);
    else
        wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,flight_state.aerodynamic_state,trim_ref);
    end
    wingaero=wingaero.f_solve_for_Cl(flight_state.get_Cl(aircraft.reference.S_ref)-delta_Cl-dCl_AE);
    if consider_thrust==1
        wingaero=wingaero.f_solve_full();
        T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)*cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
        delta_CM=0;
        delta_Cl=0;
        for i=1:length(aircraft.engines)
            aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
            delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
            delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
        end
        
    else
        delta_CM=0;
        delta_Cl=0;
    end
    
    x_AE=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 wingaero.state.alpha*pi/180  0   zeros(1,6)  zeros(1,6)]'-AeroelasticSSM.msASce*[cell2mat(aircraft.control_deflections)';cell2mat(aircraft.control_deflections)'*0;cell2mat(aircraft.control_deflections)'*0]*pi/180);
F_el=AeroelasticSSM.msCer*x_AE(1:length(AeroelasticSSM.msCer));
dCl_AE=F_el(3)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
dCm_AE=F_el(5)/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
    
        
    dCMelev=(wingaero.CM+delta_CM+dCm_AE-CM_prev)/d_elev;
    k=k+1;
    if k>15
        fprintf('Trimming not converged\n');
    end
    flight_state.aerodynamic_state=wingaero.state;
end

%wingaero=wingaero.f_solve_full();

flight_state.aerodynamic_state=wingaero.state;
for i=1:length(flight_state.aircraft_state.control_surfaces)
    if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name)
        flight_state.aircraft_state.control_deflections{i}=elev;
    end
    if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name2) 
        flight_state.aircraft_state.control_deflections{i}=-elev;
    end
    
    if ~isempty(trim_name3)
         if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name3) 
        aircraft=aircraft.f_set_control_surface(trim_name3,elev);
         end
    end
    if ~isempty(trim_name4)
        if  strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name4)
            aircraft=aircraft.f_set_control_surface(trim_name4,elev);
        end
    end
end
  x_AE=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 wingaero.state.alpha*pi/180  0   zeros(1,6)  zeros(1,6)]'-AeroelasticSSM.msASce*[cell2mat(aircraft.control_deflections)';cell2mat(aircraft.control_deflections)'*0;cell2mat(aircraft.control_deflections)'*0]*pi/180);

flight_state.aircraft_state.control_deflections=aircraft.control_deflections;
 if consider_thrust==1
     flight_state.aircraft_state.engine_thrust=T_pe;
 end
end

