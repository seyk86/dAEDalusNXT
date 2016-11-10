%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,flight_state,wingaero ] = trim_aircraft_unsteady(aircraft,flight_state,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dt=0.002;
wingaero.t_step=dt;

trim_name3=[];
trim_name4=[];
consider_thrust=0;
for i=1:length(flight_state.aircraft_state.control_deflections)
    %if (flight_state.aircraft_state.control_deflections{i}>=1E-4)
        aircraft=aircraft.f_set_control_surface(flight_state.aircraft_state.control_surfaces{i},flight_state.aircraft_state.control_deflections{i});
   % end
end

aircraft.grid_settings.wake=2;
def=0;
elev=0;
d_elev=2;
no_redeform=0;
trim_name=[];
trim_name2=[];
trim_name3=[];
trim_name4=[];
aircraft.reference.p_ref=flight_state.aircraft_state.CG_ref;
if nargin==4
    wingstructure=varargin{1};
    def=wingstructure.f_get_deflections;
end
if nargin==6
    wingstructure=varargin{1};
    elev=0;
    def=varargin{4};
end
trim_idx=1;

if isempty(aircraft.control.trim_surfaces)
    trim_name='elevator';
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
else
    trim_name=aircraft.control.trim_surfaces{1};
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
end

elev=aircraft.control_deflections{trim_idx};
if length(aircraft.control.trim_surfaces)==2
    trim_name2=aircraft.control.trim_surfaces{2};
end

if length(aircraft.control.trim_surfaces)==4
    trim_name2=aircraft.control.trim_surfaces{2};
    trim_name3=aircraft.control.trim_surfaces{3};
    trim_name4=aircraft.control.trim_surfaces{4};
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
    UVLM_settings=class_UVLM_computation_settings();
    UVLM_settings.debug=0;
    UVLM_settings.movie=0;
    wingaero=class_UVLM_solver(aircraft.name,aircraft.grid_deflected,aircraft.is_te,aircraft.panels,flight_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
    wingaero=wingaero.initialize_time_domain_solution(dt);
else
    UVLM_settings=class_UVLM_computation_settings();
    UVLM_settings.debug=0;
    UVLM_settings.movie=0;
    aircraft.reference.p_ref=-mean_axis_origin';
    wingaero=class_UVLM_solver(aircraft.name,aircraft.grid,aircraft.is_te,aircraft.panels,flight_state.aerodynamic_state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);
    wingaero=wingaero.initialize_time_domain_solution(dt);
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
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end
CM0=wingaero.CM+delta_CM;


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
    
    wingaero.grid(1,:)=aircraft.grid_deflected(1,:)/wingaero.Ma_corr;
    wingaero.grid(2,:)=aircraft.grid_deflected(2,:);
    wingaero.grid(3,:)=aircraft.grid_deflected(3,:);
    
    wingaero=wingaero.update_grid();
else
    wingaero.grid(1,:)=aircraft.grid(1,:)/wingaero.Ma_corr;
    wingaero.grid(2,:)=aircraft.grid(2,:);
    wingaero.grid(3,:)=aircraft.grid(3,:);
    wingaero=wingaero.update_grid();
end
wingaero=wingaero.f_solve_for_Cl(flight_state.get_Cl(aircraft.reference.S_ref)-delta_Cl);
if consider_thrust==1
    wingaero=wingaero.f_solve_full();
    T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)*cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
    delta_CM=0;
    delta_Cl=0;
    for i=1:length(aircraft.engines)
        aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end
    
    
CM0p=wingaero.CM+delta_CM;

dCMelev=(CM0p-CM0)/d_elev;
k=1;

while(abs(wingaero.CM+delta_CM)>1E-5)
    fprintf('CM= %f, elev=%f\n',wingaero.CM+delta_CM,elev);
    CM_prev=wingaero.CM+delta_CM;
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
        
        wingaero.grid(1,:)=aircraft.grid_deflected(1,:)/wingaero.Ma_corr;
        wingaero.grid(2,:)=aircraft.grid_deflected(2,:);
        wingaero.grid(3,:)=aircraft.grid_deflected(3,:);
        
        wingaero=wingaero.update_grid();
    else
        wingaero.grid(1,:)=aircraft.grid(1,:)/wingaero.Ma_corr;
        wingaero.grid(2,:)=aircraft.grid(2,:);
        wingaero.grid(3,:)=aircraft.grid(3,:);
        wingaero=wingaero.update_grid();
    end
    n_wake=wingaero.settings.wakelength_factor*wingaero.reference.b_ref/(wingaero.t_step*norm(wingaero.Uinf));
    wingaero=wingaero.initialize_vring_grid(wingaero.t_step);
    wingaero=wingaero.initialize_wake_grid(ceil(n_wake),wingaero.t_step);
%     wingaero=wingaero.compute_influence_coeff_matrix();
    wingaero=wingaero.f_solve_for_Cl(flight_state.get_Cl(aircraft.reference.S_ref)-delta_Cl);
    if consider_thrust==1
        wingaero=wingaero.f_solve_full();
        T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)*cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
        delta_CM=0;
        delta_Cl=0;
        for i=1:length(aircraft.engines)
            aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
            delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
            delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
        end
        
    else
        delta_CM=0;
        delta_Cl=0;
    end
        
    dCMelev=(wingaero.CM+delta_CM-CM_prev)/d_elev;
    k=k+1;
    if k>10
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

flight_state.aircraft_state.control_deflections=aircraft.control_deflections;
end

