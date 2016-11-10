%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_flight_state
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        V
        
        %> Roll rate in body COS                          double      (1x1)   [rad/s]*
        p = 0    
        %> Pitch rate in body COS                         double      (1x1)   [rad/s]*
        q = 0 
        %> Yaw rate in body COS                           double      (1x1)   [rad/s]*  
        r = 0                       
        
        %> simulation starting time
        t0=0;
        %> time of the actual state
        ts=0;               
        %> aerodynamic state
        aerodynamic_state;
        %> 
        aircraft_state;
        %> load factor
        load_factor=1;
        %> altitude
        h=0;
        
        % design speeds
        M_D;
        M_C;
        VD; %design dive speed
        VC; % design cruise speed
        VC_EAS; %design cruise speed (equivalent airspeed)
        VD_EAS; % design dive speed (equivalent airspeed)
        loadcase_index;     % loadcase index of this state (for visualization)
        % body states
        
        reference;

    end
    
    methods
        function obj=class_flight_state(H,V,alpha,beta,aircraft_state)     
           obj.h=H;
           obj.V=V;
           [rho_air,a,T,P,mu]=stdatmo(H);
           Uinf=norm(V);
           Ma=Uinf/a;
           obj.aerodynamic_state=class_aero_state(Uinf,alpha,beta,Ma,rho_air,mu);
           obj.aircraft_state=aircraft_state;
        end
        
        function Cl=get_Cl(obj,S_ref)
            Cl=2*obj.aircraft_state.weight*9.81*obj.load_factor/(obj.aerodynamic_state.rho_air*norm(obj.V)^2*S_ref);
        end
        
        function obj=f_set_V_A(obj,V_A)
            obj.aerodynamic_state.V_A=V_A; 
            obj.aerodynamic_state.V_inf=[V_A*cosd(obj.aerodynamic_state.alpha)*cosd(obj.aerodynamic_state.beta)  V_A*sind(obj.aerodynamic_state.beta) V_A*sind(obj.aerodynamic_state.alpha)];
            %obj.V_inf=[V_A*cosd(obj.alpha)*cosd(obj.beta)  V_A*sind(obj.beta) V_A*sind(obj.alpha)];
        end
        
        
    end
    
end

