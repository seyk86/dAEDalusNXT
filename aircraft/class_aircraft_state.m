%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus aeroelasticity
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se

classdef class_aircraft_state
    properties
        %> number of control surfaces
        n_ctrl;
        %> control surfaces names
        control_surfaces;
        %> control surface deflections degrees
        control_deflections; 
        %> trim weight (lift=weight)
        weight;      
        %> aircraft fueling state
        fuel_state;         % fueling state of a/c
        %> aircraft payload state
        payload_state; % cargo state of a/c
       
        engine_thrust;
        
        CG_ref = [0;0;0]       % Position of center of gravity    double  (3x1)  [m]
        
        reference;
        
        I_xyz;
    end
    
    methods    
        function obj=class_aircraft_state(aircraft,CG,varargin)
            
            obj.n_ctrl=length(aircraft.control_surfaces);
            obj.control_surfaces=aircraft.control_surfaces;
            obj.control_deflections=aircraft.control_deflections;
            
            obj.CG_ref=CG;
            
            if nargin==2
            obj.weight=aircraft.weights.MTOW;
            elseif nargin==3
                obj.weight=varargin{1};
            end
            
                R=[0.246,0.382,0.456];
            b=aircraft.reference.b_ref;
            d=37;
            OWE=obj.weight;
            obj.I_xyz=[(OWE/32.174)*(R(1)*b/2*3.28)^2 ,0,0;
            0,(OWE/32.174)*(R(2)*d/2*3.28)^2,0;
            0,0,(OWE/32.174)*(R(3)*(b+d)/2*3.28/2)^2]*1.355817962;
        end
       
        function obj = set_rudder_deflection(obj,wingno,partition_no,deflection)
            obj.aileron_deflection(1).wingno=wingno;
            obj.aileron_deflection(1).partition_no=partition_no;
            obj.aileron_deflection(1).deflection=deflection;
        end   
    end   
end

