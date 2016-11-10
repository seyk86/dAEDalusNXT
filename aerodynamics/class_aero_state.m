%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus aerodynamics
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se


classdef class_aero_state
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % aerodynamic state variables
        alpha = 0;                  % Angle of attack                                double      (1x1)   [rad]  *
        beta = 0;                   % Angle of sideslip                              double      (1x1)   [rad]  *
        V_inf=0;                    % Freestream Vector                             double      (3x1)   [m/s]  * 
        V_A=10;                   	% Aerodynamic velocity                           double      (1x1)   [m/s]  * 
        mu=1.81E-5;                 % viscosity kg/(m*s)
        Ma;
        rho_air;
        p_ref=[0 0 0];
    end
    
    methods
        
        function obj=class_aero_state(Uinf,alpha,beta,Ma,rho_air,varargin)
            obj.V_A=Uinf;
            obj.V_inf=[Uinf*cosd(alpha)*cosd(beta)  Uinf*sind(beta)*cosd(alpha) Uinf*sind(alpha)];
            obj.alpha=alpha;
            obj.beta=beta;
            obj.Ma=Ma;
            obj.rho_air=rho_air;
            
            obj.Ma=Ma;
            
            if nargin==6
               obj.mu=varargin{1};   
            end
        end
        
        function obj=set_alpha(obj,alpha)
            obj.alpha=alpha;
            obj.V_inf=[obj.V_A*cosd(obj.alpha)*cosd(obj.beta)  obj.V_A*sind(obj.beta) obj.V_A*sind(obj.alpha)];
        end
        
        function Vinf=velocity_vec(obj)
           Vinf=obj.V_A; 
           
        end
    end
end
