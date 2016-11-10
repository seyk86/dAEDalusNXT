%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_engine
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        cg_pos;
        
        m;
        
        mass_engine;
        mass_pylon;
        delta_t=1;                % thrust lever postion
        
        thrust;                 % thrust
        thrust_vec=[1 0 0];     % thrust vector ( standard [1 0 0] )
        mounting='wing_mounted';
        
    end
    
    methods
    end
    
end

