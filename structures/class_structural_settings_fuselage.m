%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_structural_settings_fuselage
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        n_partitions;                           % number of fuselage partitions
        
        t_min_sk=0.001;%1.27e-3;                       % fuselage minimum skin gauge
        
        fueled_span=0;                            % fuselage fueled length in % of length
        fuel_density=807.5;                     % fuel density (Jet-A1)
        
        delta_pressure=0;                       % cabin differential pressure
       
        materials;
    end
    
    methods
    	function obj = f_set_material(obj,material)
              obj.materials=material;
        end        
    end
end
