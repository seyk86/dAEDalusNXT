%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus structures
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se

classdef class_structural_settings_wing
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
%         n_partitions;                           % number of wing partitions
%         
%         t_min_sk=1.5e-3;                          % minimum sheet thickness for skin
%         t_min_sp=1.5e-3;                          % minimum thickness for spars
% %         
        t_min_sk=0.0005;                          % minimum sheet thickness for skin
        t_min_sp=0.0005;                          % minimum thickness for spars
        
        p_sp_fr;                                % front spar positions ( in % of chord )
        p_sp_re;                                % rear spar positions ( in % of chord )
        
        fueled_span;                            % wing fueled span in % of span
        fuel_density=807.5;                     % fuel density (Jet-A1)
        
        materials;
    end
    
    methods
    	function obj = f_set_material(obj,material)
              obj.materials=material;
        end        
    end
end

