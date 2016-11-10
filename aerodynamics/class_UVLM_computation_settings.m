%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_UVLM_computation_settings
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coeff_eps;
        % number of timesteps/wavelength
        spp;
        % waklength in multiples of b_ref
        wakelength_factor;
        % number of harmonic oscillations
        n_osc;
        
        % select debug 0/1
        debug;
        
        % variables for animation
        movie;
        filename=[];
        
        modal_data;
        n_mode;
    end
    
    methods
        
        function obj=class_UVLM_computation_settings(obj)
            obj.coeff_eps=1e-10;
            % number of timesteps/wavelength
            obj.spp=32;
            % waklength in multiples of b_ref
            obj.wakelength_factor=1;
            % number of harmonic oscillations
            obj.n_osc=10;
            
            % select debug 0/1
            obj.debug=0;
            
            % variables for animation
            obj.movie=0;
            obj.filename=[];
            
            obj.modal_data=0;
            
        end
        
        
    end
    
end

