%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_wing<class_aerosurface
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K;
        K_B;
        Cl;
        Cl_min;
    end
    
    methods
        
        function obj = class_wing(pos,symmetric,dihed,LambdaSpec,Lambda,varargin)     
                obj=obj@class_aerosurface(pos,symmetric,dihed,LambdaSpec,Lambda,varargin);
        end
        
    end
    
end

