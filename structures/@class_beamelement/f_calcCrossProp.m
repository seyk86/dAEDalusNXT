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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_calcCrossProp(obj)
%% careful.. only valid if front and rear spar have equal thickness

    
    if isa(obj.crosssection,'class_crosssection_wingbox')
         [obj.Ix,obj.Ip,obj.Iz,obj.J,obj.A,obj.A_enclosed]=obj.crosssection.calc_crosssection();
    elseif isa(obj.crosssection,'class_crosssection_wingbox_optimized_shape')
         [obj.Ix,obj.Ip,obj.Izx,obj.Iz,obj.J,obj.A,obj.A_enclosed]=obj.crosssection.calc_crosssection();
    elseif isa(obj.crosssection,'class_crosssection_C_shapeV2')
        [obj.Ix,obj.Ip,obj.Izx,obj.Iz,obj.J,obj.A,obj.A_enclosed]=obj.crosssection.calc_crosssection();
    elseif isa(obj.crosssection,'class_crosssection_fuselage')
        [obj.Ix,obj.J,obj.Iz,obj.Ip,obj.A,obj.A_enclosed]=obj.crosssection.calc_crosssection();
    end
end
