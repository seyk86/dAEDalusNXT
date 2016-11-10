%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [gmaneuverstate] =critical_g_maneuver_state(ref_state,gfactor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    gmaneuverstate=ref_state;
    gmaneuverstate.V=[ref_state.VD 0 0];
    gmaneuverstate.aerodynamic_state.V_A=ref_state.VD;
    gmaneuverstate.load_factor=gfactor;
end

