%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [gmaneuverstate] =critical_sideslip_maneuver_state(ref_state,beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    gmaneuverstate=ref_state;
    gmaneuverstate.aerodynamic_state.V_A=ref_state.VC;
    gmaneuverstate.load_factor=1;
    gmaneuverstate.aerodynamic_state.beta=beta; 
    gmaneuverstate.aerodynamic_state.V_inf=[ref_state.VC*cosd(beta) ref_state.VC*sind(beta) 0];
end

