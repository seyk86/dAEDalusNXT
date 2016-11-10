%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       This file is part of dAEDalus criticaldesign
%                       Copyright (C) 2011, Klaus Seywald
%     Author:           Klaus Seywald
%                       klaus.seywald@mytum.de
%                       seywald@kth.se

function [ailstate] =critical_aileron_roll_load_state(ref_state,gfactor,wingno,angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    ailstate=ref_state;
    ailstate.load_factor=gfactor;
    ailstate.aerodynamic_state.V_A=ref_state.VD/1.15;
    ailstate.aircraft_state.control_deflections{find(not(cellfun('isempty',strfind(ref_state.aircraft_state.control_surfaces,[wingno, '_left']))))}=angle;
    ailstate.aircraft_state.control_deflections{find(not(cellfun('isempty',strfind(ref_state.aircraft_state.control_surfaces,[wingno, '_right']))))}=angle;
end

