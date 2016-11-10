%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ ctrl_input,t_vec,maneuver_name] = generate_input_signal(aircraft,cs_name,name,type,step_length,tstart,tend,dt,amplitude)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t_vec=0:dt:tend;
maneuver_name=name;

if strcmp(type,'doublet')
    ctrl_input.time=[0 tstart tstart+dt tstart+step_length/2-dt tstart+step_length/2+dt tstart+step_length-dt tstart+step_length tend];
    ctrl_input.signals.values=zeros(length(ctrl_input.time),length(aircraft.control_surfaces));
    ctrl_input.signals.dimensions=length(aircraft.control_surfaces);
    for cs=1:length(cs_name)
        for i=1:length(aircraft.control_surfaces)
            if strcmp(aircraft.control_surfaces{i},cs_name{cs})
                ctrl_input.signals.values(:,i)=[0 0 -1 -1 1 1 0  0]*amplitude(cs);
            end
        end
    end
elseif strcmp(type,'singlet')
    ctrl_input.time=[0 tstart tstart+dt tstart+step_length/2-dt tstart+step_length/2+dt tstart+step_length-dt tstart+step_length tend];
    ctrl_input.signals.values=zeros(length(ctrl_input.time),length(aircraft.control_surfaces));
    ctrl_input.signals.dimensions=length(aircraft.control_surfaces);
    for cs=1:length(cs_name)
        for i=1:length(aircraft.control_surfaces)
            if strcmp(aircraft.control_surfaces{i},cs_name{cs})
                ctrl_input.signals.values(:,i)=[0 0 1 1 1 1 0  0]*amplitude(cs);
            end
        end
    end
    
elseif strcmp(type,'sine')
    ctrl_input.time=[0 tstart tstart+[0:dt:step_length]+dt tstart+step_length+2*dt tend];
    ctrl_input.signals.values=zeros(length(ctrl_input.time),length(aircraft.control_surfaces));
    ctrl_input.signals.dimensions=length(aircraft.control_surfaces);
    ctrl_input.signals.values(1:length(ctrl_input.time),1:ctrl_input.signals.dimensions)=0;
    for cs=1:length(cs_name)
        for i=1:length(aircraft.control_surfaces)
            if strcmp(aircraft.control_surfaces{i},cs_name{cs})
                ctrl_input.signals.values(:,i)=[0 0 amplitude(cs)*sin(2*pi/step_length*[0:dt:step_length]) 0  0];
            end
        end
    end
end
end
