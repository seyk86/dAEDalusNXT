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

function [obj ] = f_add_nodal_loads(obj)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
%% dummy for follower loads (non linear FEM)
%     for i=1:1:length(obj.beamelement)+1
%         a=obj.nodal_deflections(4+(i-1)*6);
%             Lx=[1       0       0
%             0   cos(a)  sin(a)
%             0   -sin(a) cos(a)];
%         obj.nodal_loads_def(1+(i-1)*6:3+(i-1)*6)=Lx*obj.nodal_loads(1+(i-1)*6:3+(i-1)*6);  
%     end
    
%     obj.nodal_deflections
%     obj.nodal_loads


    for i=1:1:length(obj.Q)
        obj.Q(i)=obj.Q(i)+obj.nodal_loads(i);
    end
end

