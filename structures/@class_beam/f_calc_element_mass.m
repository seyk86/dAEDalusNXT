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
function [obj] = f_calc_element_mass(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    el_ndof=obj.el_ndof;
    %% calculate element stiffnes matrices and load vectors
    if el_ndof==6
         for i=1:obj.nel   
            %Calculate Element mass Matrices
            obj.beamelement(i)=obj.beamelement(i).elM_6dof();    
            obj.beamelement(i)=obj.beamelement(i).elM_lumped_6dof(); 
         end
    else
        % number of ndof not implemented
    end
end
