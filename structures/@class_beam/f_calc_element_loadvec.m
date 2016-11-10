%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [obj] =f_calc_element_loadvec(obj,add_eigenmass,add_fuelmass)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    el_ndof=obj.el_ndof;
    %% calculate element stiffnes matrices and load vectors
    if el_ndof==3
        for i=1:obj.nel  
            %Calculate Element Load Vector
            obj.beamelement(i)=obj.beamelement(i).elq_3dof();
        end
    elseif el_ndof==6
         for i=1:obj.nel      
            %Calculate Element Load Vector
            obj.beamelement(i)=obj.beamelement(i).elq_6dof(add_eigenmass,add_fuelmass);
         end
    else
        % number of ndof not implemented
    end

end

