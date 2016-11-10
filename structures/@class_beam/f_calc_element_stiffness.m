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
function [obj] = f_calc_element_stiffness(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    el_ndof=obj.el_ndof;
    %% calculate element stiffnes matrices and load vectors
    if el_ndof==3
        for i=1:obj.nel  
            %Calculate Element Load Vector
            obj.beamelement(i)=obj.beamelement(i).lin_elK_3dof();
        end
    elseif el_ndof==6
         for i=1:obj.nel   
            %Calculate Element Stiffness Matrices
            obj.beamelement(i)=obj.beamelement(i).lin_elK_6dof();      
            
%             %---- Simon: enforce beam property symmetry
%             if i>obj.nel/2
%                obj.beamelement(i)=obj.beamelement(i-((i-(obj.nel/2)-1)*2+1)).lin_elK_6dof();
%             else
%                obj.beamelement(i)=obj.beamelement(i).lin_elK_6dof();      
%             end
%             %----
         end
    else
        % number of ndof not implemented
    end
end

