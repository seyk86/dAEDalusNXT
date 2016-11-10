%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_boundary_condition.m
%> @brief File containing the class for structural boundary conditions
%>        This file is part of dAEDalus structures, Copyright (C) 2011,
%>        Klaus Seywald
%>   dAEDalus is published under the terms of the GNU General Public
%>   License by the Free Software Foundation;
%>   either version 2 or any later version.
%>
%>   You should have received a copy of the GNU General Public
%>   License along with dAEDalus; see the file GNU GENERAL 
%>   PUBLIC LICENSE.TXT.  If not, write to the Free Software 
%>   Foundation, 59 Temple Place -Suite 330, Boston, MA
%>   02111-1307, USA.
%>
% ======================================================================
%> @brief class for beam boundary conditition
%> This class is a container for structural boundary conditions
% ======================================================================

classdef class_boundary_condition
    properties
        %> index of node where boundary condition should be set
        node_idx;      
        %> define what DOF should be set for example [1 1 1 0 0 0] 
        dof;            
        %> define prescribed values [0 0 0 0 0 0]
        prescribed_val;     
    end
    
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> Initializes the boundary condition
        %>
        %> @param node_idx Node at Boundary Condition
        %> @param dof boolean vector setting DOF or not: e.g. [1 1 1 0 0 0] 
        %> @param prescribed_val define prescribed values [0 0 0 0 0 0]
        %>
        %> @return instance of the class_boundary_condition class
        % =================================================================
          function obj = class_boundary_condition(node_idx,dof,prescribed_val)
              obj.node_idx=node_idx;
              obj.dof=dof;
              obj.prescribed_val=prescribed_val;
          end
    end
end

