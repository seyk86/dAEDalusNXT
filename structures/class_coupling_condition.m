%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_coupling_condition.m
%> @brief File containing the class for structural coupling conditions
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
%> @brief class for coupling conditions of beam systems
%> This class is a container for coupling conditions of beam systems
% ======================================================================

classdef class_coupling_condition
    properties
        %> index of beam where coupling should be set, example:[1 2] couple beam 1 and 2
        beam_idx;  
        %> index of node where coupling should be set, example:[40 1] couple node 40 of beam 1 to node 1 of beam 2
        node_idx;
        %> define what DOF should be coupled for example [1 1 1 0 0 0] 
        dof;        
    end
    
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> Initializes the coupling condition
        %>
        %> @param node_idx vector with nodes to be coupled e.g. [51 32]
        %> @param beam_idx vector with beam indices to be coupled e.g.[2 3]
        %> @param dof define what DOF should be coupled [0 0 0 0 0 0]
        %>
        %> @return instance of the class_boundary_condition class
        % =================================================================
         function obj = class_coupling_condition(node_idx,beam_idx,dof)
              obj.node_idx=node_idx;
              obj.beam_idx=beam_idx;
              obj.dof=dof;
         end
    end
end

