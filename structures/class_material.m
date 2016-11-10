%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_material.m
%> @brief File containing the class for material information
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
%> @brief class containing material properties
% ======================================================================


classdef class_material
    
    properties
        %> young's modulus
        E;
        %> shear modulus
        G;
        %> material density
        rho;
        %> allowable yield stress
        sigma_allowable;
        %> allowable shear stress
        tau_allowable;
    end
    
    methods
        
    % =================================================================
    %> @brief Class constructor
    %>
    %> Initializes the coupling condition
    %>
    %> @param material_string name of the material
    %>
    %> @return instance of class_material
    % =================================================================
    function obj =class_material(material_string)
        
        if strcmp('aluminum',material_string)
            %Al
            obj.E=6.89E10;
            obj.G=2.6E10;
            obj.rho=2800;
            obj.sigma_allowable=250E6;
            obj.tau_allowable=250E6;  
        end
        
        if strcmp('CFRP',material_string)
            %%TODO: to be defined correctly
            obj.E=5.5E10;
            obj.G=2E10;
            obj.rho=1600;
            obj.sigma_allowable=250E6;
            obj.tau_allowable=250E6;  
        end
        
        if strcmp('AGARD',material_string)
            %%TODO: to be defined correctly
            obj.E=6.0E9;
            obj.G=2.9*0.4392E9;
            obj.rho=381.98;
            obj.sigma_allowable=250E6;
            obj.tau_allowable=250E6;  
        end
    end
   
    end
end

