%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_crosssection_wingbox.m
%> @brief File containing the class for a finite wingbox crosssection
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
%> @brief class for a finite fuselage crosssection
%>
%> contains all GEOMETRICAL and STRUCTURAL information for a finite fuselage crosssection
% ======================================================================

classdef class_crosssection_fuselage
    
    
    properties

        %> fuselage radius                                      [m]
        r;
        %> wetted crosssection in fuselage segment              [m2]
        A_fuel;
        
        %% environmental factors     
        
        %> cabin differential pressure                          [Pa]
        delta_pressure=0;
        
        %% simple fuselage crosssection geometry
        %% Skin
        %> equivalent skin thickness (accounts stringers)       [m]
        t_sk_eq=0;  
        %> loadcase index for the equivalent skin               [m]
        t_sk_eq_lc_idx=0;
               
        %> equivalent frame thickness                           [m]
        t_fr_eq=0;  
        %> loadcase index for the equivalent skin               [m]
        t_fr_eq_lc_idx=0;
        
        %> minimum allowable skin thickness
        t_min_sk;
        
        %% material properties at cross-section
        
        %> youngs modulus of the equivalent skin
        E_sk_eq;
        %> shear modulus of the equivalent skin
        G_sk_eq; 
        %> density of the equivalent skin
        rho_sk_eq;
        %> youngs modulus of the frame
        E_fr_eq;
        %> shear modulus of the frame
        G_fr_eq;
        %> density of the frame
        rho_fr_eq;
        
        %> tensile yield strength of the equivalent skin
        tensile_yield_sk;
        %> tensile yield strength of the frame
        tensile_yield_fr;
        %> safety factor
        safety_factor=1.5;
        %> differential pressure safety factor
        safety_factor_pressure=1.33;
        
        
        filling_ratio=1; % ratio of max payload in element

    end
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> @return instance of the class_crosssection_fuselage
        % =================================================================
         function obj = class_crosssection_fuselage()
                obj.r=0;
                obj.A_fuel=0;  
         end
        % =================================================================
        %> @brief set geometry
        %>
        %> Sets the geometry parameters of the fuselage crosssection
        %>
        %> @param r fuselage crosssection radius
        %>
        %> @return instance of the class_crosssection_fuselage
        % =================================================================
        function obj=setGeometry(obj,r,varargin)
            obj.r=r;
            if nargin==3
            obj.t_sk_eq=varargin{1};
            end
        end
        
        function obj=setMaterial(obj,E_sk_eq,G_sk_eq,rho_sk_eq,E_fr_eq,G_fr_eq,rho_fr_eq,Ult_Tstrength_sk,Ult_Tstrength_fr,t_min_sk)
            obj.E_sk_eq=E_sk_eq;
            obj.G_sk_eq=G_sk_eq; 
            obj.rho_sk_eq=rho_sk_eq;
            %> Frames
            obj.E_fr_eq=E_fr_eq;
            obj.G_fr_eq=G_fr_eq;
            obj.rho_fr_eq=rho_fr_eq;
            
%             obj.tensile_yield_sk=0.66*Ult_Tstrength_sk;
            obj.tensile_yield_sk=Ult_Tstrength_sk;
            obj.tensile_yield_fr=Ult_Tstrength_fr;
            
            obj.t_min_sk=t_min_sk;
        end
        
        function [r]=get_dimensions(obj)
            r=obj.r;
        end
        
        
        % =================================================================
        %> @brief f_calc_dm
        %>
        %> calculates the mass dm of a fuselage crosssection
        %>
        %> @return instance of the class_crosssection_fuselage
        % =================================================================
        function dm=f_calc_dm(obj)
            dm=2*pi*obj.r*(obj.t_sk_eq*obj.rho_sk_eq+obj.t_fr_eq*obj.rho_fr_eq);
        end
        
        % =================================================================
        %> @brief f_self_design_crosssection
        %>
        %> function to perform the crosssectional self-design ( this
        %> function contains the sizing rules for the crosssection)
        %>
        %> @return instance of the class_crosssection_fuselage
        % =================================================================
        function obj=f_self_design_crosssection(obj,Mbx,Mbz,Mt,Qx,Qz,loadcase_idx)
            
            % Skin sizing
            % Pressure Case
            t_p(1)=obj.safety_factor*obj.safety_factor_pressure*obj.delta_pressure*obj.r/obj.tensile_yield_sk;
            t_p(2)=0.5*t_p(1);
            
            % Kp accounts that not all the shell material resistes hoop
            % stress
            Kp=1.576;
            t_p(1)=Kp*t_p(1);
                        
            % Bending Moment Case
            t_b=abs(obj.safety_factor*Mbx)/(pi*obj.tensile_yield_sk*obj.r^2);
            
            % Shear Force Case
            t_s=abs(obj.safety_factor*Qz)/(pi*obj.r*0.55*obj.tensile_yield_sk);
            
            % Torsional Moment Case
            t_t=abs(obj.safety_factor*Mt)/(2*pi*0.55*obj.tensile_yield_sk*obj.r^2);
            
            % Parameter depending on the shell geometry that relates the equivalent
            % thickness to the specified minimum material thickness
            k_mg=1.8352;
%             k_mg=1;

%             t=max([t_p(1) t_p(2)+t_b k_mg*obj.t_min_sk]);
            
            if Mbx<0 && t_b>t_p(2)
                t=max([t_p(1) t_b k_mg*obj.t_min_sk]); % not overestimate
            else
                t=max([t_p(1) t_p(2)+t_b k_mg*obj.t_min_sk]);
            end
   
            obj.t_sk_eq=t;
            obj.t_sk_eq_lc_idx=loadcase_idx;
           
            % Frames sizing
            obj.t_fr_eq=obj.t_sk_eq/3;
            obj.t_fr_eq_lc_idx=loadcase_idx;
            
        end
        
        % =================================================================
        %> @brief calc_crosssection
        %>
        %> calculates the crosssectional parameterx Ix, Iy, Iy, A,
        %> Aenclosed
        %>
        %> @return Ix
        %> @return Iy
        %> @return Iz
        %> @return A
        %> @return Aenclosed
        % =================================================================
        function [Ix,Iy,Iz,J,A,Aenclosed]=calc_crosssection(obj)
            Iz=pi*obj.t_sk_eq*obj.r^3;
            Ix=pi*obj.t_sk_eq*obj.r^3;
            Iy=Ix+Iz;
            J=pi*0.5*((obj.r+obj.t_sk_eq)^4-(obj.r-obj.t_sk_eq)^4);
            A=2*pi*obj.t_sk_eq*obj.r;        % material crosssectional area
            Aenclosed=(pi*obj.r^2)-A;        % CHECK! Only centre wingbox fuselage has fuel % used to compute fuel volume.
        end
    end   
end

