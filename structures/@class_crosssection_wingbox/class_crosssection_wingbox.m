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
%> @brief class for a finite wingbox crosssection
%>
%> contains all GEOMETRICAL and STRUCTURAL information for a finite wingbox crosssection
% ======================================================================

classdef class_crosssection_wingbox
    
    
    properties

        %> wingbox height                                       [m]
        h;
        %> wingbox width                                        [m]
        w;
        %> chordlength of wingsection                           [m]
        c; 
        %> enclosed area                                        [m2]   
        A_enclosed; 
        %> wetted crosssection in wing segment                  [m2]
        A_fuel; 
        
        %% simple wingbox geometry
        %% Skin'
        %> thickness of upper skin                              [m]
        t_sk_up;  
        %> loadcase index for the upper skin                    [m]
        t_sk_up_lc_idx=0;
        %> thickness of lower skin                              [m]
        t_sk_lo;    
        %> loadcase index for the lower skin                    [m]
        t_sk_lo_lc_idx=0;

        %> thickness of upper stringers                         [m]
        t_st_up;    
        %> thickness of lower stringers                         [m]
        t_st_lo;     
        
        %> front spar thickness                                 [m]
        t_sp_fr;
        %> rear spar thickness                                  [m]
        t_sp_re;
        %> loadcase index for the front spar                    [m]
        t_sp_fr_lc_idx=0;    
        %> loadcase index for the rear spar                     [m]
        t_sp_re_lc_idx=0; 
        
        %> minimum allowable spar thickness
        t_min_sp;
        %> minimum allowable skin thickness
        t_min_sk;
        
        %% material properties at cross-section
        
        %> youngs modulus of the upper skin
        E_sk_up;
        %> shear modulus of the upper skin
        G_sk_up; 
        %> density of the upper skin
        rho_sk_up;
        %> youngs modulus of the lower skin
        E_sk_lo;
        %> shear modulus of the lower skin
        G_sk_lo;
        %> density of the lower skin
        rho_sk_lo;
        %> density of the spars
        rho_sp;
        
        %> tensile yield strength of the spars
        tensile_yield_sp;
        %> tensile yield strength of the lower skin
        tensile_yield_sk_l;
        %> tensile yield strength of the upper skin
        tensile_yield_sk_u;
        %> safety factor
        safety_factor=1.5;
         fueling_factor=1;
    end
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
         function obj = class_crosssection_wingbox()
                obj.h=0;
                obj.w=0; 
                obj.c=0; 
                obj.A_enclosed=0; 
                obj.A_fuel=0;
                obj.t_sp_fr=0;    
                obj.t_sp_re=0;    
                obj.t_sk_up=0;   
                obj.t_sk_lo=0;    
         end
        % =================================================================
        %> @brief set geometry
        %>
        %> Sets the geometry parameters of the wingbox crosssection
        %>
        %> @param c wing crosssection chordlength
        %> @param h wingbox height
        %> @param w wingbox width
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj=setGeometry(obj,c,h,w)
            obj.c=c;
            obj.h=h;
            obj.w=w; 
            obj.A_enclosed=h*w;
        end
  
        
        function obj=setMaterial(obj,E_sk_up,G_sk_up,rho_sk_up,E_sk_lo,G_sk_lo,rho_sk_lo,rho_sp,Ult_Tstrength_sp,Ult_Tstrength_sk_u,Ult_Tstrength_sk_l,t_min_sp,t_min_sk)
            obj.E_sk_up=E_sk_up;
            obj.G_sk_up=G_sk_up; 
            obj.rho_sk_up=rho_sk_up;
            %b) Lower Skin
            obj.E_sk_lo=E_sk_lo;
            obj.G_sk_lo=G_sk_lo;
            obj.rho_sk_lo=rho_sk_lo;
            %d) spars
            obj.rho_sp=rho_sp;
            obj.tensile_yield_sp=Ult_Tstrength_sp;
            obj.tensile_yield_sk_u=Ult_Tstrength_sk_u;
            obj.tensile_yield_sk_l=Ult_Tstrength_sk_l;
            
            obj.t_min_sk=t_min_sk;
            obj.t_min_sp=t_min_sp;
        end
        
        function [h,w, c]=get_dimensions(obj)
                h=obj.h;
                w=obj.w;
                c=obj.c;
        end
        
        
        % =================================================================
        %> @brief f_calc_dm
        %>
        %> calculates the mass dm of a wingbox crosssection
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function dm=f_calc_dm(obj)
            dm=obj.rho_sk_up*obj.t_sk_up*obj.w + obj.rho_sk_lo*obj.t_sk_lo*obj.w+(obj.t_sp_fr+obj.t_sp_re)*obj.h*obj.rho_sp;
        end
        
        % =================================================================
        %> @brief f_self_design_crosssection
        %>
        %> function to perform the crosssectional self-design ( this
        %> function contains the sizing rues for the crosssection)
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj=f_self_design_crosssection(obj,Mbx,Mby,Mt,Qx,Qz,loadcase_idx,overwrite)
        
            nt=0.8; %bending efficiency
            
            sigma_lim=obj.tensile_yield_sk_l;
            
            t_b(1) = (abs(obj.safety_factor*Mbx)/(nt*obj.h*obj.w*sigma_lim));  % Equivalent upper skin thickness including stringer and skin thickness
            t_b(2) = (abs(obj.safety_factor*Mbx)/(nt*obj.h*obj.w*sigma_lim));  % Equivalent lower skin thickness including stringer and skin thickness

            t(1)= max(t_b(1),obj.t_min_sk);
            t(2)= max(t_b(2),obj.t_min_sk);
          
            % The absolute upper and lower wing skin thickness
            
            if ~overwrite
                if t(1)>obj.t_sk_up
                    obj.t_sk_up=t(1);
                    obj.t_sk_up_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_up=t(1);
                obj.t_sk_up_lc_idx=loadcase_idx;
            end
            

            if ~overwrite
                if t(2)>obj.t_sk_lo
                    obj.t_sk_lo=t(2);
                    obj.t_sk_lo_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_lo=t(2);
                obj.t_sk_lo_lc_idx=loadcase_idx;
            end
            

       
        
            % Sizing webs
            tQz = 3/2*abs(Qz)/(2*obj.h);                    % Shear flow in the webs due to the shear force
            tMt = abs(Mt)/(2*obj.h*obj.w);                  % Shear flow in the webs due to the torsional moments
            tWeb = tQz+tMt;                                 % Total shear flow
            
            tQx = abs(Qx)/(2*obj.w);                        % Shear flow in the skin due to the shear force
            tMt = abs(Mt)/(2*obj.h*obj.w);                  % Shear flow in the skin due to the torsional moments
            tSkin=tQx+tMt;                                  

            shear_ult=sigma_lim*0.55;                       % limit shear stress
                                                            % load better value
                                                            % from material
                                                            % database later
            t_w = obj.safety_factor*tWeb/shear_ult;         % Spar thickness at each node
            t_w = max(t_w,obj.t_min_sp);
        
            t_w_skin=obj.safety_factor*tSkin/shear_ult;    
            
            if(t_w_skin>t(1))
                
                if ~overwrite
                    if t_w_skin>obj.t_sk_up
                        obj.t_sk_up=t_w_skin;
                        obj.t_sk_up_lc_idx=loadcase_idx;
                    end
                else
                    obj.t_sk_up=t_w_skin;
                    obj.t_sk_up_lc_idx=loadcase_idx;
                end
                
                if ~overwrite
                    if t_w_skin>obj.t_sk_up
                        obj.t_sk_lo=t_w_skin;
                        obj.t_sk_lo_lc_idx=loadcase_idx;
                    end
                else
                    obj.t_sk_lo=t_w_skin;
                    obj.t_sk_lo_lc_idx=loadcase_idx;
                end
            end
               
            if ~overwrite
                if t_w>obj.t_sp_fr
                    obj.t_sp_fr=t_w;
                    obj.t_sp_fr_lc_idx=loadcase_idx;
                end
            else
                obj.t_sp_fr=t_w;
                obj.t_sp_fr_lc_idx=loadcase_idx;
            end
            
            if ~overwrite
                if t_w>obj.t_sp_re
                    obj.t_sp_re=t_w;
                    obj.t_sp_re_lc_idx=loadcase_idx;
                end
            else
                obj.t_sp_re=t_w;
                obj.t_sp_re_lc_idx=loadcase_idx;
            end
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
            
            %% careful.. only valid if front and rear spar have equal
            %% thickness
            %% TODOOO FIX
            f=1;
            
            tu=obj.t_sk_up;
            tl=obj.t_sk_lo;
            tw=obj.t_sp_fr+obj.t_sp_re;
            
            Iz=((obj.h*(obj.w*f)^3)-((obj.h-tu-tl)*(((obj.w*f)-(tw))^3)))/12;
            Ix=(((obj.w*f)*obj.h^3)-((((obj.w*f)-(tw))*(obj.h-tu-tl)^3)))/12;
            Ai=((obj.w)-tw/2)*(obj.h-tl);
            Iy=4*(1.0*Ai)^2/(2*(obj.h-tu)/(tw/2)+2*((obj.w)-tw/2)/(tu));
            % J=obj.rho_sp*tw*obj.h*(0.5*obj.w)^2+obj.rho_sk_up*tu*obj.w*(0.5*obj.h)^2;
            % J=Iy*obj.rho_sp;
            J=Ix+Iz;%(obj.rho_sp*tl*obj.w*(obj.h/2)^2*2+obj.rho_sp*tw*obj.h*(obj.w/2)^2+obj.rho_sp*(tw*obj.h)*(obj.h^2+obj.t_sp_fr^2)/12+2*obj.rho_sp*(tu*obj.w)*(obj.w^2+tu^2)/12);
            %J=obj.rho_sp*(tw*obj.h+2*tl*obj.w)*(obj.w^2+obj.h^2)/12;
            %y=4*Ai^2/(2*(obj.h-2*tu)/(tw/2)+2*(obj.w-tw)/tu);
            A=(obj.w*f)*obj.h-((obj.w*f)-tw)*(obj.h-tu-tl); % material crosssectional area
            
%             Ix=obj.w*obj.h^3/12;
%             Iz=obj.w^3*obj.h/12;
%             J=0.333*obj.w*obj.h^3;
%             J=0.85*(Ix+Iz);
%             Ix=1.288*obj.w*obj.h^3/12;
            %fueling_factor=0.82;
            
         %   fueling_factor=0.42;

            Aenclosed=obj.fueling_factor*((obj.w*f)-tw)*(obj.h-tu-tl);     % used to compute fuel volume
        end
    end   
end

