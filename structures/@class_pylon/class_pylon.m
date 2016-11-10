%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_pylon<class_beam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% Element positioning and type
        %> twist at each node position
        epsilon;   
        %> dihedral at each node position
        nu;
        %> sweep at each node position
        phi;                             
        %> flag if it is a rigid beam
        is_rigid=0;
        %% general values
        %> Fuselage volume
        
    end
    
    methods (Access=public)
        
         function obj = class_pylon(nel,crosssection,varargin)
            %call class_beam constructor 
            obj=obj@class_beam(nel,crosssection,varargin);
         end
         
         function struct_deflections=f_get_deflections(obj)
            struct_deflections.def=obj.nodal_deflections'; 
         end
         
         function obj=f_add_engine(obj,engine)
             obj.n_engines=obj.n_engines+1;
             if(obj.n_engines==1)
                obj.engine=engine;
             else
                obj.engine(obj.n_engines)=engine;
             end
         end
         
         function obj=f_add_gear(obj,gear)
             obj.n_gears=obj.n_gears+1;
             if(obj.n_gears==1)
                obj.gear=gear;
             else
                obj.gear(obj.n_gears)=gear;
             end
         end
         
         function obj=f_set_state(obj,state)
             obj.load_factor=state.load_factor;
             obj.loadcase_index=state.loadcase_index;
             obj.landingimpact=0;
         end
         
         % plot_critical_case_idx(wing,varargin);
            
         % initialize beam with estimated distributed mass
         % obj = f_init_wingmass(obj,weights);
         
         % calculate estimated wing mass
         obj = f_calc_mass(obj,weights);
         
         % initialize beam manually 
         obj = f_init_stdBeam(obj,E,G,m);
         
%          obj = f_init_coords(le,phi,nu,epsilon);
%          
%          obj = f_init_geomProp(obj,r,t);
         
         obj = f_init_std_structure(obj,nu,phi,Ix,Iy,Iz,J,A,le);
         
         % assemble system matrices for solving
         obj = f_assemble(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces);
         
         % initialize elements with material properties
         obj=f_init_material_properties(obj,structure);

         
%          function structmesh = get_struct_mesh(obj)
%             structmesh.y=0:obj.wing_frontview_length/obj.nel:obj.wing_frontview_length;
%          end     
         
%          function plot_externalforces(wing,add_eigenmass,add_fuelmass,engine,gear,varargin)
%             hold on
% 
%            
%             midp=zeros(3,1);
%             
%             norm_def=wing.nodal_deflections;%/max(wing.nodal_deflections);
%             
%             x_dist=zeros(length(wing.beamelement)+1,1);
%             y_dist=zeros(length(wing.beamelement)+1,1);
%             z_dist=zeros(length(wing.beamelement)+1,1);
%             
%                 mid_def=zeros(3,1);
%             mid_rot=zeros(3,1);
%             
%             for i=1:1:length(wing.beamelement)
%                 for j=1:3
%                     midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
%                     mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
%                     mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
%                 end
%                 
%                 midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
%               
%                 midp=midp_defl;
%                 h=wing.beamelement(i).crosssection.h;
%                 w=wing.beamelement(i).crosssection.w;
%                 
%                 le=wing.beamelement(i).le;
%                 nu=wing.beamelement(i).nu;
%                 twist=wing.beamelement(i).epsilon;
%                 z_dist(i)=h/2;
%                 x_dist(i)=w/2;
%                 z_dist(i+1)=h/2;
%                 x_dist(i+1)=w/2;
%                 plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu+mid_rot(1),twist+mid_rot(2),mid_rot(3)],[1 1 1 1 1 1 1 1],0.3,1);
%             end
%             
%             axis equal
%             grid on
%             
%             qx=cell2mat({wing.beamelement(:).qx});
%             qy=cell2mat({wing.beamelement(:).qy});
%             qz=cell2mat({wing.beamelement(:).qz});
%             halfspan=sum(cell2mat({wing.beamelement(:).le}));
%             
%             max_q=max([abs(qx),abs(qy),abs(qz)]);
%             
%             scale=2*max_q/halfspan;
%             
%             optargin = size(varargin,2);
%             if optargin==1
%                 scale=varargin{1};
%             end
%             
%             qmrot=zeros(3,length(wing.beamelement));
%             qfrot=zeros(3,length(wing.beamelement));
%             
%             for i=1:length(wing.beamelement)    
%                 if add_eigenmass
%                     T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
%                     qmrot(1:3,i)=T*[0,0,wing.beamelement(i).qm]';%wing.beamelement(i).qm];           % distributed loading due to eigenmass in x,y,z coordinates
%                 end
%                 
%                 if add_fuelmass
%                     T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
%                     qfrot(:,i)=T*[0;0;wing.beamelement(i).qf];           % distributed loading due to eigenmass in x,y,z coordinates
%                 end
%             end
% 
%             xxx=subplot(1,1,1);
%             %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_co
%             %ords(:,3),'-k','LineWidth',2)
%             wing.node_coords(:,1)=wing.node_coords(:,1)+wing.nodal_deflections(1:6:end);
%             wing.node_coords(:,2)=wing.node_coords(:,2)+wing.nodal_deflections(2:6:end);
%             wing.node_coords(:,3)=wing.node_coords(:,3)+wing.nodal_deflections(3:6:end);
%             
%             
%             wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
%             plotglobalforce(xxx,qz,scale,wing.node_coords(:,:),'z',0,'b',0); 
%             wing.node_coords(:,3)=wing.node_coords(:,3)-2.*z_dist;
%             plotglobalforce(xxx,qmrot(3,:),scale,wing.node_coords(:,:),'z',0.45,'g',0); 
%             %wing.node_coords(:,3)=wing.node_coords(:,3)+2*z_dist;
%             wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)+qmrot(3,:)'./scale;
%             plotglobalforce(xxx,qfrot(3,:),scale,wing.node_coords(:,:),'z',0.9,'m',0); 
%             wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)-qmrot(3,:)'./scale;
%             wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
%             
%             wing.node_coords(:,1)=wing.node_coords(:,1)+x_dist;
%             plotglobalforce(xxx,qx,scale,wing.node_coords(:,:),'x',0.3,'c',0); 
%             axis equal
%             grid on   
%             wing.node_coords(:,1)=wing.node_coords(:,1)-x_dist;
%             
%             wing.node_coords(:,1)=wing.node_coords(:,1)-wing.nodal_deflections(1:6:end);
%             wing.node_coords(:,2)=wing.node_coords(:,2)-wing.nodal_deflections(2:6:end);
%             wing.node_coords(:,3)=wing.node_coords(:,3)-wing.nodal_deflections(3:6:end);
%             %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_coords(:,3),'-x','LineWidth',2)
%             %plot3(engine.cgpos(1),engine.cgpos(2),engine.cgpos(3),'c+');
%             %plot3(gear.pos(1),gear.pos(2),gear.pos(3),'go');
%           end
%          
         % ================================================================
         %> @brief plot deformations of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
         function plot_deformations(fuselage,varargin)

         end         

          
         
        % ================================================================
         %> @brief plot geometry of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
          function plot_geometry(fuselage,varargin)

          end  
         
          function plot_t(fuselage,varargin)
              
          end
          
          
          function write_structure_tecplot(fuselage,fileID,beam_nr)

		end
        

 
    end
    
end

