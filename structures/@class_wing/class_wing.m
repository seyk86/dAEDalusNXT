%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_wing.m
%> @brief File containing the class for a finite element wing model
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
%> @brief class for a finite element wing beam model
%> This class is the base class for a finite element beam, the follwing
%> coordinate system is used
%> Coordinate Definition (according to Airbus Wing Geometry Definitions ):
%>   y= along starboard wing
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
% ======================================================================

%wing is a subclass of beam
classdef class_wing < class_beam
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties

        % this class has all properties of class_beam + following
        
        %> vector of nodal_deflections in all DOF's in quarter chord line
        nodal_deflections_c4;           
        %> distance between quarter chord line and shear center
        dist_c4_sc; 
        %> twist at each node position
        epsilon;   
        
        %>
        nu;
        
        %>
        phi;
        
        deltanu;                        
        %> number of engines
        n_engines=0;
        %> number of wing mounted gears
        n_gears=0; 
        %> array of engine classes    
        engine;
        %> array of engine classes
        gear; 
        
        %> number of ribs
        nrib;
        %> rib spacing on wing runlength
        rib_spacing;                   
        
        landingimpact=0;
        
        %% general values
        %> wing area
        Awing=1.0;
        %> wingbox volume
        V_wingbox=0;
        %> wing surface estimate
        S_wing=0;
        
        %> wingbox total mass
        m_wingbox_total=0.0;
        %> rib mass
        m_ribs=0.0;
        %> fuel mass total
        m_fuel_total=0.0;
        %> fuel volume
        fuel_volume_total=0.0;
        %> fueling state
        fueling_state;

        
        %> flag for wingmounted gears
        wingmountedgears=0;
        %> flag for wingmounted engines
        wingmountedengines=0;
        
        % euler angles
        % gamma=0.0;  % a/c rotation in pitch (for g-load direction)    [rad]
        % add later: direction of g (now in fixed in -z direction)
        
        %> local y coordinate runlength of wing from root to tip
        wing_frontview_length=0;       

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=public)
         % constructor
         function obj = class_wing(nel,crosssection,varargin)
            %call class_beam constructor 
            obj=obj@class_beam(nel,crosssection,varargin);
            obj.deltanu=zeros((nel*6+6)/3,1);
         end
         
         % returns deflections at quarter chord line of wing ( for
         % aerodynamic mesh morphing )
         function struct_deflections=f_get_deflections_c4(obj)
            struct_deflections.def=obj.nodal_deflections_c4'; 
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
         
         plot_critical_case_idx(wing,varargin);
         
         
         wingstr= f_set_aeroloads(wingstr,wingaero);
            
         % initialize beam with estimated distributed mass
         obj = f_init_wingmass(obj,weights);
         
         obj= f_calc_mass(obj,weights);
         
         % initialize beam with standard b
         obj = f_init_stdBeam(obj,A,E,G,Iy,Iz,Ip);
         
         % assemble system matrices for solving
         %obj=f_assemble(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces);

         % initialize elements with material properties
         obj=f_init_material_properties(obj,structure);
         
         function obj=f_set_fueling_state(obj,fueling_state)
             
         end
         
         %wing structural layout

         
         
         function structmesh= get_struct_mesh(obj)
            structmesh.y=0:obj.wing_frontview_length/obj.nel:obj.wing_frontview_length;
         end
         % calculate estimated wing mass
        % obj = f_calc_mass(obj,structure,weights );
         
         
         function plot_externalforces(wing,add_eigenmass,add_fuelmass,engine,gear,varargin)
            hold on

           
            midp=zeros(3,1);
            
            norm_def=wing.nodal_deflections;%/max(wing.nodal_deflections);
            
            x_dist=zeros(length(wing.beamelement)+1,1);
            y_dist=zeros(length(wing.beamelement)+1,1);
            z_dist=zeros(length(wing.beamelement)+1,1);
            
                mid_def=zeros(3,1);
            mid_rot=zeros(3,1);
            
            for i=1:1:length(wing.beamelement)
                for j=1:3
                    midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                    mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
                    mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
                end
                
                midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
              
                midp=midp_defl;
                h=wing.beamelement(i).crosssection.h;
                w=wing.beamelement(i).crosssection.w;
                
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                twist=wing.beamelement(i).epsilon;
                z_dist(i)=h/2;
                x_dist(i)=w/2;
                z_dist(i+1)=h/2;
                x_dist(i+1)=w/2;
                plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu+mid_rot(1),twist+mid_rot(2),mid_rot(3)],[1 1 1 1 1 1 1 1],0.3,1);
            end
            
            axis equal
            grid on
            
            qx=cell2mat({wing.beamelement(:).qx});
            qy=cell2mat({wing.beamelement(:).qy});
            qz=cell2mat({wing.beamelement(:).qz});
            halfspan=sum(cell2mat({wing.beamelement(:).le}));
            
            max_q=max([abs(qx),abs(qy),abs(qz)]);
            
            scale=2*max_q/halfspan;
            
            optargin = size(varargin,2);
            if optargin==1
                scale=varargin{1};
            end
            
            qmrot=zeros(3,length(wing.beamelement));
            qfrot=zeros(3,length(wing.beamelement));
            
            for i=1:length(wing.beamelement)    
                if add_eigenmass
                    T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
                    qmrot(1:3,i)=T*[0,0,wing.beamelement(i).qm]';%wing.beamelement(i).qm];           % distributed loading due to eigenmass in x,y,z coordinates
                end
                
                if add_fuelmass
                    T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
                    qfrot(:,i)=T*[0;0;wing.beamelement(i).qf];           % distributed loading due to eigenmass in x,y,z coordinates
                end
            end

            xxx=subplot(1,1,1);
            %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_co
            %ords(:,3),'-k','LineWidth',2)
            wing.node_coords(:,1)=wing.node_coords(:,1)+wing.nodal_deflections(1:6:end);
            wing.node_coords(:,2)=wing.node_coords(:,2)+wing.nodal_deflections(2:6:end);
            wing.node_coords(:,3)=wing.node_coords(:,3)+wing.nodal_deflections(3:6:end);
            
            
            wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
            plotglobalforce(xxx,qz,scale,wing.node_coords(:,:),'z',0,'b',0); 
            wing.node_coords(:,3)=wing.node_coords(:,3)-2.*z_dist;
            plotglobalforce(xxx,qmrot(3,:),scale,wing.node_coords(:,:),'z',0.45,'g',0); 
            %wing.node_coords(:,3)=wing.node_coords(:,3)+2*z_dist;
            wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)+qmrot(3,:)'./scale;
            plotglobalforce(xxx,qfrot(3,:),scale,wing.node_coords(:,:),'z',0.9,'m',0); 
            wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)-qmrot(3,:)'./scale;
            wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
            
            wing.node_coords(:,1)=wing.node_coords(:,1)+x_dist;
            plotglobalforce(xxx,qx,scale,wing.node_coords(:,:),'x',0.3,'c',0); 
            axis equal
            grid on   
            wing.node_coords(:,1)=wing.node_coords(:,1)-x_dist;
            
            wing.node_coords(:,1)=wing.node_coords(:,1)-wing.nodal_deflections(1:6:end);
            wing.node_coords(:,2)=wing.node_coords(:,2)-wing.nodal_deflections(2:6:end);
            wing.node_coords(:,3)=wing.node_coords(:,3)-wing.nodal_deflections(3:6:end);
            %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_coords(:,3),'-x','LineWidth',2)
            %plot3(engine.cgpos(1),engine.cgpos(2),engine.cgpos(3),'c+');
            %plot3(gear.pos(1),gear.pos(2),gear.pos(3),'go');
          end
         
         % ================================================================
         %> @brief plot deformations of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
         function plot_deformations(wing,varargin)
            hold on
            midp=zeros(3,1);
            mid_def=zeros(3,1);
            mid_rot=zeros(3,1);
            norm_def=wing.nodal_deflections;%/max(wing.nodal_deflections);
            
            col=0.1;
            tr=0.5;

            optargin = size(varargin,2);
            if optargin==1
               xcol=varargin{1};
               xcol=cell2mat(xcol);
               col=xcol(1);
               tr=xcol(2);
            end

            
            for i=1:1:length(wing.beamelement)
                for j=1:3
                    midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                    mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
                    mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
                end

                mid_str=0;%abs((norm_def(2+6*(i-1))+norm_def(6+2+6*(i-1))))/2;
                
                mid_rot=wing.beamelement(i).T(1:3,1:3)'*mid_rot;
                midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
                h=wing.beamelement(i).crosssection.h;
                w=wing.beamelement(i).crosssection.w;
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                twist=wing.beamelement(i).epsilon;
                sweep=wing.beamelement(i).phi;
                
                %plotcube([midp(1),midp(2),midp(3)],[w,le,h],[twist,0,nu],[1 1 1 1 1 1 1 1],transparency,1);
                plotcube([midp_defl(1),midp_defl(2),midp_defl(3)],[w,le+mid_str,h],[mid_rot(1)+nu,mid_rot(2)+twist,mid_rot(3)+sweep],[col col col col col col col col],tr,1);
            end
            axis equal
            grid on
         end         

          
         
        % ================================================================
         %> @brief plot geometry of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
         function plot_geometry(wing,varargin)
            hold on

            col=1;
            tr=0.5;
            
            optargin = size(varargin,2);
            if optargin==1
               xcol=varargin{1};
               xcol=cell2mat(xcol);
               col=xcol(1);
               tr=xcol(2);
            end

            midp=zeros(3,1);
            for i=1:1:length(wing.beamelement)
                for j=1:3
                    midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                end
                h=wing.beamelement(i).crosssection.h;
                w=wing.beamelement(i).crosssection.w;%/cos(wing.beamelement(i).phi);
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                twist=wing.beamelement(i).epsilon;
                sweep=wing.beamelement(i).phi;
                %le=wing.beamelement(i).le*cos(sweep);
                plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu,twist,sweep],[col col col col col col col col],tr,1);
            end
            axis equal
            grid on
         end  
         
         function plot_structure(wing,varargin)
            mycolormap=jet(256);
            hold on
            midp=zeros(3,1);
            mid_def=zeros(3,1);
            mid_rot=zeros(3,1);
            norm_def=wing.nodal_deflections;
                for i=1:1:length(wing.beamelement) 
                    t_sk_up(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_up});
                    t_sk_lo(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_lo});
                    t_sp_fr(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_fr});
                    t_sp_re(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_re});
                end
            
            optargin = size(varargin,2);
            if optargin==2
                t_min=varargin{2};
                t_max=varargin{1};
            else
                 t_max=0;
                t_min=100;
                tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re]
                t_max=max([tk t_max])
                t_min=min([tk t_min])
            end
            
            for i=1:1:length(wing.beamelement)
                
                for j=1:3
                    midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                                        mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
                    mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
                end
                
                mid_rot=wing.beamelement(i).T(1:3,1:3)'*mid_rot;
                midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
                
                
                h=wing.beamelement(i).crosssection.h;
                w=wing.beamelement(i).crosssection.w;
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                sweep=wing.beamelement(i).phi*0;
                le=wing.beamelement(i).le*cos(wing.beamelement(i).phi);
                twist=wing.beamelement(i).epsilon;

                col_i=min(round((t_sk_lo(i)-t_min)*255/(t_max-t_min)+1),255);
                colidx(i,1)=col_i;
                col=ones(8,1)*mycolormap(col_i,:);
                
                a=-wing.beamelement(i).nu;
                b=-wing.beamelement(i).epsilon;
                %lower skin
                plotcube([midp(1),midp(2)-h/2*sin(a),midp(3)-h/2*cos(a)],[w,le,t_sk_lo(i)],[nu,twist,sweep],ones(8,1)*col_i/255,1,1);
               
                
                % front spar
                col_i=min(round((t_sp_fr(i)-t_min)*255/(t_max-t_min)+1),255);
                colidx(i,2)=col_i;
                col=ones(8,1)*mycolormap(col_i,:);
                
                % rot=wing.nodal_deflections((i-1)*wing.el_ndof+4:(i-1)*wing.el_ndof+6);
                

                c=0;
              
                Lx=[1       0       0
                0   cos(a)  sin(a)
                0   -sin(a) cos(a)];

                Ly=[cos(b) 0 -sin(b)
                0      1    0
                sin(b)  0   cos(b)];


                M=Ly*Lx;
                
                diff=[-w/2 0 0];
                
                xxx=M*diff';

                plotcube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_fr(i),le,h],[nu,twist,sweep],ones(8,1)*col_i/255,1,1);
                %rear spar
                
                diff=[w/2 0 0];
                xxx=M*diff';
                col_i=min(round((t_sp_re(i)-t_min)*255/(t_max-t_min)+1),255);
                plotcube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_re(i),le,h],[nu,twist,sweep],ones(8,1)*col_i/255,1,1);
 
                %upper skin
                col_i=min(round((t_sk_up(i)-t_min)*255/(t_max-t_min)+1),255);
                col=ones(8,1)*mycolormap(col_i,:);
                plotcube([midp(1),midp(2)+h/2*sin(a),midp(3)+h/2*cos(a)],[w,le,t_sk_up(i)],[nu,twist,sweep],ones(8,1)*col_i/255,1,1);

            end
            axis equal
            grid on
            
            if ~isempty(wing.engine)
                for ii=1:length(wing.engine)
                     %plot3(wing.engine(ii).cg_pos(1),wing.engine(ii).cg_pos(2),wing.engine(ii).cg_pos(3),'kx','MarkerSize',20);

                            %           find closest node:
                            dist=zeros(wing.nel,1);
                            for i=1:wing.nel
                                dist(i)=sqrt(sum((wing.node_coords(i,2)-wing.engine(ii).cg_pos(2)').^2));
                            end
                            [Y,I] = min(dist);
                            hold on
                     
                     plot3([wing.engine(ii).cg_pos(1) wing.node_coords(I,1)],[wing.engine(ii).cg_pos(2) wing.node_coords(I,2)],[wing.engine(ii).cg_pos(3) wing.node_coords(I,3)],'-k','LineWidth',2);
                     plot3(wing.engine(ii).cg_pos(1),wing.engine(ii).cg_pos(2),wing.engine(ii).cg_pos(3),'-kx','MarkerSize',15,'LineWidth',2);
                end
            end
            
           if ~isempty(wing.gear)
                for ii=1:length(wing.gear)
                     %plot3(wing.engine(ii).cg_pos(1),wing.engine(ii).cg_pos(2),wing.engine(ii).cg_pos(3),'kx','MarkerSize',20);

                            %           find closest node:
                            dist=zeros(wing.nel,1);
                            for i=1:wing.nel
                                dist(i)=sqrt(sum((wing.node_coords(i,2)-wing.gear(ii).pos(2)').^2));
                            end
                            [Y,I] = min(dist);
                            hold on
                     
                     plot3([wing.gear(ii).pos(1) wing.node_coords(I,1)],[wing.gear(ii).pos(2) wing.node_coords(I,2)],[wing.gear(ii).pos(3) wing.node_coords(I,3)],'-k','LineWidth',2);
                     plot3(wing.gear(ii).pos(1),wing.gear(ii).pos(2),wing.gear(ii).pos(3),'-ko','MarkerSize',15,'LineWidth',2);
                end
           end
         end

         function obj=write_structure_tecplot(wing,fileID,beam_nr)
            
             midp=zeros(3,1);
             mid_def=zeros(3,1);
             mid_rot=zeros(3,1);
             norm_def=wing.nodal_deflections;
             for i=1:1:length(wing.beamelement)
                 t_sk_up(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_up});
                 t_sk_lo(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_lo});
                 t_sp_fr(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_fr});
                 t_sp_re(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_re});
             end

             for i=1:1:length(wing.beamelement)
                 
                 for j=1:3
                     midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                     mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
                     mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
                 end
                 
                 mid_rot=wing.beamelement(i).T(1:3,1:3)'*mid_rot;
                 %midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
                 
                 h=wing.beamelement(i).crosssection.h;
                 w=wing.beamelement(i).crosssection.w;
                 le=wing.beamelement(i).le;
                 nu=wing.beamelement(i).nu;
                 sweep=wing.beamelement(i).phi;
                 le=wing.beamelement(i).le;%*cos(wing.beamelement(i).phi);
                 twist=wing.beamelement(i).epsilon;
                
                 
                 a=-wing.beamelement(i).nu;
                 b=-wing.beamelement(i).epsilon;
                 % Lower Skin
                 fprintf(fileID,'ZONE T="%s_%s_1" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                 [c1, c2, c3, c4, c5, c6, c7, c8] = get_cube([midp(1),midp(2)-h/2*sin(a),midp(3)-h/2*cos(a)],[w,le,t_sk_lo(i)],[nu,twist,sweep]);
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',t_sk_lo(i),t_sk_lo(i),t_sk_lo(i),t_sk_lo(i),t_sk_lo(i),t_sk_lo(i),t_sk_lo(i),t_sk_lo(i));
                 
                 % front spar

                 Lx=[1       0       0
                     0   cos(a)  sin(a)
                     0   -sin(a) cos(a)];
                 
                 Ly=[cos(b) 0 -sin(b)
                     0      1    0
                     sin(b)  0   cos(b)];
                 
                 
                 M=Ly*Lx;
                 
                 diff=[-w/2 0 0];
                 
                 xxx=M*diff';
                 
                 fprintf(fileID,'ZONE T="%s_%s_2" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                 [c1, c2, c3, c4, c5, c6, c7, c8] = get_cube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_fr(i),le,h],[nu,twist,sweep]);
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',t_sp_fr(i),t_sp_fr(i),t_sp_fr(i),t_sp_fr(i),t_sp_fr(i),t_sp_fr(i),t_sp_fr(i),t_sp_fr(i));
                 
                 %rear spar
                 
                 diff=[w/2 0 0];
                 xxx=M*diff';
                 
                 fprintf(fileID,'ZONE T="%s_%s_3" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                 [c1, c2, c3, c4, c5, c6, c7, c8] = get_cube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_re(i),le,h],[nu,twist,sweep]);
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',t_sp_re(i),t_sp_re(i),t_sp_re(i),t_sp_re(i),t_sp_re(i),t_sp_re(i),t_sp_re(i),t_sp_re(i));

                 %upper skin
                 fprintf(fileID,'ZONE T="%s_%s_4" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                 [c1, c2, c3, c4, c5, c6, c7, c8] = get_cube([midp(1),midp(2)+h/2*sin(a),midp(3)+h/2*cos(a)],[w,le,t_sk_up(i)],[nu,twist,sweep]);
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',t_sk_up(i),t_sk_up(i),t_sk_up(i),t_sk_up(i),t_sk_up(i),t_sk_up(i),t_sk_up(i),t_sk_up(i));
             end
         end
             
         function obj=write_fuel_tecplot(wing,fileID,beam_nr)
            
             midp=zeros(3,1);
             mid_def=zeros(3,1);
             mid_rot=zeros(3,1);
             norm_def=wing.nodal_deflections;
             for i=1:1:length(wing.beamelement)
                 t_sk_up(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_up});
                 t_sk_lo(i)=cell2mat({wing.beamelement(i).crosssection.t_sk_lo});
                 t_sp_fr(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_fr});
                 t_sp_re(i)=cell2mat({wing.beamelement(i).crosssection.t_sp_re});
             end

             for i=1:1:length(wing.beamelement)
                 
                 for j=1:3
                     midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                     mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
                     mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
                 end
                 
                 mid_rot=wing.beamelement(i).T(1:3,1:3)'*mid_rot;
                 %midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
                 
                 h=wing.beamelement(i).crosssection.h;
                 w=wing.beamelement(i).crosssection.w;
                 le=wing.beamelement(i).le;
                 nu=wing.beamelement(i).nu;
                 sweep=wing.beamelement(i).phi;
                 is_fueled=wing.beamelement(i).is_fueled;
                 le=wing.beamelement(i).le;%*cos(wing.beamelement(i).phi);
                 twist=wing.beamelement(i).epsilon;
                 
                 a=-wing.beamelement(i).nu;
                 b=-wing.beamelement(i).epsilon;
                 % Lower Skin
                 fprintf(fileID,'ZONE T="%s_%s_1" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                 [c1, c2, c3, c4, c5, c6, c7, c8] = get_cube([midp(1),midp(2),midp(3)],[w-t_sp_fr(i)-t_sp_re(i),le,h-t_sk_lo(i)-t_sk_up(i)],[nu,twist,sweep]);
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3));
                 fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',is_fueled,is_fueled,is_fueled,is_fueled,is_fueled,is_fueled,is_fueled,is_fueled);
             end  
         end

         
         function plot_internalforces(wing,varargin)
            delta_x=0;
                delta_y=0;
                delta_z=0;
                med_x=0;
                med_y=0;
                med_z=0;
                scale=0;
                delta=0;
                deltaxy=0;
                
            aQx=wing.node_loadings_loc(1:6:end);
            aQy=wing.node_loadings_loc(2:6:end);
            aQz=wing.node_loadings_loc(3:6:end);
            aMx=wing.node_loadings_loc(4:6:end);
            aMy=wing.node_loadings_loc(5:6:end);
            aMz=wing.node_loadings_loc(6:6:end);

            optargin = size(varargin,2);
            if optargin==15
                delta_x=varargin{1};
                delta_y=varargin{2};
                delta_z=varargin{3};
                med_x=varargin{4};
                med_y=varargin{5};
                med_z=varargin{6};
                scale=0.6*sqrt((delta_x)^2+(delta_z)^2+(delta_y)^2);
                delta=max([delta_y,delta_z]);
                deltaxy=max([delta_y,delta_x]);
                scale_bendingmoment=varargin{7};
                scale_torsionmoment=varargin{8};
                scale_force=varargin{9};
                maxMx=varargin{10};
                maxMt=varargin{11};
                maxMz=varargin{12};
                maxQx=varargin{13};
                maxQy=varargin{14};
                maxQz=varargin{15};
            else 
                max_x=max(wing.node_coords(:,1));
                min_x=min(wing.node_coords(:,1));
            
                max_y=max(wing.node_coords(:,2));
                min_y=min(wing.node_coords(:,2));
            
                max_z=max(wing.node_coords(:,3));
                min_z=min(wing.node_coords(:,3)); 
             
                delta_x=max_x-min_x;
                delta_z=max_z-min_z;
                delta_y=max_y-min_y;

                med_y=(max_y+min_y)/2;
                med_z=(max_z+min_z)/2;
                med_x=(max_x+min_x)/2;
                
                 scale=0.6*sqrt((delta_x)^2+(delta_z)^2+(delta_y)^2);
                delta=max([delta_y,delta_z]);
                deltaxy=max([delta_y,delta_x]);
                
                     max_Mb=max([abs(aMx)' abs(aMz)' abs(aMy)']);
                    %scale torsional moments
                 max_Mt=max_Mb;
                 max_Q =max([abs(aQx)' abs(aQy)' abs(aQz)']);
            
                    scale_bendingmoment=(2*max_Mb)/scale;
                    scale_torsionmoment=(2*max_Mt)/scale;
                                scale_force=(2*max_Q)/scale;
                                
                maxMx=max(abs(aMx));
                maxMt=max(abs(aMy));
                maxMz=max(abs(aMz));
                maxQx=max(abs(aQx));
                maxQy=max(abs(aQy));
                maxQz=max(abs(aQz));
            end
             
            % plot Mx 
            handle1=subplot(3,2,2);
            plotforce(handle1,aMx,scale_bendingmoment,wing.node_coords,'z',1,'r',1);
           
            title('Mx')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            ylabel('Y')
            zlabel('Z')
            text(0,delta,delta/1.5,['Mx_{max}=' sprintf('%4.3gMNm',maxMx/10^6)],'BackgroundColor','white');
            axis equal
            
            % plot Qz 
            handle2=subplot(3,2,1);
            plotforce(handle2,aQz,scale_force,wing.node_coords,'z',1,'r',1);
            title('Qz')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            ylabel('Y')
            zlabel('Z')
            text(0,delta,delta/1.5,['Qz_{max}=' sprintf('%4.3gMN',maxQz/10^6)],'BackgroundColor','white');
            axis equal
            
            handle4=subplot(3,2,3);
            plotforce(handle4,aQy,scale_force,wing.node_coords,'z',1,'r',1);
            title('Qy')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            ylabel('Y')
            zlabel('Z')
            text(0,delta,delta/1.5,['Qy_{max}=' sprintf('%4.3gMN',maxQy/10^6)],'BackgroundColor','white');
            axis equal
                
            handle5=subplot(3,2,4);
            plotforce(handle5,aMy,scale_torsionmoment,wing.node_coords,'z',1,'r',1);
            title('Mt')
            view(-90,0)
            ylabel('Y')
            zlabel('Z')
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            text(0,delta,delta/1.5,['Mt_{max}=' sprintf('%4.3gMNm',maxMt/10^6)],'BackgroundColor','white');
            axis equal
            
            handle3=subplot(3,2,6);
            plotforce(handle3,aMz,scale_bendingmoment,wing.node_coords,'x',1,'r',1);
            title('Mz')
            view(-90,90)
            ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
            xlim([med_x-0.3*deltaxy,med_x+0.3*deltaxy]);
            ylabel('Y')
            xlabel('X')
              text(0,delta,-delta/2,['Mz_{max}=' sprintf('%4.3gMNm',maxMz/10^6)],'BackgroundColor','white');
            axis equal
            
            handle6=subplot(3,2,5);
            plotforce(handle6,aQx,scale_force,wing.node_coords,'x',1,'r',1);
            title('Qx')
            view(-90,90)
            ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
            xlim([med_x-0.3*deltaxy,med_x+0.3*deltaxy]);
            ylabel('Y')
            xlabel('X')
            text(0,delta,-delta/2,['Qx_{max}=' sprintf('%4.3gMN',maxQx/10^6)],'BackgroundColor','white');
            axis equal
         end   
    end  
end

