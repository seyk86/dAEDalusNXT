%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_fuselage_geometry
 
    properties
        
        name;
        
        isExternalFEM = 0;  % 0 if fuselage should be selfdesigned in Daedalus, 
                            % 1 if fuselage is from external FEM model
        pathNodeCoords;     % if fuselage is from external FEM model, gives path
                            % of node coordinates matlab file. Nodes from
                            % nose to tail
        pathMassMatrix;         % if fuselage is from external FEM model, gives path
                                % of mass matrix matlab file.
        pathStiffnessMatrix;    % if fuselage is from external FEM model, gives path
                                % of stiffness matrix matlab file.
        
        fuselage_segments;
        % total length
        l;
        % effective diameter
        d_eff;
        % surface wetted area
        S_wet;
        
        % center points for the defined grid
        center_coords;
        % shell width in each node defined in center_coords
        shell_width;
        % shell height in each node defined in center_coords
        shell_height;
        % coordinates of the points that form an elliptical segment
        shell_coords

        % volumetric grid parameters       
        grid_vol;
        panels_vol;
        is_te_vol;
        opposite_te_vol;
        te_idx_vol;
        %cruzifix grid parameters
        grid_flat;
        grid_deflected
        panels_flat;
        is_te_flat;
        opposite_te_flat;
        te_idx_flat;
        
        % parasite and upsweep drag coefficient
        CD_f;
        % parasite drag force
        D_f;
        
        % for 2D grid
        grid_start_idx_flat;
        panel_start_idx_flat;
        
        % for 3D grid
        grid_start_idx_vol;
        panel_start_idx_vol;
        
        
        
    end
    
    methods
        
        function obj=class_fuselage_geometry(varargin)
            if nargin==1
                obj=obj.read_xml_definition(varargin{1});
                obj=obj.compute_shell_coords();
                obj=obj.complete_params();
            elseif nargin==6
                pos_n=varargin{1};
                obj.l=varargin{2};
                w=varargin{3};
                h=varargin{4};
                l_n=varargin{5};
                l_t=varargin{6};
                l_c=obj.l-l_n-l_t;
                pos_c=pos_n+[l_n 0 -l_n*sin(-7*pi/180)];
                obj.fuselage_segments=class_fuselagesegment.empty(0,3);
                obj.fuselage_segments(2)=class_fuselagesegment('Cabin',pos_c,l_c,w,h);
                obj.fuselage_segments(1)=obj.fuselage_segments(2).set_nose(obj.fuselage_segments(2),l_n);
                obj.fuselage_segments(3)=obj.fuselage_segments(2).set_tail(obj.fuselage_segments(2),l_t);
                obj=obj.compute_shell_coords();
                obj=obj.complete_params();
            end    
        end
        
        function obj=compute_shell_coords(obj,varargin)
            
            center_coords=[];
            shell_width=[];
            shell_height=[];
            shell_coords=[];
            
            iNodes = 1;
            for i=1:length(obj.fuselage_segments)
                if nargin==2
                    le_max=varargin{1};
                    n=ceil(obj.fuselage_segments(i).l/le_max);
                    if obj.isExternalFEM==0
                        obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_shell_coords(n);
                    elseif obj.isExternalFEM==1
                        obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_shell_coords(n,iNodes,obj.pathNodeCoords);
                    end
                else
                    obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_shell_coords();
                end   
                if isempty(center_coords)
                    center_coords=[center_coords obj.fuselage_segments(i).center_coords];
                    shell_width=[shell_width obj.fuselage_segments(i).shell_width];
                    shell_height=[shell_height obj.fuselage_segments(i).shell_height];
                else
                    center_coords=[center_coords obj.fuselage_segments(i).center_coords(:,2:end,:)];
                    shell_width=[shell_width obj.fuselage_segments(i).shell_width(:,2:end,:)];
                    shell_height=[shell_height obj.fuselage_segments(i).shell_height(:,2:end,:)];
                end
                shell_coords=[shell_coords obj.fuselage_segments(i).shell_coords];
                iNodes = size(center_coords,2);
            end
            obj.center_coords=center_coords;
            obj.shell_width=shell_width;
            obj.shell_height=shell_height;
            obj.shell_coords=shell_coords;
            
        end
        
        
        function obj=update_idx(obj)
            start=1;
            for i=1:length(obj.fuselage_segments)
                
                obj.fuselage_segments(i).grid_start_idx_flat=start+(obj.grid_start_idx_flat-1);
                start=start+length(obj.fuselage_segments(i).grid_flat);
                %obj.fuselage_segments(i).panel_start_idx=obj.wing_segments(i).panel_start_idx+(obj.panel_start_idx-1);
            end
        end
        
        
        function obj=complete_params(obj)
            
            l=0;
            w_max=[];
            h_max=[];
            
            for i=1:length(obj.fuselage_segments)
                l=l+obj.fuselage_segments(i).l;
                w_max=[w_max max(obj.fuselage_segments(i).w_f,obj.fuselage_segments(i).w_r)];
                h_max=[h_max max(obj.fuselage_segments(i).h_f,obj.fuselage_segments(i).h_r)];
            end
            w_max=max(w_max);
            h_max=max(h_max);
            
            obj.l=l;
            obj.d_eff=2*(w_max*h_max)^0.5;           
                
        end 
        
        
        function obj=compute_grid_flat(obj,z_coord,x_max)
            grid=[];
            panels=[];
            te_idx=[];
            len=[];
            
            n_pan_horz=ceil((abs(max(obj.grid_vol(2,:)))+abs(min(obj.grid_vol(2,:))))/x_max);
            n_pan_vert=ceil((abs(max(obj.grid_vol(3,:)))+abs(min(obj.grid_vol(3,:))))/x_max);
            
            for i=1:length(obj.fuselage_segments)
                obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_grid_flat(z_coord,n_pan_horz,n_pan_vert);
                len=length(grid);
                grid=[grid obj.fuselage_segments(i).grid_flat];
                te_idx=[te_idx,obj.fuselage_segments(i).te_idx_flat+len];
                panels=[panels,obj.fuselage_segments(i).panels_flat+len];
            end
            
            obj.grid_flat=grid;
            obj.panels_flat=panels;
            obj.te_idx_flat=te_idx;
            
%             inv_pan_fus(1,:)=obj.panels(1,:);
%             inv_pan_fus(2,:)=obj.panels(4,:);
%             inv_pan_fus(3,:)=obj.panels(3,:);
%             inv_pan_fus(4,:)=obj.panels(2,:);
%             
%             obj.grid_flat=obj.grid;
%             obj.panels_flat=inv_pan_fus;
%             obj.is_te_flat=zeros(1,length(obj.panels));
%             obj.opposite_te_flat=zeros(1,length(obj.grid));
%             obj.te_idx_flat=zeros(1,length(obj.panels));
        end
        
        
            
        function obj=compute_grid(obj)
            grid=[];
            panels=[];
            te_idx=[];
            len=[];
            
            for i=1:length(obj.fuselage_segments)
                obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_grid();
                len=length(grid);
                grid=[grid obj.fuselage_segments(i).grid];
                te_idx=[te_idx,obj.fuselage_segments(i).te_idx+len];
                panels=[panels,obj.fuselage_segments(i).panels+len];
            end
            
            obj.grid_vol=grid;
            obj.panels_vol=panels;
            obj.te_idx_vol=te_idx;
            
             inv_pan_fus(1,:)=obj.panels_vol(1,:);
             inv_pan_fus(2,:)=obj.panels_vol(4,:);
             inv_pan_fus(3,:)=obj.panels_vol(3,:);
             inv_pan_fus(4,:)=obj.panels_vol(2,:);
            
           
            obj.panels_vol=inv_pan_fus;
            obj.is_te_vol=zeros(1,length(obj.panels_vol));
            obj.opposite_te_vol=zeros(1,length(obj.grid_vol));
            obj.te_idx_vol=zeros(1,length(obj.panels_vol));
        end
        
        function grid_deflected=compute_deflected_grid_flat(obj,panels,grid_deflected,deflections,z_coord,x_max)
            start_idx=1;
            n_pan_horz=ceil((abs(max(obj.grid_vol(2,:)))+abs(min(obj.grid_vol(2,:))))/x_max);
            n_pan_vert=ceil((abs(max(obj.grid_vol(3,:)))+abs(min(obj.grid_vol(3,:))))/x_max);
            if length(deflections)==6*size(obj.center_coords,2) %% if is a fuselage
                for i=1:length(obj.fuselage_segments)
                    end_idx=start_idx+size(obj.fuselage_segments(i).center_coords,2)*6-1;
                    grid_deflected=obj.fuselage_segments(i).compute_grid_deflected(panels,grid_deflected,deflections(:,start_idx:end_idx),z_coord,n_pan_horz,n_pan_vert);
                    start_idx=end_idx-5;
                end
            else % else it is a nacelle
                 for i=1:length(obj.fuselage_segments)
                    grid_deflected=obj.fuselage_segments(i).compute_grid_deflected(panels,grid_deflected,deflections(:,end-5:end),z_coord,n_pan_horz,n_pan_vert);
                 end
            end
        end
        
        
        function plot_fuselage_surface(fuselage_surface)
            
            for i=1:length(fuselage_surface.fuselage_segments)
                fuselage_surface.fuselage_segments(i).plot_elliptical_segment();
            end
            axis equal
            axis tight
            grid on

        end
        
        function obj=plot_grid(obj)
            hold on
            for i=1:length(obj.panels_flat)
                handle= fill3(obj.grid_flat(1,obj.panels_flat(:,i)), obj.grid_flat(2,obj.panels_flat(:,i)),obj.grid_flat(3,obj.panels_flat(:,i)),'b');
                alpha(handle,0.4);
            end
            axis equal
            axis tight
            grid on
        end
        
        function obj=plot_grid_vol(obj)
            hold on
            for i=1:length(obj.panels_vol)
                handle= fill3(obj.grid_vol(1,obj.panels_vol(:,i)), obj.grid_vol(2,obj.panels_vol(:,i)),obj.grid_vol(3,obj.panels_vol(:,i)),'b');
                alpha(handle,1);
            end
            axis equal
            axis tight
            grid on
        end
        
        function obj=complete_fuselage(obj,varargin)
            obj.fuselage_segments(1,3)=class_fuselagesegment;
            obj.fuselage_segments(2)=obj.fuselage_segments(1);
            obj.fuselage_segments(1)=obj.fuselage_segments(2).set_nose(obj.fuselage_segments(2));
            obj.fuselage_segments(3)=obj.fuselage_segments(2).set_tail(obj.fuselage_segments(2));
            obj=obj.complete_params();
        end
        
        function obj=compute_wetted_area(obj)
            obj.S_wet=0;
            for i=1:length(obj.fuselage_segments)
                obj.fuselage_segments(i)=obj.fuselage_segments(i).compute_wetted_area();
                obj.S_wet=obj.S_wet+obj.fuselage_segments(i).S_wet;
            end
        end
        
        function obj=compute_friction_drag(obj,state,S_ref)
            
            % flat plate CF for laminar or turbulent flow
            Re_l=state.rho_air*norm(state.V_inf)*obj.l/state.mu;
            if Re_l>5e5
                CF=0.455/log10(Re_l)^2.58;
            else
                CF=1.328/sqrt(Re_l);
            end
            
            % form factor for bodies: http://adg.stanford.edu/aa241/drag/BODYFORMFACTOR.HTML
            d=obj.d_eff/obj.l;
            D=sqrt(1-(1-state.Ma^2)*d^2);
            a=2*(1-state.Ma^2)*d^2/D^3*(atanh(D)-D);
            du_maxU0=a/(2-a)/sqrt(1-state.Ma^2);
            C=2.3;
            k=(1+C*du_maxU0)^2;
            
            % aft fuselage upsweep drag: http://adg.stanford.edu/aa241/drag/upsweepdrag.html
            sweep_t=0;
            for i=1:length(obj.fuselage_segments)
                if strcmp(obj.fuselage_segments(i).name,'Tail')
                    sweep_t=abs(obj.fuselage_segments(i).sweep);
                end
            end
            CD_upsweep=0.075*0.75*sin(sweep_t*pi/180)*obj.d_eff^2*pi/4/S_ref;
            
            obj.CD_f=k*CF*obj.S_wet/S_ref+CD_upsweep;
            obj.D_f=0.5*state.rho_air*norm(state.V_inf)^2*k*CF*obj.S_wet;
            
        end
        
        function obj=read_xml_definition(obj,xmlstruct)
            if strcmp(xmlstruct.tag,'FUSELAGE')
                obj.name=xmlstruct.attribs(1).value;
                for i=1:length(xmlstruct.child)
                    if strcmp(xmlstruct.child(i).tag,'EXTERNALFEM')
                        obj.isExternalFEM = 1;
                        for j=1:length(xmlstruct.child(i).child)
                            if strcmp(xmlstruct.child(i).child(j).tag,'PATHNODECOORDS')
                                obj.pathNodeCoords = xmlstruct.child(i).child(j).value;
                            elseif strcmp(xmlstruct.child(i).child(j).tag,'PATHMASSMATRIX')
                                obj.pathMassMatrix = xmlstruct.child(i).child(j).value;
                            elseif strcmp(xmlstruct.child(i).child(j).tag,'PATHSTIFFNESSMATRIX')
                                obj.pathStiffnessMatrix = xmlstruct.child(i).child(j).value;
                            end
                        end
                    end
                    if strcmp(xmlstruct.child(i).tag,'SEGMENT')
                        if (length(xmlstruct.child(i).child(1).child)==3)
                            obj.fuselage_segments=class_fuselagesegment(xmlstruct.child(i));
                        elseif  strcmp(xmlstruct.child(i).child(1).child.tag,'ATTACH_TO')
                            if strcmp(xmlstruct.child(i).child(1).child.value,'previous')
                                 attach_pos=obj.fuselage_segments(end).get_rear_center();
                                 xmlstruct.child(i).child(1).child(1).tag='X';
                                 xmlstruct.child(i).child(1).child(2).tag='Y';
                                 xmlstruct.child(i).child(1).child(3).tag='Z';
                                 xmlstruct.child(i).child(1).child(1).value=num2str(attach_pos(1));
                                 xmlstruct.child(i).child(1).child(2).value=num2str(attach_pos(2));
                                 xmlstruct.child(i).child(1).child(3).value=num2str(attach_pos(3));
                                 obj.fuselage_segments(i)=class_fuselagesegment(xmlstruct.child(i));
                            end
                        end
                    end 
                end
            else
                fprintf('Unknown Data Format');
            end
            if strcmp(xmlstruct.tag,'NACELLE')
                obj.name=xmlstruct.attribs(1).value;
                for i=1:length(xmlstruct.child)
                    if strcmp(xmlstruct.child(i).tag,'SEGMENT')
                        if (length(xmlstruct.child(i).child(1).child)==3)
                            obj.fuselage_segments=class_fuselagesegment(xmlstruct.child(i));
                        elseif  strcmp(xmlstruct.child(i).child(1).child.tag,'ATTACH_TO')
                            if strcmp(xmlstruct.child(i).child(1).child.value,'previous')
                                 attach_pos=obj.fuselage_segments(end).get_rear_center();
                                 xmlstruct.child(i).child(1).child(1).tag='X';
                                 xmlstruct.child(i).child(1).child(2).tag='Y';
                                 xmlstruct.child(i).child(1).child(3).tag='Z';
                                 xmlstruct.child(i).child(1).child(1).value=num2str(attach_pos(1));
                                 xmlstruct.child(i).child(1).child(2).value=num2str(attach_pos(2));
                                 xmlstruct.child(i).child(1).child(3).value=num2str(attach_pos(3));
                                 obj.fuselage_segments(i)=class_fuselagesegment(xmlstruct.child(i));
                            end
                        end
                    end 
                end
            else
                fprintf('Unknown Data Format');
            end
            
        end
        
    end
    
end
