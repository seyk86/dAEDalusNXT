%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_fuselagesegment
    
    properties
        
        name;
        
        % length
        l;
        % front width
        w_f;
        % rear width
        w_r;
        % front height
        h_f;
        % rear height
        h_r;
        % center line sweep (angle formed from the x axis in x-z plane, signal defined acording to rotation in y axis)
        sweep;
        % surface wetted area
        S_wet;
        
        % center position of the segment's front part
        pos;
        % edges center point coordinates
        xyz;
        % center points for the defined grid
        center_coords;
        % shell width in each node defined in center_coords
        shell_width;
        % shell height in each node defined in center_coords
        shell_height;
        % coordinates of the points that form an elliptical segment
        shell_coords;
        % taper ratio calculated in each node defined in center_coords
        taper;
        % number of circle divisions used for computing shell_coords
        n_circ;
        
        % number of length divisions used for computing grid
        n_length;
        % grid coordinates
        grid;
        % panels edge points
        panels;
        % number of the point representing the segment trailing edge
        te_idx;
        % grid coordinates
        grid_flat;
        % panels edge points
        panels_flat;
        % number of the point representing the segment trailing edge
        te_idx_flat;
        
        % for 2D grid
        grid_start_idx_flat;
        panel_start_idx_flat;
        
        % for 3D grid
        grid_start_idx_vol;
        panel_start_idx_vol;
        
        %> number of beam elements; will be determined by grid spacing if
        % zero
        nBeamelements = 0;
    end
    
    methods
        
        function obj=class_fuselagesegment(varargin)
            if nargin==1
                obj=obj.read_xml_definition(varargin{1});
                obj=obj.compute_shell_coords();
            elseif nargin==5
                obj.name=varargin{1};
                obj.pos=varargin{2};
                obj.l=varargin{3};
                obj.w_f=varargin{4};
                obj.w_r=varargin{4};
                obj.h_f=varargin{5};
                obj.h_r=varargin{5};
                obj.sweep=0;
                obj=obj.compute_segment_coordinates();
                obj=obj.compute_shell_coords();
            else
                obj.l=0;
                obj.w_f=0;
                obj.w_r=0;
                obj.h_f=0;
                obj.h_r=0;
                obj.sweep=0;
                obj.pos=zeros(1,3);
                obj.xyz=zeros(3,2);
            end
        end
        
        function pos=get_rear_center(obj)
            pos=obj.pos+[obj.l 0 -obj.l*sin(obj.sweep*pi/180)];
        end
        
        function obj=compute_segment_coordinates(obj)
            
            p1=obj.pos;
            p2=obj.pos+[obj.l 0 -obj.l*sin(obj.sweep*pi/180)];
            
            obj.xyz=[p1' p2'];
            
        end
        
        function obj=compute_shell_coords(obj,varargin)
            
            obj.center_coords=[];
            obj.shell_coords=[];
            obj.shell_width=[];
            obj.shell_height=[];
            
            switch nargin
                case 1
                    r_min=min([obj.w_f/obj.w_r,obj.w_r/obj.w_f,obj.h_f/obj.h_r,obj.h_r/obj.h_f]);
                    n=floor(sqrt(obj.l/r_min));
                    obj.n_circ=18;
                    isExternalFEM = 0;
                case 2
                    n=varargin{1};
                    obj.n_circ=18;
                    isExternalFEM = 0;
                case 3
                    n=varargin{1};
                    obj.n_circ=varargin{2};
                    isExternalFEM = 0;
                case 4
                    n=varargin{1};
                    iNodes = varargin{2};
                    pathNodeCoords = varargin{3};
                    obj.n_circ=18;
                    isExternalFEM = 1;
            end
            
            if obj.nBeamelements>0
                if obj.nBeamelements<n
                    disp('WARNING: The number of beamelements selected for this fuselage segment results in a lower resolution than defined in the grid spacing!')
                end
                n = obj.nBeamelements;
            end
            
            length_grid=0:1/n:1;
            obj.n_length=n;
            ncirc=obj.n_circ+1;
            aux_points=linspace(0,2*pi,ncirc);
            
            obj.center_coords=zeros(3,n+1);
            obj.shell_coords=zeros(3,ncirc*(n+1));
            
            if isExternalFEM == 0
                for i=1:n+1
                    obj.center_coords(:,i)=obj.xyz(:,1)*(1-length_grid(i))+obj.xyz(:,2)*length_grid(i);
                end
            else
                nodeCoords = importdata(pathNodeCoords);
                obj.center_coords = nodeCoords(iNodes:iNodes+n,:)';
            end
            
            obj.taper=ones(2,n+1);
            x=obj.center_coords(1,:)-ones(1,n+1)*obj.center_coords(1,1);
            if obj.w_f~=obj.w_r
                % taper ratio of the type taper=sqrt((x-c)/a)/w
                a_w=-obj.l/(obj.w_f^2-obj.w_r^2);
                c_w=obj.l*obj.w_f^2/(obj.w_f^2-obj.w_r^2);
                obj.taper(1,:)=((x-ones(1,n+1)*c_w)./(ones(1,n+1)*a_w)).^0.5/obj.w_f;
            end
            if obj.h_f~=obj.h_r
                a_h=-obj.l/(obj.h_f^2-obj.h_r^2);
                c_h=obj.l*obj.h_f^2/(obj.h_f^2-obj.h_r^2);
                obj.taper(2,:)=((x-ones(1,n+1)*c_h)./(ones(1,n+1)*a_h)).^0.5/obj.h_f;
            end
            
            obj.shell_width=obj.w_f*obj.taper(1,:);
            obj.shell_height=obj.h_f*obj.taper(2,:);
            
            for i=1:n+1
                shell_coords(:,1+ncirc*(i-1):ncirc*i)=...
                    [ones(1,ncirc)*obj.center_coords(1,i);...
                    ones(1,ncirc)*obj.center_coords(2,i)+obj.shell_width(i)*cos(aux_points);...
                    ones(1,ncirc)*obj.center_coords(3,i)+obj.shell_height(i)*sin(aux_points)];
            end
            
            
            
            obj.shell_coords=shell_coords;
            
            
            
        end
        
        function obj=compute_grid_flat(obj,z_coord,n_pan_horz,n_pan_vert)
            n=obj.n_length;
            nel=size(obj.center_coords,2);
            nel_horz=0;
            nel_vert=0;
            nflat=5;
            offset=0;
            for i=1:n+1
                rlsg=real(sqrt(1-(z_coord-obj.center_coords(3,i))^2/obj.shell_width(i)^2))*obj.shell_height(i);
                y_coord=rlsg;%+obj.center_coords(2,i);
                
                if rlsg~=0
                    y_spacing=(obj.center_coords(2,i)-y_coord):2*y_coord/(n_pan_horz-1):(y_coord+obj.center_coords(2,i));
                    shell_coords_flat_horz(:,1+n_pan_horz*(i-offset-1):n_pan_horz*(i-offset))=...
                        [ones(1,n_pan_horz)*obj.center_coords(1,i);...
                        y_spacing ;...
                        ones(1,n_pan_horz)*z_coord];
                    nel_horz=nel_horz+1;
                else
                    offset=offset+1;
                end

                
            end
            
            for i=1:n+1
                z_spacing=(obj.center_coords(3,i)-obj.shell_height(i)):2*obj.shell_height(i)/(n_pan_vert-1):(obj.center_coords(3,i)+obj.shell_height(i));
                shell_coords_flat_vert(:,1+n_pan_vert*(i-1):n_pan_vert*i)=...
                    [ones(1,n_pan_vert)*obj.center_coords(1,i);...
                    ones(1,n_pan_vert)*obj.center_coords(2,i);...
                    z_spacing];
                nel_vert=nel_vert+1;
            end
            
            
            
            

            
            grid_flat_horz=zeros(3,n_pan_horz*(nel_horz));
            grid_flat_vert=zeros(3,n_pan_vert*(nel_vert));
            te_idx_flat_horz=zeros(1,n_pan_horz*(nel_horz));
            te_idx_flat_vert=zeros(1,n_pan_vert*(nel_vert));
            % te_idx_flat=zeros(1,nflat*(nel_horz+nel_vert));
            panels_flat_horz=zeros(4,(n_pan_horz-1)*(nel_horz-1));
            panels_flat_vert=zeros(4,(n_pan_vert-1)*(nel_vert-1));
            
            for j=1:n_pan_horz-1
                for i=1:nel_horz-1
                    panels_flat_horz(:,i+(nel_horz-1)*(j-1))=[i+nel_horz*(j-1); i+nel_horz*j; (i+1)+nel_horz*j; (i+1)+nel_horz*(j-1)];
                end
            end
            
            for j=1:n_pan_vert-1
                for i=1:nel_vert-1
                    panels_flat_vert(:,i+(nel_vert-1)*(j-1))=[i+nel_vert*(j-1); i+nel_vert*j; (i+1)+nel_vert*j; (i+1)+nel_vert*(j-1)];
                end
            end
            %
            %             for j=nflat/2+1:nflat-1
            %                 for i=1:nel_vert-1
            %                     panels_flat(:,i+(nel_vert-1)*(j-1))=[i+nel_vert*(j-1); i+nel_vert*j; (i+1)+nel_vert*j; (i+1)+nel_vert*(j-1)];
            %                 end
            %             end
            
            for j=1:n_pan_horz
                for i=1:nel_horz
                    grid_flat_horz(:,i+nel_horz*(j-1))=shell_coords_flat_horz(:,j+n_pan_horz*(i-1));
                end
                %                if(obj.taper(i,1)<1)
                te_idx_flat_horz(1+nel_horz*(j-1):nel_horz*j)=ones(1,nel_horz)*nel_horz*j;
                %      else
                %    te_idx_flat_horz(1+nel_horz*(j-1):nel_horz*j)=ones(1,nel_horz)*nel_horz*j;
                %     end
            end
            
            for j=1:n_pan_vert
                for i=1:nel_vert
                    grid_flat_vert(:,i+nel_vert*(j-1))=shell_coords_flat_vert(:,j+n_pan_vert*(i-1));
                end
                te_idx_flat_vert(1+nel_vert*(j-1):nel_vert*j)=ones(1,nel_vert)*nel_vert*j;
            end
            %
            
            obj.panels_flat=[panels_flat_horz panels_flat_vert+length(grid_flat_horz)];
            obj.grid_flat=[grid_flat_horz grid_flat_vert];
            obj.te_idx_flat=[te_idx_flat_horz te_idx_flat_vert];
        end
        
        
        function obj=compute_grid(obj)
            
            nel=size(obj.center_coords,2);
            ncirc=obj.n_circ+1;
            shell_coords=obj.shell_coords;
            
            grid=zeros(3,ncirc*nel);
            te_idx=zeros(1,ncirc*nel);
            panels=zeros(4,(ncirc-1)*(nel-1));
            
            for j=1:ncirc-1
                for i=1:nel-1
                    panels(:,i+(nel-1)*(j-1))=[i+nel*(j-1); i+nel*j; (i+1)+nel*j; (i+1)+nel*(j-1)];
                end
            end
            
            for j=1:ncirc
                for i=1:nel
                    grid(:,i+nel*(j-1))=shell_coords(:,j+ncirc*(i-1));
                end
                te_idx(1+nel*(j-1):nel*j)=ones(1,nel)*nel*j;
            end
            
            obj.panels=panels;
            obj.grid=grid;
            obj.te_idx=te_idx;
            
        end
        
        function grid_deflected=compute_grid_deflected(obj,panels,grid_deflected,deflections,z_coord,n_pan_horz,n_pan_vert)
            if length(deflections)==6 %% then it's a nacelle
                for i=1:length(obj.center_coords)
                   deflections=[deflections deflections]; 
                end
            end
            
            
            n=obj.n_length;
            nel=size(obj.center_coords,2);
            nel_horz=0;
            nel_vert=0;
            nflat=5;
            offset=0;
            for i=1:n+1
                torsion_angle=deflections(5+6*(i-1));
                rlsg=real(sqrt(1-(z_coord-obj.center_coords(3,i))^2/obj.shell_width(i)^2))*obj.shell_height(i);
                y_coord=rlsg;%+obj.center_coords(2,i);
                
                if rlsg~=0
                    y_spacing=(obj.center_coords(2,i)+deflections(2+6*(i-1))-y_coord):2*y_coord/(n_pan_horz-1):(y_coord+obj.center_coords(2,i)+deflections(2+6*(i-1)));
                    shell_coords_flat_horz(:,1+n_pan_horz*(i-offset-1):n_pan_horz*(i-offset))=...
                        [ones(1,n_pan_horz)*(obj.center_coords(1,i)+deflections(1+6*(i-1)));...
                        y_spacing ;...
                        ones(1,n_pan_horz)*(z_coord+deflections(3+6*(i-1)))];
                    nel_horz=nel_horz+1;
                else
                    offset=offset+1;
                end
                
                rot_matrix=[1 0 0;
                        0 sin(torsion_angle) cos(torsion_angle);
                        0 cos(torsion_angle) sin(torsion_angle)];
            end
            
            for i=1:n+1
                z_spacing=(obj.center_coords(3,i)+deflections(3+6*(i-1))-obj.shell_height(i)):2*obj.shell_height(i)/(n_pan_vert-1):(obj.center_coords(3,i)+deflections(3+6*(i-1))+obj.shell_height(i));
                shell_coords_flat_vert(:,1+n_pan_vert*(i-1):n_pan_vert*i)=...
                    [ones(1,n_pan_vert)*(obj.center_coords(1,i)+deflections(1+6*(i-1)));...
                    ones(1,n_pan_vert)*(obj.center_coords(2,i)+deflections(2+6*(i-1)));...
                    z_spacing];
                nel_vert=nel_vert+1;
            end
            
            grid_flat_horz=zeros(3,n_pan_horz*(nel_horz));
            grid_flat_vert=zeros(3,n_pan_vert*(nel_vert));
 
            for j=1:n_pan_horz
                for i=1:nel_horz
                    grid_flat_horz(:,i+nel_horz*(j-1))=shell_coords_flat_horz(:,j+n_pan_horz*(i-1));
                end
            end
            
            for j=1:n_pan_vert
                for i=1:nel_vert
                    grid_flat_vert(:,i+nel_vert*(j-1))=shell_coords_flat_vert(:,j+n_pan_vert*(i-1));
                end
             
            end  
            
            grid_deflected(:,obj.grid_start_idx_flat:obj.grid_start_idx_flat+length(obj.grid_flat)-1)=[grid_flat_horz grid_flat_vert];     
        end
        
        
        function obj=compute_wetted_area(obj)
            A=0;
            for i=1:length(obj.shell_coords)-(obj.n_circ+1)
                if mod(i,(obj.n_circ+1))
                    v1=obj.shell_coords(:,i+obj.n_circ+1)-obj.shell_coords(:,i);
                    v2=obj.shell_coords(:,i+1)-obj.shell_coords(:,i);
                    A1=norm(cross(v1,v2))/2;
                    
                    v3=obj.shell_coords(:,i+obj.n_circ+1)-obj.shell_coords(:,i+obj.n_circ+2);
                    v4=obj.shell_coords(:,i+1)-obj.shell_coords(:,i+obj.n_circ+2);
                    A2=norm(cross(v3,v4))/2;
                else
                    A1=0;
                    A2=0;
                end
                A=A+A1+A2;
            end
            obj.S_wet=A;
        end
        
        function obj=set_nose(obj,cabin_segment,varargin)
            
            obj.name='Nose';
            if nargin==3
                obj.l=varargin{1};
            else
                % Considering a typical fineness ratio of 1.5 (http://adg.stanford.edu/aa241/fuselayout/fuseplanform.html)
                obj.l=1.5*2*(cabin_segment.w_f*cabin_segment.h_f)^0.5;
            end
            obj.sweep=-7;
            obj.w_f=0.05*cabin_segment.w_f;
            obj.w_r=cabin_segment.w_f;
            obj.h_f=0.05*cabin_segment.h_f;
            obj.h_r=cabin_segment.h_f;
            
            obj.pos=cabin_segment.pos-[obj.l 0 -obj.l*sin(obj.sweep*pi/180)];
            
            obj=obj.compute_segment_coordinates();
            obj=obj.compute_shell_coords();
            
        end
        
        function obj=set_tail(obj,cabin_segment,varargin)
            
            obj.name='Tail';
            if nargin==3
                obj.l=varargin{1};
            else
                % Considering a typical fineness ratio of 2 (http://adg.stanford.edu/aa241/fuselayout/fuseplanform.html)
                obj.l=2*2*(cabin_segment.w_r*cabin_segment.h_r)^0.5;
            end
            %             obj.sweep=-14;
            obj.sweep=-3;
            obj.w_f=cabin_segment.w_r;
            obj.w_r=0.3*cabin_segment.w_r;
            obj.h_f=cabin_segment.h_r;
            obj.h_r=0.3*cabin_segment.h_r;
            
            obj.pos=get_rear_center(cabin_segment);
            
            obj=obj.compute_segment_coordinates();
            obj=obj.compute_shell_coords();
            
        end
        
        function obj=plot_elliptical_segment(obj)
            
            caxis([0 1])
            
            ncirc=obj.n_circ+1;
            
            for i=1:obj.n_length
                
                hold on
                
                coords_front=obj.shell_coords(:,1+ncirc*(i-1):ncirc*i)';
                coords_rear=obj.shell_coords(:,1+ncirc*i:ncirc*(i+1))';
                
                surface([coords_front(:,1) coords_rear(:,1)],[coords_front(:,2) coords_rear(:,2)],[coords_front(:,3) coords_rear(:,3)],0.9*ones(obj.n_circ,2),'FaceAlpha',0.3);
                
            end
            axis equal
            grid on
        end
        
        function obj=plot_grid(obj)
            hold on
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                alpha(handle,0.4);
            end
            axis equal
            axis tight
            grid on
        end
        
        function obj=read_xml_definition(obj,xmlstruct)
            if strcmp(xmlstruct.tag,'SEGMENT')
                obj.name=xmlstruct.attribs(1).value;
                obj.pos(1)=str2double(xmlstruct.child(1).child(1).value);
                obj.pos(2)=str2double(xmlstruct.child(1).child(2).value);
                obj.pos(3)=str2double(xmlstruct.child(1).child(3).value);
                obj.l=str2double(xmlstruct.child(2).value);
                obj.w_f=str2double(xmlstruct.child(3).value);
                obj.w_r=str2double(xmlstruct.child(4).value);
                obj.h_f=str2double(xmlstruct.child(5).value);
                obj.h_r=str2double(xmlstruct.child(6).value);
                obj.sweep=str2double(xmlstruct.child(7).value);
                
                if length(xmlstruct.child) == 8
                    obj.nBeamelements=str2double(xmlstruct.child(8).value);
                end
                
                obj=obj.compute_segment_coordinates();
                
            else
                fprintf('Unknown Data Format');
            end
        end
        
        
    end
    
end
