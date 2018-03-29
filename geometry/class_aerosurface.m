classdef class_aerosurface
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        name;
        isExternalFEM = 0;  % 0 if wing should be selfdesigned in Daedalus, 
                            % 1 if wing is from external FEM model
        pathNodeCoords;     % if wing is from external FEM model, gives path
                            % of node coordinates matlab file. Nodes of 
                            % right half wing from root to tip
        pathMassMatrix;         % if wing is from external FEM model, gives path
                                % of mass matrix matlab file. Mass matrix
                                % of whole wing
        pathStiffnessMatrix;    % if wing is from external FEM model, gives path
                                % of stiffness matrix matlab file. Stiffness matrix
                                % of whole wing
        
        wing_segments;  % array of wing segments
        
        % for 2D grid
        grid_start_idx;
        panel_start_idx;
        
        % for 3D grid
        grid_start_idx_vol;
        panel_start_idx_vol;
        
        
        c4_coords;
        c4_rl_coords;
        c4_forces;
        c4_moments;
        
        wingbox_coords;
        wingbox_rl_coords;
        wingbox_c4;
        wingbox_height;
        air;
        
        wingbox_rl_coords_mid
        c4_forces_structmesh;
        c4_moments_structmesh;
        
        
        beam_forces_structmesh;
        beam_moments_structmesh;
        
        pos;
        symmetric;
        
        S_ref;      %   wing reference area
        S_wet;      %   wetted area
        S;          %   wing area
        AR;         %   aspect ratio
        TR;         %   taper ratio
        c4sweep;    %   sweep of quarter chord line
        Theta;      %   twist
        b;          %   segment span
        Lamda;      %   quarter chord sweep
        
        c_r;
        c_t;
        c_mac;
        
        CD_0;
        D_f;
        CD_f;
        
        Re;
        grid;
        grid_flat;
        is_te;
        
        grid_vol_lo;
        grid_vol_up;
        grid_deflected;
        panels;
        te_idx;
        
        % volume grid
        grid_vol;
        panels_vol;
        is_te_vol;
        opposite_te_vol;
        te_idx_vol;
        
        grid_wake;
        panels_wake;
        
        T;
        T_sort;
        db;
        db_sort;
        
        panel_to_beam_element;
        
        grid_3D=1;
    end
    
    methods
        % constructor
        function obj = class_aerosurface(varargin)
            if nargin==1
                obj=obj.read_xml_definition(varargin{1});
            elseif nargin==19
                pos=varargin{1};
                symmetric=varargin{2};
                dihed=varargin{3};
                LambdaSpec=varargin{4};
                Lambda=varargin{5};
                property1=varargin{6};
                property2=varargin{8};
                property3=varargin{10};
                property4=varargin{12};
                property5=varargin{14};
                property6=varargin{16};
                property7=varargin{18};
                value1=varargin{7};
                value2=varargin{9};
                value3=varargin{11};
                value4=varargin{13};
                value5=varargin{15};
                value6=varargin{17};
                value7=varargin{19};
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
            elseif nargin==21
                pos=varargin{1};
                symmetric=varargin{2};
                dihed=varargin{3};
                LambdaSpec=varargin{4};
                Lambda=varargin{5};
                property1=varargin{6};
                property2=varargin{8};
                property3=varargin{10};
                property4=varargin{12};
                property5=varargin{14};
                property6=varargin{16};
                property7=varargin{18};
                property8=varargin{20};
                value1=varargin{7};
                value2=varargin{9};
                value3=varargin{11};
                value4=varargin{13};
                value5=varargin{15};
                value6=varargin{17};
                value7=varargin{19};
                value8=varargin{21};
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8);
            else
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda);
            end
            
            %             obj.c_mac=cell2mat(obj.wing_segments.c_mac);
            %             obj.S=cell2mat(obj.wing_segments.S);
            obj.symmetric=obj.wing_segments(1).symmetric;
            %             property1=varargin{1};
            %             property2=varargin{3};
            %             property3=varargin{5};
            %
            %
            %             if(strcmp(property1,'SurfaceArea'))&&(strcmp(property2,'AR'))&&(strcmp(property3,'TR'))
            %
            %                 obj.S=varargin{2};
            %                 obj.S_ref=obj.S;
            %                 obj.AR=varargin{4};
            %                 obj.TR=varargin{6};
            %                 obj.S_wet=2*obj.S;
            %                 obj.b=sqrt(obj.S*obj.AR);
            %                 obj.c_r=2*obj.S/((1+obj.TR)*obj.b);
            %                 obj.c_t=obj.TR*obj.c_r;
            %
            %                 obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
            %             end
            
        end
        
        
        function obj =add_segment(obj,symmetric,dihed,LambdaSpec,Lambda,varargin)
            if nargin==19
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec, Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
            elseif nargin==20
                
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                
                new_segment=class_wingsegment(property8,symmetric, dihed,LambdaSpec, Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
                
            elseif nargin==21
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                value8=varargin{16};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8);
            elseif nargin==23
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                property9=varargin{17};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                value8=varargin{16};
                value9=varargin{18};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8,property9,value9);
            end
            obj.wing_segments=[obj.wing_segments new_segment];
            obj.c_mac=0;
            obj.b=0;
            for i=1:length(obj.wing_segments)
                obj.c_mac=obj.c_mac+obj.wing_segments(i).c_mac*obj.wing_segments(i).b;
                obj.b=obj.b+obj.wing_segments(i).b;
            end
            obj.c_mac=obj.c_mac/obj.b;
            obj.S=obj.S+obj.wing_segments(end).S;
        end
        
        function obj =split_segment(obj,segment_idx,f_1,f_2,varargin)
            
            xyz=obj.wing_segments(segment_idx).xyz;
            p1n=xyz(:,1)+f_1*(xyz(:,2)-xyz(:,1));
            p4n=xyz(:,4)+f_1*(xyz(:,3)-xyz(:,4));
            
            p2n=xyz(:,1)+f_2*(xyz(:,2)-xyz(:,1));
            p3n=xyz(:,4)+f_2*(xyz(:,3)-xyz(:,4));
            symm=obj.symmetric;
            if nargin==4
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);
                end
            elseif nargin==5
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s2=s2.add_control_surface(varargin{1});
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);
                    s2=s2.add_control_surface(varargin{1});
                end
            elseif nargin==6
                x=varargin{1};
                y=varargin{2};
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,x,y,profile_r,profile_t);
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,x,y,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);
                end
            end
            if f_2==1
                if segment_idx+1<=length(obj.wing_segments)
                    obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 obj.wing_segments(segment_idx+1:end)];
                else
                    obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 ];
                end
            elseif f_1==0
                obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s2 s3 obj.wing_segments(segment_idx+1:end)];
            else
                obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 s3 obj.wing_segments(segment_idx+1:end)];
            end
            obj.c_mac=0;
            for i=1:length(obj.wing_segments)
                obj.c_mac=obj.c_mac+obj.wing_segments(i).c_mac/length(obj.wing_segments);
            end
        end
        
        function obj = compute_friction_drag(obj,state,S_ref)
            obj.CD_f=0;
            obj.D_f=0;
            for i=1:1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_friction_drag(state,S_ref);
                obj.D_f=obj.D_f+obj.wing_segments(i).D_f;
                obj.CD_f=obj.CD_f+obj.wing_segments(i).CD_f;
            end
        end
        
        function obj= set_cs_deflections(obj,varargin)
            if nargin==4
                idx=varargin{2};
                delta=varargin{3};
                
                for i=1:length(idx)
                    if strcmp(varargin{1},'te')
                        obj.wing_segments(idx(i))=obj.wing_segments(idx(i)).set_cs_deflections('te',delta(i,:));
                    elseif strcmp(varargin{1},'le')
                        obj.wing_segments(idx(i))=obj.wing_segments(idx(i)).set_cs_deflections('le',delta(i,:));
                    end
                end
                
            elseif nargin==7
                idx1=varargin{2};
                delta1=varargin{3};
                
                idx2=varargin{5};
                delta2=varargin{6};
                
                if strcmp(varargin{1},'te') && strcmp(varargin{4},'le')
                    for i=1:length(idx1)
                        obj.wing_segments(idx1(i))=obj.wing_segments(idx1(i)).set_cs_deflections('te',delta1(i,:));
                    end
                    
                    for i=1:length(idx2)
                        obj.wing_segments(idx2(i))=obj.wing_segments(idx2(i)).set_cs_deflections('le',delta2(i,:));
                    end
                    
                elseif strcmp(varargin{1},'le') && strcmp(varargin{4},'te')
                    for i=1:length(idx1)
                        obj.wing_segments(idx1(i))=obj.wing_segments(idx1(i)).set_cs_deflections('le',delta1(i,:));
                    end
                    
                    for i=1:length(idx2)
                        obj.wing_segments(idx2(i))=obj.wing_segments(idx2(i)).set_cs_deflections('te',delta2(i,:));
                    end
                end
            end
        end
        
        function obj = plot(obj)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).plot_segment();
                %obj.wing_segments(i).compute_grid();
                hold on
            end
        end
        
        function obj=update_idx(obj)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).grid_start_idx=obj.wing_segments(i).grid_start_idx+(obj.grid_start_idx-1);
                obj.wing_segments(i).panel_start_idx=obj.wing_segments(i).panel_start_idx+(obj.panel_start_idx-1);
            end
        end
        
        function obj=compute_bspline_grid(obj,x_max,y_max)
            
            
        end
        
        function obj=compute_grid(obj,x_max,y_max,wake)
            
            grid=[];
            grid_flat=[];
            panels=[];
            grid_len_b4=[];
            grid_vol_up=[];
            grid_vol_lo=[];
            is_te=[];
            
            te_idx=[];
            
            if (wake==1)||(wake==2)
                grid_wake=[];
                panels_wake=[];
            end
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).left_control_surfaces();
            end
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_grid(x_max,y_max,wake);
                grid_len_b4=length(grid);
                grid=[grid obj.wing_segments(i).grid];
                grid_flat=[grid_flat obj.wing_segments(i).grid_flat];
                
                grid_vol_up=[grid_vol_up obj.wing_segments(i).grid_vol_upper];
                grid_vol_lo=[grid_vol_lo obj.wing_segments(i).grid_vol_lower];
                
                
                is_te=[is_te obj.wing_segments(i).is_te];
                
                panel_len_b4=length(panels);
                te_idx=[te_idx, obj.wing_segments(i).te_idx+grid_len_b4];
                panels=[panels, obj.wing_segments(i).panels+grid_len_b4];
                panel_len_aft=length(panels);
                obj.wing_segments(i).grid_start_idx=grid_len_b4+1;
                obj.wing_segments(i).panel_start_idx=panel_len_b4+1;
                
                if (wake==1)||(wake==2)
                    grid_wake_len_b4=length(grid_wake);
                    grid_wake=[grid_wake obj.wing_segments(i).grid_wake];
                    panels_wake=[panels_wake obj.wing_segments(i).panels_wake+grid_wake_len_b4];
                end 
            end
            
            if obj.symmetric==1
                grid_mirror=[];
                grid_flat_mirror=[];
                
                grid_vol_mirror_up=[];
                grid_vol_mirror_lo=[];
                is_te_mirror=[];
                
                if (wake==1)||(wake==2)
                    grid_wake_mirror=[];
                    panels_wake_mirror=[];
                end
                
                for i=1:length(obj.wing_segments)
                    obj.wing_segments(i)=obj.wing_segments(i).right_control_surfaces();
                    obj.wing_segments(i)=obj.wing_segments(i).compute_grid(x_max,y_max,wake);
                    grid_mirror=[grid_mirror obj.wing_segments(i).grid];
                    grid_flat_mirror=[grid_flat_mirror obj.wing_segments(i).grid_flat];
                    
                    is_te_mirror=[is_te_mirror obj.wing_segments(i).is_te];
                    grid_vol_mirror_up=[grid_vol_mirror_up obj.wing_segments(i).grid_vol_upper];
                    grid_vol_mirror_lo=[grid_vol_mirror_lo obj.wing_segments(i).grid_vol_lower];
                    
                    if (wake==1)||(wake==2)
                        grid_wake_mirror=[grid_wake_mirror obj.wing_segments(i).grid_wake];
                    end
                end
                
                grid_len_b4=length(grid_mirror);
                obj.grid=[grid [grid_mirror(1,:);-grid_mirror(2,:);grid_mirror(3,:)]];
                obj.grid_flat=[grid_flat [grid_flat_mirror(1,:);-grid_flat_mirror(2,:);grid_flat_mirror(3,:)]];

                obj.is_te=[is_te is_te_mirror];
                obj.grid_vol_up=[grid_vol_up [grid_vol_mirror_up(1,:);-grid_vol_mirror_up(2,:);grid_vol_mirror_up(3,:)]];
                obj.grid_vol_lo=[grid_vol_lo [grid_vol_mirror_lo(1,:);-grid_vol_mirror_lo(2,:);grid_vol_mirror_lo(3,:)]];
                if (wake==1)||(wake==2)
                    grid_wake_len_b4=length(grid_wake_mirror);
                    obj.grid_wake=[grid_wake [grid_wake_mirror(1,:);-grid_wake_mirror(2,:);grid_wake_mirror(3,:)]];
                end
                
                obj.panels=[panels,[panels(2,:);panels(1,:);panels(4,:);panels(3,:)]+grid_len_b4];
                
                if (wake==1)||(wake==2)
                   obj.panels_wake=[panels_wake,[panels_wake(2,:);panels_wake(1,:);panels_wake(4,:);panels_wake(3,:)]+grid_wake_len_b4]; 
                end
                
                obj.te_idx=[te_idx,te_idx+grid_len_b4];
%                 for i=1:length(obj.wing_segments)
%                     obj.wing_segments(i)=obj.wing_segments(i).mirror_control_surfaces();
%                 end
                for i=1:length(obj.wing_segments)
                    obj.wing_segments(i)=obj.wing_segments(i).left_control_surfaces();
                end
            else
                
                obj.grid=grid;
                obj.grid_flat=grid_flat;
                obj.is_te=is_te;
                obj.grid_vol_up=grid_vol_up;
                obj.grid_vol_lo=grid_vol_lo;
                
                obj.panels=panels;
                obj.te_idx=te_idx;
                if (wake==1)||(wake==2)
                    obj.grid_wake=grid_wake;
                    obj.panels_wake=panels_wake;
                end
            end
            
            obj.grid_vol=[obj.grid_vol_up obj.grid_vol_lo];
            inv_pan(1,:)=obj.panels(4,:);
            inv_pan(2,:)=obj.panels(3,:);
            inv_pan(4,:)=obj.panels(1,:);
            inv_pan(3,:)=obj.panels(2,:);
            obj.panels_vol=[obj.panels inv_pan+length(obj.grid_vol_up)]; 
            
            obj.is_te_vol=[obj.is_te obj.is_te*0];
            obj.te_idx_vol=[obj.te_idx obj.te_idx+length(obj.grid)];
            op_idx=1:1:length(obj.is_te);
            op_idx=op_idx+length(obj.is_te);
            obj.opposite_te_vol=[op_idx obj.is_te*0];
        end
        
        function panel_to_beam_element=compute_force_interpolation_matrix(obj,panel_to_beam_element,panels,grid)
            % for debugging
            observer=1;
            
            % for debugging
            if observer==1
                figure
                % obj.plot_grid
                hold on
            end
            
            % start with beam index 1
            beam_idx=1;
            
            % do for all wing segments
            for k=1:length(obj.wing_segments)
                next_segment=0;
                prev_segment=0;
                % do for all beam elements in this wing segment
                for kk=1:length(obj.wing_segments(k).wingbox_coords(1,:,1))-1
                    %get current beamelement points; inner: B1; outer: B2
                    B1=0.5*obj.wing_segments(k).wingbox_coords(:,kk,1)+0.5*obj.wing_segments(k).wingbox_coords(:,kk,2);
                    B2=0.5*obj.wing_segments(k).wingbox_coords(:,kk+1,1)+0.5*obj.wing_segments(k).wingbox_coords(:,kk+1,2);
                    %get 4 corner points of wing segment, 1:front left
                    %2:front right, 3: rear right, 4: rear left
                    P1=obj.wing_segments(k).xyz(:,1);
                    P2=obj.wing_segments(k).xyz(:,2);
                    P3=obj.wing_segments(k).xyz(:,3);
                    P4=obj.wing_segments(k).xyz(:,4);
                    
                    segment=[P1 P2 P3 P4];
                    beamelement=[B1 B2];
                    
                    %% U1
                    % U1 is the intersection point of the leading edge of
                    % the wingsegment and a plane normal to the
                    % beamelement through B1
                    
                    % beam node 1
                    
                    diff_1=B2-B1;
                    
                    [nono,comp_idx_2]=max(diff_1);  %comp_idx_2 = Main Direction
                    diff_1(comp_idx_2)=0;
                    [nono,comp_idx_1]=max(diff_1);  %comp_idx_1 = Second Main Direction
                    if comp_idx_1==comp_idx_2
                        if comp_idx_1==3
                            comp_idx_2=1;
                        else
                            comp_idx_1=comp_idx_2+1;
                        end
                    end
                    n1=cross((B1-P2),(B1-P1));  % normal vector of plane through leading edge of segment and B1 
                    n1=n1/norm(n1);
                    rb1=cross((B2-B1),n1); %normal vector to beam element in plane normal to n1
                    rb1=rb1/norm(rb1);
                    rp=P2-P1;
                    %lambda=(B1(1)-P1(1)-B1(2)*rp(1)/rp(2)+P1(2)*rp(1)/rp(2))/(rb1(2)*rp(1)/rp(2)-rb1(1));
                    if (((P1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2))+(B1(comp_idx_2)-P1(comp_idx_2))*(P2(comp_idx_1)-P1(comp_idx_1))-B1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2)))~=0) && ((rb1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2))-rb1(comp_idx_2)*(P2(comp_idx_1)-P1(comp_idx_1)))~=0))
                        lambda=(P1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2))+(B1(comp_idx_2)-P1(comp_idx_2))*(P2(comp_idx_1)-P1(comp_idx_1))-B1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2)))/(rb1(comp_idx_1)*(P2(comp_idx_2)-P1(comp_idx_2))-rb1(comp_idx_2)*(P2(comp_idx_1)-P1(comp_idx_1)));
                    else
                        lambda=0;
                    end
                    U1=B1+lambda*rb1;
                    %Correction S.Binder
                    %works better when the intersection point is computed between a plane normal
                    %to the beam element and through B1 and the leading edge
                    U1=P1+(dot((B1-P1),(B2-B1))/dot((P2-P1),(B2-B1)))*((P2-P1));
                    
                    
                    % U12=B1+lambda*rb1;
                    search_idx=2;
                    if P2(2)-P1(2)>0
                        search_idx=2;
                    elseif P2(3)-P1(3)>0
                        search_idx=3;
                    elseif P2(1)-P1(1)>0
                        search_idx=1;
                    end
                    lambda=(U1(search_idx)-P1(search_idx))/(P2(search_idx)-P1(search_idx));
                    errU1_online=norm(U1-P1-lambda*(P2-P1));
                    U1=P1+lambda*(P2-P1);
                    if lambda>1+100*eps
                        next_segment=1;
                    elseif lambda<0.0-100*eps
                        prev_segment=1;
                    end
                    %% U2
                    % U2 is the intersection point of the leading edge of
                    % the wingsegment and a plane normal to the
                    % beamelement through B2
                    
                    n1=cross((B1-P3),(B1-P4));
                    n1=n1/norm(n1);
                    rb1=cross((B2-B1),n1);
                    rb1=rb1/norm(rb1);
                    rp=P4-P3;
                    %lambda=(B1(1)-P4(1)-B1(2)*rp(1)/rp(2)+P4(2)*rp(1)/rp(2))/(rb1(2)*rp(1)/rp(2)-rb1(1));
                    
                    if (((P4(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2))+(B1(comp_idx_2)-P4(comp_idx_2))*(P3(comp_idx_1)-P4(comp_idx_1))-B1(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2)))~=0) &&((rb1(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2))-rb1(comp_idx_2)*(P3(comp_idx_1)-P4(comp_idx_1)))~=0))
                        lambda=(P4(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2))+(B1(comp_idx_2)-P4(comp_idx_2))*(P3(comp_idx_1)-P4(comp_idx_1))-B1(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2)))/(rb1(comp_idx_1)*(P3(comp_idx_2)-P4(comp_idx_2))-rb1(comp_idx_2)*(P3(comp_idx_1)-P4(comp_idx_1)));
                    else
                        lambda=0;
                    end
                    
                    U2=B1+lambda*rb1;
                    
                    %Correction S.Binder
                    %works better when the intersection point is computed between a plane normal
                    %to the beam element and through B2 and the leading edge
                    U2=P3+(dot((B1-P3),(B2-B1))/dot((P3-P4),(B2-B1)))*((P3-P4));
                    %U22=B1+lambda2*rb1;
                    errU2=norm(U2-B1-lambda*rb1);
                    search_idx=2;
                    if P2(2)-P1(2)>0
                        search_idx=2;
                    elseif P2(3)-P1(3)>0
                        search_idx=3;
                    elseif P2(1)-P1(1)>0
                        search_idx=1;
                    end
                    lambda=(U2(search_idx)-P4(search_idx))/(P3(search_idx)-P4(search_idx));
                    errU2_online=norm(U2-P4-lambda*(P3-P4));
                    U2=P4+lambda*(P3-P4);
                    if lambda>1+100*eps
                        next_segment=1;
                    elseif lambda<0.0-100*eps
                        prev_segment=1;
                    end
                    %% U3
                    % U3 is the intersection point of the trailing edge of
                    % the wingsegment and a plane normal to the
                    % beamelement through B2
                    
                    n2=cross((B2-P2),(B2-P1));
                    n2=n2/norm(n2);
                    rb2=cross((B2-B1),n2);
                    rb2=rb2/norm(rb2);
                    rp=P2-P1;
                    
                    if (((B2(comp_idx_1)-P1(comp_idx_1)-B2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)+P1(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2))~=0) && ((rb2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)-rb2(comp_idx_1))~=0))
                        lambda=(B2(comp_idx_1)-P1(comp_idx_1)-B2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)+P1(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2))/(rb2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)-rb2(comp_idx_1));
                    else
                        lambda=0;
                    end
                    U3=B2+lambda*rb2;
                    
                    %Correction S.Binder
                    %works better when the intersection point is computed between a plane normal
                    %to the beam element and through B2 and the trailing edge
                    U3=P1+(dot((B2-P1),(B2-B1))/dot((P2-P1),(B2-B1)))*((P2-P1));
                    errU3=norm(U3-B2-lambda*rb2);
                    search_idx=2;
                    if P2(2)-P1(2)>0
                        search_idx=2;
                    elseif P2(3)-P1(3)>0
                        search_idx=3;
                    elseif P2(1)-P1(1)>0
                        search_idx=1;
                    end
                    
                    lambda=(U3(search_idx)-P1(search_idx))/(P2(search_idx)-P1(search_idx));
                    U3=P1+lambda*(P2-P1);
                    errU3_online=norm(U3-P1-lambda*(P2-P1));
                    if lambda>1+100*eps
                        next_segment=1;
                    elseif lambda<0.0-100*eps
                        prev_segment=1;
                    end
                    
                    %% U4
                    % U4 is the intersection point of the trailing edge of
                    % the wingsegment and a plane normal to the
                    % beamelement through B1
                    
                    n2=cross((B2-P3),(B2-P4));
                    n2=n2/norm(n2);
                    rb2=cross((B2-B1),n2);
                    rb2=rb2/norm(rb2);
                    rp=P4-P3;
                    if (((B2(comp_idx_1)-P4(comp_idx_1)-B2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)+P4(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2))~=0) && ((rb2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)-rb2(comp_idx_1))~=0))
                        lambda=(B2(comp_idx_1)-P4(comp_idx_1)-B2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)+P4(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2))/(rb2(comp_idx_2)*rp(comp_idx_1)/rp(comp_idx_2)-rb2(comp_idx_1));
                    else
                        lambda=0;
                    end
                        
                    U4=B2+lambda*rb2;
                    
                    %Correction S.Binder
                    %works better when the intersection point is computed between a plane normal
                    %to the beam element and through B1 and the trailing edge
                    U4=P3+(dot((B2-P3),(B2-B1))/dot((P3-P4),(B2-B1)))*((P3-P4));
                    %                     errU4=norm(U4-B2-lambda*rb2);
                    %                     if errU4>1E-10
                    %                         disp('geometry error');
                    %                     end
                    % check if point is in between P4 and P3
                    search_idx=2;
                    if P2(2)-P1(2)>0
                        search_idx=2;
                    elseif P2(3)-P1(3)>0
                        search_idx=3;
                    elseif P2(1)-P1(1)>0
                        search_idx=1;
                    end
                    lambda=(U4(search_idx)-P4(search_idx))/(P3(search_idx)-P4(search_idx));
                    U4=P4+lambda*(P3-P4);
                    errU4_online=norm(U4-P4-lambda*(P3-P4));
                    if lambda>1+100*eps
                        next_segment=1;
                    elseif lambda<0.0-100*eps
                        prev_segment=1;
                    end
                    
%                     obj.plot_grid
%                     fill3(segment(1,:),segment(2,:),segment(3,:),0.4);
%                     axis equal
%                     
%                     
%                                         plot3(beamelement(1,:),beamelement(2,:),beamelement(3,:),'r')
%                                         plot3([beamelement(1,1) beamelement(1,1)+n1(1)],[beamelement(2,1) beamelement(2,1)+n1(2)],[beamelement(3,1) beamelement(3,1)+n1(3)],'r')
%                                         plot3([beamelement(1,1) beamelement(1,1)+rb1(1)],[beamelement(2,1) beamelement(2,1)+rb1(2)],[beamelement(3,1) beamelement(3,1)+rb1(3)],'r')
%                     
%                                         plot3([beamelement(1,2) beamelement(1,2)+n2(1)],[beamelement(2,2) beamelement(2,2)+n2(2)],[beamelement(3,2) beamelement(3,2)+n2(3)],'r')
%                                         plot3([beamelement(1,2) beamelement(1,2)+rb2(1)],[beamelement(2,2) beamelement(2,2)+rb2(2)],[beamelement(3,2) beamelement(3,2)+rb2(3)],'r')
%                     
%                                         plot3(U1(1),U1(2),U1(3),'bx')
%                                         plot3(U2(1),U2(2),U2(3),'bx')
%                                         plot3(U3(1),U3(2),U3(3),'bx')
%                                         plot3(U4(1),U4(2),U4(3),'bx')
                    
                    block=[U1 U3 U4 U2];
                    
                    if next_segment==1
                        if k+1<=length(obj.wing_segments)
                            lim_block=cut_line_quad(obj.wing_segments(k).xyz(:,2:3),block);
                            if size(lim_block,2)>=2
                                
                                lp1=norm(lim_block(:,1)-obj.wing_segments(k).xyz(:,2))/norm(obj.wing_segments(k).xyz(:,3)-obj.wing_segments(k).xyz(:,2));
                                lp2=norm(lim_block(:,2)-obj.wing_segments(k).xyz(:,2))/norm(obj.wing_segments(k).xyz(:,3)-obj.wing_segments(k).xyz(:,2));
                                
                                if lp1>0.01
                                    P1_next=obj.wing_segments(k+1).xyz(:,1)+lp1*(obj.wing_segments(k+1).xyz(:,4)-obj.wing_segments(k+1).xyz(:,1));
                                    P2_next=obj.wing_segments(k+1).xyz(:,1)+lp2*(obj.wing_segments(k+1).xyz(:,4)-obj.wing_segments(k+1).xyz(:,1));
                                    
                                    n_next=cross((P2_next-obj.wing_segments(k+1).xyz(:,1)),(obj.wing_segments(k+1).xyz(:,2)-P2_next));
                                    n_next=n_next/norm(n_next);
                                    rb_next=cross((B2-B1),n_next);
                                    rb_next=rb_next/norm(rb_next);
                                    rp_next=obj.wing_segments(k+1).xyz(:,2)-obj.wing_segments(k+1).xyz(:,1);
                                    lambda=(P2_next(1)-obj.wing_segments(k+1).xyz(1,1)-P2_next(2)*rp_next(1)/rp_next(2)+obj.wing_segments(k+1).xyz(2,1)*rp_next(1)/rp_next(2))/(rb_next(2)*rp_next(1)/rp_next(2)-rb_next(1));
                                    P3_next=P2_next+lambda*rb_next;
                                    lambda=(P1_next(1)-obj.wing_segments(k+1).xyz(1,1)-P1_next(2)*rp_next(1)/rp_next(2)+obj.wing_segments(k+1).xyz(2,1)*rp_next(1)/rp_next(2))/(rb_next(2)*rp_next(1)/rp_next(2)-rb_next(1));
                                    P4_next=P1_next+lambda*rb_next;
                                    block_next=[P1_next P2_next P3_next P4_next];
                                else
                                    P1_next=obj.wing_segments(k+1).xyz(:,1);
                                    P2_next=obj.wing_segments(k+1).xyz(:,1)+lp2*(obj.wing_segments(k+1).xyz(:,4)-obj.wing_segments(k+1).xyz(:,1));
                                    n_next=cross((P2_next-obj.wing_segments(k+1).xyz(:,1)),(obj.wing_segments(k+1).xyz(:,2)-P2_next));
                                    n_next=n_next/norm(n_next);
                                    rb_next=cross((B2-B1),n_next);
                                    rb_next=rb_next/norm(rb_next);
                                    rp_next=obj.wing_segments(k+1).xyz(:,2)-obj.wing_segments(k+1).xyz(:,1);
                                    lambda=(P2_next(1)-obj.wing_segments(k+1).xyz(1,1)-P2_next(2)*rp_next(1)/rp_next(2)+obj.wing_segments(k+1).xyz(2,1)*rp_next(1)/rp_next(2))/(rb_next(2)*rp_next(1)/rp_next(2)-rb_next(1));
                                    P3_next=P2_next+lambda*rb_next;
                                    block_next=[P1_next P2_next P3_next];
                                end
                                %TODO CHECK if comp idx is needed for main direction in initialize_beam_to_panel_for_block
                                panel_to_beam_element=obj.initialize_beam_to_panel_for_block(beam_idx,k+1,block_next,panels,grid,panel_to_beam_element,0,comp_idx_2);
                                h=fill3(block_next(1,:),block_next(2,:),block_next(3,:),'r');
                                set(h,'facealpha',.25);
                            end
                        end
                        next_segment=0;
                    end
                    
                    if prev_segment==1
                        if k>=2
                            lim_block=cut_line_quad(obj.wing_segments(k).xyz(:,[1 4]),block);
                            
                            if size(lim_block,2)>=2
                                lp1=norm(lim_block(:,2)-obj.wing_segments(k).xyz(:,1))/norm(obj.wing_segments(k).xyz(:,4)-obj.wing_segments(k).xyz(:,1));
                                lp2=norm(lim_block(:,1)-obj.wing_segments(k).xyz(:,1))/norm(obj.wing_segments(k).xyz(:,4)-obj.wing_segments(k).xyz(:,1));
                                
                                if lp2<0.99
                                    P1_prev=obj.wing_segments(k-1).xyz(:,2)+lp1*(obj.wing_segments(k-1).xyz(:,3)-obj.wing_segments(k-1).xyz(:,2));
                                    P2_prev=obj.wing_segments(k-1).xyz(:,2)+lp2*(obj.wing_segments(k-1).xyz(:,3)-obj.wing_segments(k-1).xyz(:,2));
                                    n_prev=cross((P1_prev-obj.wing_segments(k-1).xyz(:,3)),(obj.wing_segments(k-1).xyz(:,4)-P1_prev));
                                    n_prev=n_prev/norm(n_prev);
                                    rb_prev=cross((B2-B1),n_prev);
                                    rb_prev=rb_prev/norm(rb_prev);
                                    rp_prev=obj.wing_segments(k-1).xyz(:,4)-obj.wing_segments(k-1).xyz(:,3);
                                    lambda=(P2_prev(1)-obj.wing_segments(k-1).xyz(1,3)-P2_prev(2)*rp_prev(1)/rp_prev(2)+obj.wing_segments(k-1).xyz(2,3)*rp_prev(1)/rp_prev(2))/(rb_prev(2)*rp_prev(1)/rp_prev(2)-rb_prev(1));
                                    P3_prev=P2_prev+lambda*rb_prev;
                                    lambda=(P1_prev(1)-obj.wing_segments(k-1).xyz(1,3)-P1_prev(2)*rp_prev(1)/rp_prev(2)+obj.wing_segments(k-1).xyz(2,3)*rp_prev(1)/rp_prev(2))/(rb_prev(2)*rp_prev(1)/rp_prev(2)-rb_prev(1));
                                    P4_prev=P1_prev+lambda*rb_prev;
                                    block_prev=[P1_prev P2_prev P3_prev P4_prev];
                                else
                                    P1_prev=obj.wing_segments(k-1).xyz(:,2)+lp1*(obj.wing_segments(k-1).xyz(:,3)-obj.wing_segments(k-1).xyz(:,2));
                                    P2_prev=obj.wing_segments(k-1).xyz(:,3);
                                    n_prev=cross((P1_prev-obj.wing_segments(k-1).xyz(:,3)),(obj.wing_segments(k-1).xyz(:,4)-P1_prev));
                                    n_prev=n_prev/norm(n_prev);
                                    rb_prev=cross((B2-B1),n_prev);
                                    rb_prev=rb_prev/norm(rb_prev);
                                    rp_prev=obj.wing_segments(k-1).xyz(:,4)-obj.wing_segments(k-1).xyz(:,3);
                                    lambda=(P1_prev(1)-obj.wing_segments(k-1).xyz(1,3)-P1_prev(2)*rp_prev(1)/rp_prev(2)+obj.wing_segments(k-1).xyz(2,3)*rp_prev(1)/rp_prev(2))/(rb_prev(2)*rp_prev(1)/rp_prev(2)-rb_prev(1));
                                    P3_prev=P1_prev+lambda*rb_prev;
                                    
                                    block_prev=[P1_prev P2_prev P3_prev];
                                end
                                %TODO CHECK if comp idx is needed for main direction in initialize_beam_to_panel_for_block
                                panel_to_beam_element=obj.initialize_beam_to_panel_for_block(beam_idx,k-1,block_prev,panels,grid,panel_to_beam_element,0,comp_idx_2);
                                h=fill3(block_prev(1,:),block_prev(2,:),block_prev(3,:),'r');
                                set(h,'facealpha',.25);
                            end
                        end
                        prev_segment=0;
                    end
                    %fill3(block(1,:),block(2,:),block(3,:),0.4);
                    
                    %TODO CHECK if comp idx is needed for main direction in initialize_beam_to_panel_for_block
                    panel_to_beam_element=obj.initialize_beam_to_panel_for_block(beam_idx,k,block,panels,grid,panel_to_beam_element,0,comp_idx_2);
                    beam_idx=beam_idx+1;
                end
            end

            
  
            for j=obj.panel_start_idx:1:size(obj.panels,2)
                sum=0;
                for i=2:5:size(panel_to_beam_element,2)
                    sum=sum+panel_to_beam_element(j,i);
                end
                if sum==0 && panel_to_beam_element(j,1)~=0
                    disp('WARNING: panel has zero area')
                    panel_to_beam_element(j,2)=1;
                    if j>1
                        panel_to_beam_element(j,1)=panel_to_beam_element(j-1,1);
                    end
                end
            end
            
            for j=1:size(panels,2)
                sum=0;
                for i=2:5:size(panel_to_beam_element,2)
                     sum=sum+panel_to_beam_element(j,i);
                end
                if sum>eps
                    if sum<0.98
                        h=fill3(grid(1,panels(:,j)),grid(2,panels(:,j)),grid(3,panels(:,j)),'y');
                        set(h,'facealpha',.25);
                    end
                    %                      beam_splits=0;
                    for i=2:5:size(panel_to_beam_element,2)
                        %                         if obj.panel_to_beam_element(j,i)~=0
                        %                             beam_splits=beam_splits+1;
                        %                         end
                        %                      end
                        %                      for k=1:beam_splits
                        panel_to_beam_element(j,i)=panel_to_beam_element(j,i)/sum;
                    end
                end
            end
            
                %plotting legend
                [~,h1]=legend('1','2','3','4');
               
                set(h1(6), 'facea', 0.4)
                set(h1(7), 'facea', 0.25)
                set(h1(7), 'faceColor', 'y')
                set(h1(8), 'facea', 0.25)
                set(h1(8), 'faceColor', 'r')
                
            if obj.symmetric==1
                n_beam=beam_idx-2;
                obj.panel_start_idx;
                offset=length(obj.panels)/2;
                for j=1:size(obj.panels,2)/2
                    panel_to_beam_element(j+offset+obj.panel_start_idx-1,:)=panel_to_beam_element(j+obj.panel_start_idx-1,:);
                    panel_to_beam_element(j+offset+obj.panel_start_idx-1,[4 9 14])=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,[4 9 14]);
                    for i=1:5:size(panel_to_beam_element,2)
                        if not(panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)==0)
                            panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)+beam_idx;
%                             panel_to_beam_element(j+offset+obj.panel_start_idx-1,i+2)=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,i+2);
%                             panel_to_beam_element(j+obj.panel_start_idx-1,i+2)=-panel_to_beam_element(j+obj.panel_start_idx-1,i+2);
                        end
                    end
                end   
                for j=1:size(obj.panels,2)/2
                    for i=1:5:size(panel_to_beam_element,2)
                        if not(panel_to_beam_element(j+obj.panel_start_idx-1,i)==0)
                            com_idx=panel_to_beam_element(j+obj.panel_start_idx-1,i)+beam_idx-1;%+(n_beam-2*(panel_to_beam_element(j+obj.panel_start_idx-1,i)-1));
                            panel_to_beam_element(j+obj.panel_start_idx-1,i)=com_idx; 
                        end
                    end
                end
            end  
        end
        
        function obj=T_matrix(obj,panel_to_beam_element,grid,panels,beam)
            
            if obj.symmetric==1
                f=2;
            else
                f=1;
            end
           %  obj.T=zeros(length(obj.wingbox_coords(1,:,1))*6*f-6,size(panel_to_beam_element,1));
           %  obj.db=zeros(length(obj.wingbox_coords(1,:,1))*6*f-6,size(panel_to_beam_element,1));
            obj.T=zeros(length(beam.beamelement)*6+6,size(panel_to_beam_element,1));
            obj.db=zeros(length(beam.beamelement)*6+6,size(panel_to_beam_element,1));
            for j=obj.panel_start_idx:1:(obj.panel_start_idx+size(obj.panels,2))-1
                for i=1:5:size(panel_to_beam_element,2)
                    if not(panel_to_beam_element(j,i)==0)
                        diff=norm(panel_to_beam_element(j,i+2:i+4));
                        %                       idx=j-obj.panel_start_idx+1;
                        %                       r=(0.75*obj.grid(:,obj.panels(1,idx))+0.25*obj.grid(:,obj.panels(4,idx)))-(0.75*obj.grid(:,obj.panels(2,idx))+0.25*obj.grid(:,obj.panels(3,idx)));
                        %phi=180/pi*beam.beamelement(panel_to_beam_element(j,i)).phi;
                        r=(0.75*grid(:,panels(1,j))+0.25*grid(:,panels(4,j)))-(0.75*grid(:,panels(2,j))+0.25*grid(:,panels(3,j)));
                        
                        s=abs(r(2));
                        Lz=eye(3,3);%beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3);
                        
                        F=[0;
                           0;
                           -1/2*panel_to_beam_element(j,i+1)]*s;
                       
                        Ft=Lz'*F;
                        
                        
                        Lz=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3);
                        M=[-1/24*panel_to_beam_element(j,i+1)*s^2;
                            -diff*1/2*panel_to_beam_element(j,i+1)*s;
                            0];

                        Mt1=Lz'*M;
                        
                        M=[+1/24*panel_to_beam_element(j,i+1)*s^2;
                            -diff*1/2*panel_to_beam_element(j,i+1)*s;
                            0];
                        
                        Mt2=Lz'*M;

                        %beamlength=norm(obj.wingbox_coords(:,panel_to_beam_element(j,i)+1)-obj.wingbox_coords(:,panel_to_beam_element(j,i)));
                        obj.T((panel_to_beam_element(j,i)-1)*6+1,j)=obj.T((panel_to_beam_element(j,i)-1)*6+1,j)+Ft(1);
                        obj.T((panel_to_beam_element(j,i)-1)*6+2,j)=obj.T((panel_to_beam_element(j,i)-1)*6+2,j)+Ft(2);
                        obj.T((panel_to_beam_element(j,i)-1)*6+3,j)=obj.T((panel_to_beam_element(j,i)-1)*6+3,j)+Ft(3);
                        obj.T((panel_to_beam_element(j,i)-1)*6+4,j)=obj.T((panel_to_beam_element(j,i)-1)*6+4,j)+Mt1(1);
                        obj.T((panel_to_beam_element(j,i)-1)*6+5,j)=obj.T((panel_to_beam_element(j,i)-1)*6+5,j)+Mt1(2);
                        obj.T((panel_to_beam_element(j,i)-1)*6+6,j)=obj.T((panel_to_beam_element(j,i)-1)*6+6,j)+Mt1(3);    

                        obj.T((panel_to_beam_element(j,i))*6+1,j)=obj.T((panel_to_beam_element(j,i))*6+1,j)+Ft(1);
                        obj.T((panel_to_beam_element(j,i))*6+2,j)=obj.T((panel_to_beam_element(j,i))*6+2,j)+Ft(2);
                        obj.T((panel_to_beam_element(j,i))*6+3,j)=obj.T((panel_to_beam_element(j,i))*6+3,j)+Ft(3);
                        obj.T((panel_to_beam_element(j,i))*6+4,j)=obj.T((panel_to_beam_element(j,i))*6+4,j)+Mt2(1);
                        obj.T((panel_to_beam_element(j,i))*6+5,j)=obj.T((panel_to_beam_element(j,i))*6+5,j)+Mt2(2);
                        obj.T((panel_to_beam_element(j,i))*6+6,j)=obj.T((panel_to_beam_element(j,i))*6+6,j)+Mt2(3);
                        
                        obj.db((panel_to_beam_element(j,i)-1)*6+1,j)=0;
                        obj.db((panel_to_beam_element(j,i)-1)*6+2,j)=0;
                        obj.db((panel_to_beam_element(j,i)-1)*6+3,j)=0;

                        dind=Lz'*[0  -1/2*panel_to_beam_element(j,i+1)*4*pi  0]';
                        %  obj.db((panel_to_beam_element(j,i)-1)*6+4,j)=obj.db((panel_to_beam_element(j,i)-1)*6+4,j)-1/7*panel_to_beam_element(j,i+1)*4*pi*sind(phi)*cosd(phi);
                        %  obj.db((panel_to_beam_element(j,i)-1)*6+5,j)=obj.db((panel_to_beam_element(j,i)-1)*6+5,j)+1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*cosd(phi)*cosd(dihed);
                        obj.db((panel_to_beam_element(j,i)-1)*6+4,j)=obj.db((panel_to_beam_element(j,i)-1)*6+4,j)+dind(1);%-1/2*panel_to_beam_element(j,i+1)*4*pi*sind(phi);
                        obj.db((panel_to_beam_element(j,i)-1)*6+5,j)=obj.db((panel_to_beam_element(j,i)-1)*6+5,j)+dind(2);%-1/4*panel_to_beam_element(j,i+1)*4*pi*cosd(phi);
                        obj.db((panel_to_beam_element(j,i)-1)*6+6,j)=obj.db((panel_to_beam_element(j,i)-1)*6+6,j)+dind(3);%-1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*sin(dihed);
                        
                        obj.db((panel_to_beam_element(j,i))*6+1,j)=0;
                        obj.db((panel_to_beam_element(j,i))*6+2,j)=0;
                        obj.db((panel_to_beam_element(j,i))*6+3,j)=0;
                        %  obj.db((panel_to_beam_element(j,i))*6+4,j)=obj.db((panel_to_beam_element(j,i))*6+4,j)-1/7*panel_to_beam_element(j,i+1)*4*pi*sind(phi)*cosd(phi);
                        %  obj.db((panel_to_beam_element(j,i))*6+5,j)=obj.db((panel_to_beam_element(j,i))*6+5,j)+1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*cosd(phi)*cosd(dihed);
                        obj.db((panel_to_beam_element(j,i))*6+4,j)=obj.db((panel_to_beam_element(j,i))*6+4,j)+dind(1);%-1/2*panel_to_beam_element(j,i+1)*4*pi*sind(phi);
                        obj.db((panel_to_beam_element(j,i))*6+5,j)=obj.db((panel_to_beam_element(j,i))*6+5,j)+dind(2);%-1/4*panel_to_beam_element(j,i+1)*4*pi*cosd(phi);
                        obj.db((panel_to_beam_element(j,i))*6+6,j)=obj.db((panel_to_beam_element(j,i))*6+6,j)+dind(3);%-1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*sind(dihed);
                    end
                end
            end
        end
        
        function panel_to_beam_element=initialize_beam_to_panel_for_block(obj,beam_idx,k,block,panels,grid,panel_to_beam_element,mirr,comp_idx)
            observer=1;
            plotBtp=0; % debug plot by simon 
            if mirr==1
                grid_offset=length(obj.panels)/2;
            else
                grid_offset=0;
            end
           
            for i=1:obj.wing_segments(k).n_span
                
                n_pan=sum([obj.wing_segments(k).n_le_panels obj.wing_segments(k).n_chord  obj.wing_segments(k).n_te_panels]);

                for j=1:n_pan
                    
                    idx=obj.wing_segments(k).panel_start_idx+n_pan*(i-1)-1+j+grid_offset;
                    % check if current panel lies within block
                    % if ymin of panel is smaller or equal to ymax of block  OR  ymax of panel is smaller or equal to min 
                    if ~(min(grid(comp_idx,panels(:,idx)))> 100*eps+max(block(comp_idx,:))) || ~(max(grid(comp_idx,panels(:,idx)))> min(block(comp_idx,:)-100*eps))
                        
                        if size(block,2)==3
                            in=points_in_triangle(block,grid(:,panels(:,idx)));
                        else
                            in=points_in_quad(block,grid(:,panels(:,idx))); %which points of current panel lie in block
                        end
                        
                        if in(1)==1 && in(2)==1 && in(3)==1 && in(4)==1     %if panel lies completely in block
                            if plotBtp==1
                                observer=0;
                                figure;
                                h=fill3(block(1,:),block(2,:),block(3,:),'b');
                                set(h,'facealpha',.25);
                                hold on;
                                axis equal;
                                plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],1),obj.wingbox_coords(2,[beam_idx, beam_idx+1],1),obj.wingbox_coords(3,[beam_idx, beam_idx+1],1),'--ks')
                                plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],2),obj.wingbox_coords(2,[beam_idx, beam_idx+1],2),obj.wingbox_coords(3,[beam_idx, beam_idx+1],2),'--ks')
                                plot3((obj.wingbox_coords(1,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(1,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(2,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(2,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(3,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(3,[beam_idx, beam_idx+1],2))/2,'-ks')
                                h=fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'g');
                                set(h,'facealpha',.25);
                            end
                            
                            
                            rb=0.25*obj.wingbox_coords(:,beam_idx+1,1)+0.25*obj.wingbox_coords(:,beam_idx+1,2)+0.25*obj.wingbox_coords(:,beam_idx,1)+0.25*obj.wingbox_coords(:,beam_idx,2);
                            r=0.5*(0.75*grid(:,panels(1,idx))+0.25*grid(:,panels(4,idx)))+0.5*(0.75*grid(:,panels(2,idx))+0.25*grid(:,panels(3,idx)));
                            vecb=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2))-(0.5*obj.wingbox_coords(:,beam_idx,1)+0.5*obj.wingbox_coords(:,beam_idx,2));
                            theta=cross(vecb,r-rb)/(norm(r-rb)*norm(vecb));
                            s=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)) + ((vecb / norm(vecb)) * (r - (0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)))') * (vecb / norm(vecb));
%                                     if theta(3)<0
%                                         dist=-norm(r-rb);
%                                     else
%                                         dist=norm(r-rb);
%                                     end
                            dist=r-s;
                            if observer==1
                                fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),[0.2 0.7 0.4])
                            end
                            panel_to_beam_element(idx,1)=beam_idx;
                            panel_to_beam_element(idx,2)=1.000;
                            panel_to_beam_element(idx,3:5)=dist;
                        else   %panel needs to be cut and only portion within block goes to panel to beam element matrix
                                [newPoints,faceNewPoint]=cut_polygons(block,grid(:,panels(:,idx)));
                            if sum(faceNewPoint)>=1 
                                
                                %plot3(new_polygon(1,:),new_polygon(2,:),new_polygon(3,:),'bo');
                                if size(block,2)==3
                                    idOrigPanelPointsInBlock=points_in_triangle_poly(block,grid(:,panels(:,idx)));
                                else
                                    idOrigPanelPointsInBlock=points_in_quad_poly(block,grid(:,panels(:,idx))); 
                                end

                                %red_polygon is the polygon inside the
                                %block consisting of newPoints and
                                %origPanelPointsInBlock
                                nRedPolygonPoints=size(newPoints,2)+sum(idOrigPanelPointsInBlock);
                                red_polygon=zeros(3,nRedPolygonPoints);
                                iRedPolygonPoints=1;
                                for iFace=1:4  %for every face of orig Panel
                                    if idOrigPanelPointsInBlock(iFace)==1
                                        red_polygon(:,iRedPolygonPoints)=grid(:,panels(iFace,idx));
                                        iRedPolygonPoints=iRedPolygonPoints+1;
                                    end
                                    for jNewPoints=1:length(faceNewPoint)
                                        if iFace==faceNewPoint(jNewPoints)
                                            red_polygon(:,iRedPolygonPoints)=newPoints(:,jNewPoints);
                                            iRedPolygonPoints=iRedPolygonPoints+1;
                                        end
                                    end
                                end
                                %fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),0.1)
                                if size(red_polygon,1)>=3 && size(red_polygon,2)>=3
                                    if observer==1
                                        handle=fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),[0.2 0.7 0.4]);
                                        alpha(handle,0.4);
                                    end
                                    if plotBtp==1
                                        observer=0;
                                        figure;
                                        h=fill3(block(1,:),block(2,:),block(3,:),'b');
                                        set(h,'facealpha',.25);
                                        hold on;
                                        axis equal;
                                        plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],1),obj.wingbox_coords(2,[beam_idx, beam_idx+1],1),obj.wingbox_coords(3,[beam_idx, beam_idx+1],1),'--ks')
                                        plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],2),obj.wingbox_coords(2,[beam_idx, beam_idx+1],2),obj.wingbox_coords(3,[beam_idx, beam_idx+1],2),'--ks')
                                        plot3((obj.wingbox_coords(1,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(1,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(2,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(2,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(3,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(3,[beam_idx, beam_idx+1],2))/2,'-ks')
                                        handle=fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),'y');
                                        scatter3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'dr')
                                        alpha(handle,0.25);
                                    end
                                        
                                    panel_area=1/2*norm(cross(grid(:,panels(3,idx))-grid(:,panels(1,idx)),grid(:,panels(4,idx))-grid(:,panels(2,idx))));
                                    moment_line=[(0.75*grid(:,panels(1,idx))+0.25*grid(:,panels(4,idx))) (0.75*grid(:,panels(2,idx))+0.25*grid(:,panels(3,idx)))];
                                    red_moment_line=cut_line_polygon(red_polygon,moment_line);
                                    %plot3(red_moment_line(1,:),red_moment_line(2,:),red_moment_line(3,:),'-ko');
                                
                                    if isnan(panel_area)
                                        panel_area;
                                    end
                                    if isnan(abs(planar_polygon_area(red_polygon)))
                                        red_polygon;
                                        abs(planar_polygon_area(red_polygon));
                                    end
                                    
                                    rb=0.25*obj.wingbox_coords(:,beam_idx+1,1)+0.25*obj.wingbox_coords(:,beam_idx+1,2)+0.25*obj.wingbox_coords(:,beam_idx,1)+0.25*obj.wingbox_coords(:,beam_idx,2);
                                    %r= sum(red_polygon,2)/size(red_polygon,2);
                                    r=0.5*red_moment_line(:,1)+0.5*red_moment_line(:,end); 
                                    vecb=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2))-(0.5*obj.wingbox_coords(:,beam_idx,1)+0.5*obj.wingbox_coords(:,beam_idx,2));
                                    theta=cross(vecb,r-rb)/(norm(r-rb)*norm(vecb));
                                    s=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)) + ((vecb / norm(vecb)) * (r - (0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)))') * (vecb / norm(vecb));
%                                     if theta(3)<0
%                                         dist=-norm(r-rb);
%                                     else
%                                         dist=norm(r-rb);
%                                     end
                                    dist=r-s;
                                    
                                    if plotBtp==1
                                        plot3(moment_line(1,:),moment_line(2,:),moment_line(3,:),'-k');
                                        plot3(red_moment_line(1,:),red_moment_line(2,:),red_moment_line(3,:),'-ro');
                                        scatter3(rb(1),rb(2),rb(3),'r','filled');
                                        scatter3(r(1),r(2),r(3),'rd','filled');
                                        scatter3(s(1),s(2),s(3),'bd','filled');
                                        line([r(1), s(1)],[r(2), s(2)],[r(3), s(3)])
                                    end
                                    if panel_to_beam_element(idx,2)==0
                                        panel_to_beam_element(idx,1)=beam_idx;
                                        if abs(planar_polygon_area(red_polygon))/panel_area<1
                                            panel_to_beam_element(idx,2)=abs(planar_polygon_area(red_polygon))/panel_area;
                                        else
                                            panel_to_beam_element(idx,2)=1;
                                        end
                                        if isnan(panel_to_beam_element(idx,2))
                                            panel_to_beam_element(idx,2)=0;
                                        end
                                        panel_to_beam_element(idx,3:5)=dist;
                                        
                                    elseif panel_to_beam_element(idx,7)==0
                                        panel_to_beam_element(idx,6)=beam_idx;
                                        if abs(planar_polygon_area(red_polygon))/panel_area<1
                                            panel_to_beam_element(idx,7)=abs(planar_polygon_area(red_polygon))/panel_area;
                                        else
                                           panel_to_beam_element(idx,7)=1; 
                                        end
                                        panel_to_beam_element(idx,8:10)=dist;
                                        
                                    else
                                        panel_to_beam_element(idx,12);
                                        
                                        panel_to_beam_element(idx,11)=beam_idx;
                                        share=abs(planar_polygon_area(red_polygon));
                                        if ~(isinf(share)) && ~(isnan(share))
                                            panel_to_beam_element(idx,12)=share/panel_area;
                                            panel_to_beam_element(idx,13:15)=dist;
                                        else
                                            panel_to_beam_element(idx,12)=0;
                                            panel_to_beam_element(idx,7)=0;
                                        end
                                    end
                                else
                                    if panel_to_beam_element(idx,1)==0
                                        panel_to_beam_element(idx,1)=beam_idx;
                                        panel_to_beam_element(idx,2)=0.000;
                                        panel_to_beam_element(idx,3:5)=0.000;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        function obj=compute_beam_forces(obj,panel_to_beam_element,panel_forces,panel_moments, beam)
            if obj.symmetric==1
                f=2;
            else
                f=1;
            end
            k=obj.panel_start_idx;
            obj.beam_forces_structmesh=zeros(3,size(obj.wingbox_coords,2)*f-1);
            obj.beam_moments_structmesh=zeros(3,size(obj.wingbox_coords,2)*f-1);
            end_idx=(obj.panel_start_idx+size(obj.panels,2)-1);
            ctr=1;
            overall_force=0;
            overall_force_loc=0;
            for j=obj.panel_start_idx:1:end_idx                             %for every panel i
                for i=1:5:size(panel_to_beam_element,2)-1                   %for all 3 corresponding beamelements j of the panel 
                    if not(panel_to_beam_element(j,i)==0)                   %only if there is parts of the panels assigned to the current beamelement j
                        beamlength=beam.beamelement(panel_to_beam_element(k,i)).le;%norm(obj.wingbox_coords(:,panel_to_beam_element(k,i)+1)-obj.wingbox_coords(:,panel_to_beam_element(k,i)));
                        %Trot=beam.beamelement(panel_to_beam_element(j,i)).f_rotM_6dof(real(beam.beamelement(panel_to_beam_element(j,i)).nu),real(beam.beamelement(panel_to_beam_element(j,i)).epsilon),0);
                        if obj.symmetric==1
                            tf=panel_forces(:,j);
                            %tf(2)=-tf(2);                                   %mirror global y component
                            force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*tf;%[0 0 dot(crp,tf)]';
                            moment_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_moments(:,j);
                            %beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)
                            force_glob=tf;
                            % force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_forces(:,j);
                            % force_loc(2)=-force_loc(2);
                        else
                            
                           tf=panel_forces(:,j);
                          %  tf(2)=tf(2);
                           force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*tf;%[0 0 dot(crp,tf)]';
                            moment_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_moments(:,j);
                            %beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)
                            force_glob=tf;
                            % force_loc=panel_forces(:,j);%Trot(1:3,1:3)*
                            %  force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_forces(:,j);
                          %   force_loc=panel_forces(:,j);%[0 0 dot(crp,panel_forces(:,j))]';
%                             
%                             
                         %    force_glob=panel_forces(:,j);
                        end
                        overall_force=overall_force+force_glob(3)*panel_to_beam_element(j,i+1);
                        overall_force_loc=overall_force_loc+force_loc(3)*panel_to_beam_element(j,i+1);
                        obj.beam_forces_structmesh(:,panel_to_beam_element(j,i))=obj.beam_forces_structmesh(:,panel_to_beam_element(j,i))+force_loc*panel_to_beam_element(j,i+1)/beamlength;
                        %obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))=obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))+sign(force_loc(3))*norm(force_loc(1:2:3))*panel_to_beam_element(j,i+1)*norm(panel_to_beam_element(j,i+2:i+4))/beamlength;
                        % transform distance to local coordinate system
                        dist = beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_to_beam_element(j,i+2:i+4)';
                        obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))=obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))+panel_to_beam_element(j,i+1)*(-force_loc(3)*dist(1)+force_loc(1)*dist(3)+moment_loc(2))/beamlength;
                    end
                end
                ctr=ctr+1;
                k=k+1;
                if k==(size(obj.panels,2)/2)+1
                   k=obj.panel_start_idx; 
                end
            end
%             overall_force
%             overall_force_loc
%             if obj.symmetric==1
%                 mm=obj.beam_forces_structmesh;
%                 obj.beam_forces_structmesh=mm(:,end-1:-1:1);
%                 mm=obj.beam_moments_structmesh;
%                 obj.beam_moments_structmesh=mm(:,end-1:-1:1);   
%             end          
        end
        
        function obj=compute_c4_forces(obj,panel_forces,panels,grid)
            
            obj.c4_forces=[];
            obj.c4_moments=[];
            
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_forces(panel_forces,panels,grid);
                obj.c4_forces=[obj.c4_forces obj.wing_segments(i).span_forces];
                obj.c4_moments=[obj.c4_moments obj.wing_segments(i).span_moments];
            end
            
            obj.c4_forces_structmesh=zeros(3,length(obj.wingbox_rl_coords_mid));
            obj.c4_moments_structmesh=zeros(3,length(obj.wingbox_rl_coords_mid));
            obj.c4_forces_structmesh(1,:)=interp1(obj.c4_rl_coords,obj.c4_forces(1,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_forces_structmesh(2,:)=interp1(obj.c4_rl_coords,obj.c4_forces(2,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_forces_structmesh(3,:)=interp1(obj.c4_rl_coords,obj.c4_forces(3,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(1,:)=interp1(obj.c4_rl_coords,obj.c4_moments(1,:),obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(2,:)=interp1(obj.c4_rl_coords,obj.c4_moments(2,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(3,:)=interp1(obj.c4_rl_coords,obj.c4_moments(3,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            
            %             figure
            %             hold on
            %             plot(c4_rl_coord,obj.c4_forces(3,:));
            %             plot(c4_rl_coord_sm,obj.c4_forces_structmesh(3,:));
        end
        
        function grid=compute_deflected_grid(obj,panels,grid,deflections_structmesh)
            % computation of  station of panels along the c4 line 
            c4_rl_coords_edge=0;
            for i=1:length(obj.wing_segments)
                x=0:1/obj.wing_segments(i).n_span:1;
                %seg_rl : stations times lenght of the c4 line of the segment
                seg_rl=obj.wing_segments(i).b*x/cosd(obj.wing_segments(i).c4_sweep);
                c4_rl_coords_edge=[c4_rl_coords_edge seg_rl(2:end)+c4_rl_coords_edge(end)];
            end
            
            for j=1:6
                def(j,:)=squeeze(deflections_structmesh(j:6:end));
                %leftdef(j,:)=interp1(squeeze(obj.mesh_struct(i).y),squeeze(obj.deflect(i).def(j:6:end)),obj.AC_geo.LiftingSurfaceList(i).Panels.s_R,'spline','extrap');
            end
            %             figure
            %             plot(obj.wingbox_rl_coords)
            %             figure
            %             plot(c4_rl_coords_edge)
            %deflections_aeromesh are the deflections at the c/4 line. for
            %determination, the deflections are interpolated onto the
            %aero stations on the c/4 line
            if obj.symmetric==1
                n_def_hw=(length(def(1,:))-1)/2; %halfwing number deflections structmesh
                %left wing
                deflections_aeromesh_left(1,:)=interp1(obj.wingbox_rl_coords, def(1,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(2,:)=interp1(obj.wingbox_rl_coords, def(2,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(3,:)=interp1(obj.wingbox_rl_coords, def(3,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(4,:)=interp1(obj.wingbox_rl_coords, def(4,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(5,:)=interp1(obj.wingbox_rl_coords,def(5,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(6,:)=interp1(obj.wingbox_rl_coords, def(6,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                %right wing
                deflections_aeromesh(1,:)=interp1(obj.wingbox_rl_coords,[ def(1,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(2,:)=interp1(obj.wingbox_rl_coords,[ def(2,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(3,:)=interp1(obj.wingbox_rl_coords,[ def(3,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(4,:)=interp1(obj.wingbox_rl_coords,[ def(4,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(5,:)=interp1(obj.wingbox_rl_coords,[ def(5,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(6,:)=interp1(obj.wingbox_rl_coords,[ def(6,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');

            else
                deflections_aeromesh(1,:)=interp1(obj.wingbox_rl_coords,def(1,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(2,:)=interp1(obj.wingbox_rl_coords,def(2,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(3,:)=interp1(obj.wingbox_rl_coords,def(3,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(4,:)=interp1(obj.wingbox_rl_coords,def(4,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(5,:)=interp1(obj.wingbox_rl_coords,def(5,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(6,:)=interp1(obj.wingbox_rl_coords,def(6,:), c4_rl_coords_edge,'lin','extrap');
            end
            
            start_idx=1;
%                         hold on
%                         plot3(obj.wingbox_coords(1,:,1),obj.wingbox_coords(2,:,1),obj.wingbox_coords(3,:,1),'cx');
%             plot3(obj.wingbox_coords(1,:,1)+def(1,24:end),obj.wingbox_coords(2,:,1)+def(2,24:end),obj.wingbox_coords(3,:,1)+def(3,24:end),'ro');
% %             
            for i=1:length(obj.wing_segments)
                end_idx=start_idx+size(obj.wing_segments(i).c4_coords,2);
                grid=obj.wing_segments(i).compute_deflected_grid(panels,grid,deflections_aeromesh(:,start_idx:end_idx));
                start_idx=end_idx;
            end
            
            
%             if obj.grid_3D==1
%                 for i=1:length(obj.wing_segments)
%                         obj.grid_vol
%                 end  
%             end
            
%                         hold on
%                         for i=1:length(obj.panels)
%                             handle= fill3(grid(1,obj.panels(:,i)), grid(2,obj.panels(:,i)),grid(3,obj.panels(:,i)),'b');
%                             alpha(handle,0.4);
%                         end
%                         axis equal
            start_idx=1;
            if obj.symmetric==1
                for i=1:length(obj.wing_segments)
                    end_idx=start_idx+size(obj.wing_segments(i).c4_coords,2);
                    grid=obj.wing_segments(i).compute_deflected_grid_left(panels,grid,deflections_aeromesh_left(:,start_idx:end_idx),length(obj.panels)/2);
                    start_idx=end_idx;
                end
               % grid(:,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2))=[grid(:,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2) [grid(1,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2);-grid(2,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2);grid(3,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2)]];
            end
          %  obj.grid_deflected=grid;
%             hold on
%                      plot3(obj.wingbox_coords(1,:,1),obj.wingbox_coords(2,:,1),obj.wingbox_coords(3,:,1),'x-');
%                      plot3(obj.wingbox_coords(1,:,2),obj.wingbox_coords(2,:,2),obj.wingbox_coords(3,:,2),'x-');
%                      plot3(0.5*obj.wingbox_coords(1,:,1)+0.5*obj.wingbox_coords(1,:,2),0.5*obj.wingbox_coords(2,:,1)+0.5*obj.wingbox_coords(2,:,2),0.5*obj.wingbox_coords(3,:,1)+0.5*obj.wingbox_coords(3,:,2),'o-');
%                      figure
%                      plot(c4_rl_coords_edge,deflections_aeromesh(1:6,:))
        end
        
        function obj=compute_wingbox_coords(obj,le_max,varargin)
            frontspar=0.25*ones(1,length(obj.wing_segments)+1);
            rearspar=0.75*ones(1,length(obj.wing_segments)+1);
            
            if nargin==4
                frontspar=varargin{1};
                rearspar=varargin{2};
            end
            
            wingbox_coords=[];
            wingbox_c4=[];
            wingbox_height=[];
            
            if(length(frontspar)==2*length(obj.wing_segments))
                incr=2;
            else
                incr=1;
            end
            ctr=1;
            iNodes = 1;
           % incr
            for i=1:length(obj.wing_segments)
      %          i
                n=ceil(obj.wing_segments(i).b/le_max);
                if obj.wing_segments(i).nBeamelements>0
                    if obj.wing_segments(i).nBeamelements<n
                        disp('WARNING: The number of beamelements selected for this wingsegment results in a lower resolution than defined in the grid spacing!')
                    end
                    n = obj.wing_segments(i).nBeamelements;
                else
                    if n<3
                        n=3;
                    end
                end
                %ctr;
%                 frontspar
%                 rearspar
                if obj.isExternalFEM==0
                    obj.wing_segments(i)=obj.wing_segments(i).compute_wingbox_coords(frontspar(ctr:ctr+1),rearspar(ctr:ctr+1),n);
                elseif obj.isExternalFEM==1
                    obj.wing_segments(i) = obj.wing_segments(i).read_wingbox_coords(obj.pathNodeCoords, n, iNodes);
                else
                    disp('ERROR: wing.isExternalFEM not correctly defined!')
                end
                if isempty(wingbox_coords)
                    wingbox_coords=[wingbox_coords obj.wing_segments(i).wingbox_coords];
                    wingbox_c4=[wingbox_c4 obj.wing_segments(i).wingbox_c4];
                    wingbox_height=[wingbox_height obj.wing_segments(i).wingbox_height'];
                else
                    wingbox_coords=[wingbox_coords obj.wing_segments(i).wingbox_coords(:,2:end,:)];
                    wingbox_c4=[wingbox_c4 obj.wing_segments(i).wingbox_c4(:,2:end)];
                    wingbox_height=[wingbox_height obj.wing_segments(i).wingbox_height(2:end,:)'];
                end
                ctr=ctr+incr;
                iNodes = iNodes+n;
            end
            obj.wingbox_coords=wingbox_coords;
            obj.wingbox_c4=wingbox_c4;
            obj.wingbox_height=wingbox_height';
            
            c4_rl_coord_sm(1)=0;
            for i=1:size(obj.wingbox_c4,2)-1
                c4_rl_coord_sm(i+1)=norm(obj.wingbox_c4(1:3,i+1)-obj.wingbox_c4(1:3,i));
            end
            
            c4_rl_coord_sm=cumsum(c4_rl_coord_sm);
            
            obj.wingbox_rl_coords=c4_rl_coord_sm; 
            
            for i=1:length(obj.wingbox_rl_coords)-1
                c4_rl_coord_sm_mid(i)=0.5* obj.wingbox_rl_coords(i)+0.5*(obj.wingbox_rl_coords(i+1));
            end
            

            obj.wingbox_rl_coords_mid=c4_rl_coord_sm_mid;
            
            obj.c4_coords=[];
            for i=1:length(obj.wing_segments)
                obj.c4_coords=[obj.c4_coords obj.wing_segments(i).c4_coords];
            end
            
            c4_rl_coord(1)=norm(obj.c4_coords(1:3,1))-norm(obj.wingbox_c4(1:3,1));
            for i=1:size(obj.c4_coords,2)-1
                c4_rl_coord(i+1)= norm(obj.c4_coords(1:3,i+1)-obj.c4_coords(1:3,i));
            end
            c4_rl_coord=cumsum(c4_rl_coord);
            obj.c4_rl_coords=c4_rl_coord;
            
            %             plot3(wingbox_coords(1,:,1),wingbox_coords(2,:,1),wingbox_coords(3,:,1),'x-');
            %             plot3(wingbox_coords(1,:,2),wingbox_coords(2,:,2),wingbox_coords(3,:,2),'x-');
            
        end
        
        function obj=plot_grid(obj)
            hold on
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                alpha(handle,0.4);
            end
        end
        
        function obj=plot_grid_vol(obj)
            hold on
            for i=1:length(obj.panels)
                if obj.is_te(i)==0
                    handle= fill3(obj.grid_vol_up(1,obj.panels(:,i)), obj.grid_vol_up(2,obj.panels(:,i)),obj.grid_vol_up(3,obj.panels(:,i)),'b');
                else
                    handle= fill3(obj.grid_vol_up(1,obj.panels(:,i)), obj.grid_vol_up(2,obj.panels(:,i)),obj.grid_vol_up(3,obj.panels(:,i)),'r');   
                end
                alpha(handle,1);
            end
            
            for i=1:length(obj.panels)
                 if obj.is_te(i)==0
                    handle= fill3(obj.grid_vol_lo(1,obj.panels(:,i)), obj.grid_vol_lo(2,obj.panels(:,i)),obj.grid_vol_lo(3,obj.panels(:,i)),'b');
                 else
                     handle= fill3(obj.grid_vol_lo(1,obj.panels(:,i)), obj.grid_vol_lo(2,obj.panels(:,i)),obj.grid_vol_lo(3,obj.panels(:,i)),'r');
                 end 
                alpha(handle,1);
            end    
        end
        
        
        function obj=plot_grid_flat(obj)
            hold on
            %fprintf(fileID,'VARIABLES = "X", "Y", "Z"\n');
            for i=1:length(obj.panels)
                handle= fill3(obj.grid_flat(1,obj.panels(:,i)), obj.grid_flat(2,obj.panels(:,i)),obj.grid_flat(3,obj.panels(:,i)),'r');
                alpha(handle,0.4);
            end
        end
        
        function obj=write_tecplot_grid(obj,fileID)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).write_tecplot_grid(fileID,obj.name,i);
            end
        end
        
        function obj=read_xml_definition(obj,xmlstruct)
            di=0;
            if strcmp(xmlstruct.tag,'WING')
                obj.name=xmlstruct.attribs(1).value;
                obj.symmetric=str2double(xmlstruct.attribs(2).value);
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
                            n=length(obj.wing_segments);
                            if n~=0
                                obj.wing_segments(n+1)=class_wingsegment(xmlstruct.child(i));
                                obj.wing_segments(n+1).symmetric=obj.symmetric;
                            else
                                obj.wing_segments=class_wingsegment(xmlstruct.child(i));
                                obj.wing_segments.symmetric=obj.symmetric;
                            end
                        else
                            if  strcmp(xmlstruct.child(i).child(1).child(1).tag,'ATTACH_TO')
                                if strcmp(xmlstruct.child(i).child(1).child(1).value,'previous')
                                    attach_pos=obj.wing_segments(end).get_tip_c4point();
                                    xmlstruct.child(i).child(1).child(1).tag='X';
                                    xmlstruct.child(i).child(1).child(2).tag='Y';
                                    xmlstruct.child(i).child(1).child(3).tag='Z';
                                    xmlstruct.child(i).child(1).child(1).value=num2str(attach_pos(1));
                                    xmlstruct.child(i).child(1).child(2).value=num2str(attach_pos(2));
                                    xmlstruct.child(i).child(1).child(3).value=num2str(attach_pos(3));
                                    n=length(obj.wing_segments);
                                    obj.wing_segments(n+1)=class_wingsegment(xmlstruct.child(i));
                                    obj.wing_segments(n+1).symmetric=obj.symmetric;
                                end
                            end
                        end
                        
                        try
                            cs_idx=11;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'CONTROL_SURFACE')
                                if strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'trailing_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                    end
                                    
                                    if (str2double(xmlstruct.child(i).child(cs_idx).child(2).value)==0) && (str2double(xmlstruct.child(i).child(cs_idx).child(3).value)==1)
                                        obj.wing_segments(i+di)=obj.wing_segments(i+di).add_control_surface(te_device);
                                    else
                                        n=length(obj.wing_segments);
                                        obj=obj.split_segment(n,str2double(xmlstruct.child(i).child(cs_idx).child(2).value),str2double(xmlstruct.child(i).child(cs_idx).child(3).value),te_device);
                                        if(str2double(xmlstruct.child(i).child(cs_idx).child(2).value)~=0)
                                            di=di+1;
                                        end
                                        if(str2double(xmlstruct.child(i).child(cs_idx).child(3).value)~=0)
                                            di=di+1;
                                        end
                                    end
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'leading_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        le_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl);
                                    catch
                                        le_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value);
                                    end
                                    
                                    obj.wing_segments(i+di)=obj.wing_segments(i+di).add_control_surface(le_device);
                                end
                            % NBEAMELEMENTS defines the number of the wingsegment's beamelements. This number (if defined) will be used instead of the grid spacing
                            elseif strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end
                        end
                        
                        try
                            cs_idx=12;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'CONTROL_SURFACE')
                                if strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'trailing_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                    end
                                    %te_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value);
                                    obj.wing_segments(i+di)=obj.wing_segments(i+di).add_control_surface(te_device);
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'leading_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        le_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl);
                                    catch
                                        le_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value);
                                    end
                                    %le_device=class_control_surface(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value);
                                    obj.wing_segments(i+di)=obj.wing_segments(i+di).add_control_surface(le_device);
                                end
                            elseif strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end  
                        end
                        
                        try
                            cs_idx=13;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
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

