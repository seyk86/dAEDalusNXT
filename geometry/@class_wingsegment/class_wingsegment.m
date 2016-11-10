%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_wingsegment
    %class wingsegment geometric definitions for a wing segment
    %   Detailed explanation goes here
    
    properties
        %> reference position (c4 position)
        pos;
        %> is this surface symmetric?
        symmetric;
        %> grid start index of this segment in the global grid
        grid_start_idx;
        %> panel start index of this segment in the global grid
        panel_start_idx;
        %> quarter chordline coordinates
        c4_coords;
        %> spanwise forces
        span_forces;
        %> spanwise moments
        span_moments;
        %> special case
        real_sweep=[];
        %> quarter chordline sweep
        c4_sweep;
        %> leading edge sweep
        le_sweep;
        %> taper ratio
        TR;
        %> span
        b;
        %> root chord
        c_r;
        %> tip chord
        c_t;
        %> root profile data points
        profile_r;
        %> name of root profile
        profile_name_r;
        %> tip profile data points
        profile_t;
        %> name of root profile
        profile_name_t;
        %> root twist
        Theta_r;
        %> tip twist
        Theta_t;
        %> dihedral at c4 line
        dihed;
        %> dihedral at leading edge line
        le_dihed;
        
        %> thickness ratio
        tc=0.2;
        %> panel edge coordinates
        xyz;
        %> part of wing fixed
        xyz_fixed;
        %> edge coordinates of trailing edge device
        xyz_te_device;
        %> edge coordinates of leading edge device
        xyz_le_device;
        %> segment area
        S;
        %> segment wetted area
        S_wet;
        %> mean aerodynamic chord
        c_mac;
        
        is_laminar;
        
        skeleton_line_r;
        skeleton_line_t;
        
        has_te_cs=0; % information flag if there is a trailing edge control surface
        has_le_cs=0; % information flag if there is a leading edge control surface
        
        le_device; % class containing the information about the leading edge control surface
        te_device; % class containing the information about the leading edge control surface
        
        c_le_device;
        c_te_device;
        
        delta_te_device;
        delta_le_device;
        
        te_max_defl;
        te_min_defl;
        
        n_span;
        n_chord;
        n_le_panels;
        n_te_panels;
        
        D_f;
        CD_f;
        
        grid;
        grid_flat;
        %for potential flow
        is_te;
        grid_vol_upper;
        grid_vol_lower;
        % grid for wake
        grid_wake;
        % panels for wake
        panels_wake;
        
        
        
        panels;
        nxt_ds_pt;
        te_idx;
        
        grid_le;
        panels_le;
        grid_te;
        panels_te;
        
        wingbox_coords;
        wingbox_c4;
        wingbox_height;
        
        %> number of beam elements; will be determined by grid spacing if
        % zero
        nBeamelements = 0;
    end
    
    methods
        
        function obj = class_wingsegment(varargin)
            
            % input parameter is edge coordinates
            if nargin==1
                obj=obj.read_xml_definition(varargin{1});
            elseif nargin==4
                obj.profile_name_r='flat_plate';
                obj.profile_name_t='flat_plate';
                obj=obj.load_airfoils();
                sz=size(varargin{1});
                if sz(1)==3 && sz(2)==4
                    obj.xyz=varargin{1};
                end
                
                obj=obj.compute_parameters_from_coords();
                obj.symmetric=varargin{2};
                obj.profile_r=varargin{3};
                obj.profile_t=varargin{4};
                % input parameter is edge coordinates and control surfaces
            elseif nargin==6
                obj.profile_name_r='flat_plate';
                obj.profile_name_t='flat_plate';
                obj=obj.load_airfoils();
                sz=size(varargin{1});
                if sz(1)==3 && sz(2)==4
                    obj.xyz=varargin{1};
                end
                
                obj=obj.compute_parameters_from_coords();
                obj.symmetric=varargin{2};
                obj.profile_r=varargin{5};
                obj.profile_t=varargin{6};
                
                if strcmp(varargin{3},'cs_te')
                    cs_te=varargin{4};
                    for i=1:length(cs_te)
                        if i<length(cs_te)
                            obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                        else
                            obj.c_te_device(i)=obj.c_r*cs_te(i);
                        end
                        obj.delta_te_device(i)=25*i;
                    end
                    obj.has_te_cs=1;
                elseif strcmp(varargin{3},'cs_le')
                    cs_le=varargin{4};
                    for i=1:length(cs_le)
                        if i<length(cs_le)
                            obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i+1));
                        else
                            obj.c_le_device(i)=obj.c_r*cs_le(i);
                        end
                        obj.delta_le_device(i)=-10*i;
                    end
                    obj.has_le_cs=1;
                end
                
                % standard initialization
            else
                obj.pos=varargin{1};
                obj.symmetric=varargin{2};
                obj.dihed=varargin{3};
                
                LambdaSpec=varargin{4};
                if strcmp(LambdaSpec,'c4-sweep')
                    obj.c4_sweep=varargin{5};
                elseif  strcmp(LambdaSpec,'real-sweep')  
                    obj.real_sweep=varargin{5};
                elseif strcmp(LambdaSpec,'le-sweep')
                    obj.le_sweep=varargin{5};
                end
                
                if nargin==5
                    obj.profile_name_r='flat_plate';
                    obj.profile_name_t='flat_plate';
                    obj=obj.load_airfoils();
                    obj.Theta_r=0;
                    obj.Theta_t=0;
                    obj.b=5;
                    obj.TR=0.5;
                    obj.c_r=3;
                    obj.S=obj.b/2*(obj.c_r*(1+obj.TR));
                    obj.S_wet=2*obj.S;
                    obj.c_t=obj.TR*obj.c_r;
                    obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
                    
                    if strcmp(LambdaSpec,'c4-sweep')
                        obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
                    elseif strcmp(LambdaSpec,'le-sweep')
                        obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
                    elseif strcmp(LambdaSpec,'real-sweep')
                        obj.le_sweep=[];
                        obj.c4_sweep=[];
                    end
                elseif nargin==19
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        
                        obj=obj.complete_params_from_stdinit();
                        obj.has_le_cs=0;
                        obj.has_te_cs=0;
                    end
                    
                elseif nargin==21
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    property8=varargin{20};
                    
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))&&(strcmp(property8,'cs_te'))
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs=varargin{21};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs)
                            if i<length(cs)
                                obj.c_te_device(i)=obj.c_r*(cs(i)-cs(i+1));
                                
                            else
                                obj.c_te_device(i)=obj.c_r*cs(i);
                            end
                            obj.delta_te_device(i)=-10*i;
                        end
                        
                        obj.has_te_cs=1;
                    end
                    
                elseif nargin==23
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    property8=varargin{20};
                    property9=varargin{22};
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))&&(strcmp(property8,'cs_te'))&&(strcmp(property9,'cs_le'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs_te=varargin{21};
                        cs_le=varargin{23};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs_te)
                            if i<length(cs_te)
                                obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                            else
                                obj.c_te_device(i)=obj.c_r*cs_te(i);
                            end
                            obj.delta_te_device(i)=-10*i;
                        end
                        
                        for i=1:length(cs_le)
                            if i==1
                                obj.c_le_device(i)=obj.c_r*cs_le(i);
                            else
                                obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i-1));
                            end
                            obj.delta_le_device(i)=0;
                        end
                        obj.has_te_cs=1;
                        obj.has_le_cs=1;
                    elseif(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'cs_te'))&&(strcmp(property7,'delta_te'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs_te=varargin{21};
                        delta_te=varargin{23};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs_te)
                            if i<length(cs_te)
                                obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                                
                            else
                                obj.c_te_device(i)=obj.c_r*cs_te(i);
                            end
                            obj.delta_te_device(i)=0;
                        end
                        
                        for i=1:length(delta_te)
                            obj.delta_te_device(i)=delta_te(i);
                        end
                        obj.has_te_cs=1;
                    end
                end
                obj=obj.compute_segment_coordinates();
            end
            % compute flap edge points
            obj=obj.compute_controlsurface_coordinates();
            obj=obj.compute_xyz_fixed();
        end
        
        function obj=add_control_surface(obj,control_surface)
            switch(control_surface.pos)
                case 0
                    obj.te_device=control_surface;
                    cs_te=control_surface.hinge;
                    
                    for i=1:length(cs_te)
                        if i<length(cs_te)
                            obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                        else
                            obj.c_te_device(i)=obj.c_r*cs_te(i);
                        end
                        obj.delta_te_device(i)=0;
                    end
                    obj.has_te_cs=1;
                    
                case 1
                    obj.le_device=control_surface;
                    cs_le=control_surface.hinge;
                    for i=1:length(cs_le)
                        if i==1
                            obj.c_le_device(i)=obj.c_r*cs_le(i);
                        else
                            obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i-1));
                        end
                        obj.delta_le_device(i)=0;
                    end
                    obj.has_le_cs=1;
            end
            
            obj=obj.compute_controlsurface_coordinates();
            obj=obj.compute_xyz_fixed();
        end
        
        
        function obj= complete_params_from_stdinit(obj)
            %                                     if strcmp(LambdaSpec,'c4-sweep')
            %                             obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
            %                         elseif strcmp(LambdaSpec,'le-sweep')
            %                             obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
            %                         end
            if isempty(obj.le_sweep)
                obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
            elseif isempty(obj.c4_sweep)
                obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
            end
            
            obj.TR=obj.c_t/obj.c_r;
            obj.S=obj.b/2*(obj.c_r*(1+obj.TR));
            obj.S_wet=2*obj.S;
            obj.c_mac =obj.c_r * 2/3 * (( 1 + obj.TR +obj.TR^2 )/ ( 1 + obj.TR ));
        end
        
        function obj = compute_friction_drag(obj,state,S_ref,varargin)
            obj.S_wet=2*obj.S;
            Re_l=state.rho_air*norm(state.V_inf)*obj.c_mac/state.mu;
            if nargin==3
                % laminar or turbulent
                if Re_l>5E5
                    cf=0.455/log10(Re_l)^2.58;
                else
                    cf=1.328/sqrt(Re_l);
                end
                
                if obj.is_laminar==1
                    cf=1.328/sqrt(Re_l);
                end
            end
            % controlsurface gap drag
            % http://adg.stanford.edu/aa241/drag/gapdrag.html
            CD_gaps=0;
            if obj.has_te_cs
                CD_gaps = 0.0002*cos(obj.c4_sweep*pi/180)^2*obj.S/S_ref;
            end
            % form factor formula
            % (adg.stanford.edu/aa241/drag/lsformfactor.html)
            if state.Ma*cos(obj.le_sweep*pi/180)>1
                C=0;
            else
                C=1.1;
            end
            
            cossw2=cos(obj.le_sweep*pi/180)^2;
            Ma2=state.Ma^2;
            form_factor=1+2*C*obj.tc*cossw2/sqrt(1-Ma2*cossw2)+C^2*cossw2*obj.tc^2*(1+5*cossw2)/(2*(1-Ma2*cossw2));
            
            obj.CD_f=form_factor*cf*obj.S_wet/S_ref+CD_gaps;
            obj.D_f=form_factor*1/2*state.rho_air*norm(state.V_inf)^2*cf*obj.S_wet;
        end
        
        function obj=load_airfoils(obj)
            airfoil_name_r=obj.profile_name_r;
            airfoil_name_t=obj.profile_name_t;
            if strcmp(airfoil_name_r,'flat_plate')
                obj.profile_name_r='flat_plate';
                obj.profile_r=[3 3;0 0;0.5 0;1 0;0 0;0.5 0;1 0];
            else
                obj.profile_name_r=['airfoil/',airfoil_name_r,'.DAT'];
                obj.profile_r=load(['airfoil/',airfoil_name_r,'.DAT']);
            end
            
            if strcmp(airfoil_name_r,'flat_plate')
                obj.profile_t=[3 3;0 0;0.5 0;1 0;0 0;0.5 0;1 0];
                obj.profile_name_t='flat_plate';
            else
                obj.profile_t=load(['airfoil/',airfoil_name_t,'.DAT']);
                obj.profile_name_t=['airfoil/',airfoil_name_t,'.DAT'];
            end
        end
        
        function obj=compute_parameters_from_coords(obj)
            dv=obj.xyz(2:3,2)-obj.xyz(2:3,1);
            obj.b=hypot(dv(1),dv(2));
            obj.c_r=sqrt(sum((obj.xyz(:,4)-obj.xyz(:,1)).^2));
            obj.c_t=sqrt(sum((obj.xyz(:,3)-obj.xyz(:,2)).^2));
            obj.S=(obj.c_t+obj.c_r)/2*obj.b;
            
            obj.pos=[obj.xyz(:,1)+0.25*(obj.xyz(:,4)-obj.xyz(:,1))]';
            pos_t=obj.xyz(:,2)+0.25*(obj.xyz(:,3)-obj.xyz(:,2));
            
            dx=pos_t(1)-obj.pos(1);
            dy=pos_t(2)-obj.pos(2);
            dz=pos_t(3)-obj.pos(3);
            obj.c4_sweep=atand(dx/obj.b);
            if dy>1e-9
                obj.dihed=asind(dz/obj.b);
            elseif dy<-1e-9
                obj.dihed=180-asind(dz/obj.b);
            else
                obj.dihed=90;
            end
            
            obj.le_sweep=atand((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b);
            
            dz=obj.xyz(3,4)-obj.xyz(3,1);
            dy=obj.xyz(2,4)-obj.xyz(2,1);
            if obj.dihed==90
                dn=-sqrt(dz^2+dy^2)*sign(dy);
            else
                dn=sqrt(dz^2+dy^2)*sign(dz);
            end
            obj.Theta_r=-asind(dn/obj.c_r);
            
            dz=obj.xyz(3,3)-obj.xyz(3,2);
            dy=obj.xyz(2,3)-obj.xyz(2,2);
            
            if obj.dihed==90
                dn=-sqrt(dz^2+dy^2)*sign(dy);
            else
                dn=sqrt(dz^2+dy^2)*sign(dz);
            end
            
            obj.le_dihed=asind((obj.xyz(3,2)-obj.xyz(3,1))/obj.b);
            
            obj.Theta_t=-asind(dn/obj.c_t);
            obj.TR=obj.c_t/obj.c_r;
            
            %obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
            obj.c_mac =obj.c_r * 2/3 * (( 1 + obj.TR +obj.TR^2 )/ ( 1 + obj.TR ));
            %c_mac= obj.c_r-(2*(obj.c_r-obj.c_t)*(0.5*obj.c_r+obj.c_t)/(3*(obj.c_r+obj.c_t)))
            %obj.c_mac=c_mac;
        end
        
        function obj=set_cs_deflections(obj,varargin)
            if nargin==3
                if strcmp(varargin{1},'te')
                    obj.delta_te_device=varargin{2};
                elseif strcmp(varargin{1},'le')
                    obj.delta_le_device=varargin{2};
                end
            elseif nargin==5
                if strcmp(varargin{1},'te') && strcmp(varargin{3},'le')
                    obj.delta_te_device=varargin{2};
                    obj.delta_le_device=varargin{4};
                elseif strcmp(varargin{3},'le') && strcmp(varargin{1},'te')
                    obj.delta_te_device=varargin{4};
                    obj.delta_le_device=varargin{2};
                end
            else
                
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        function obj=f_deflect_control_surface(obj,name,deflection,varargin)
            side=[];
            if nargin==4
                side=varargin{1};
            end
            
            if strcmp(name,obj.te_device.name)
                if ~isempty(side)
                    if strcmp(side,'left')
                        obj.te_device.delta=deflection;
                        obj.te_device.delta_l_r(1)=deflection;
                    elseif strcmp(side,'right')
                        obj.te_device.delta_l_r(2)=deflection;
                    end
                end
                obj.te_device.delta=deflection;
            elseif strcmp(name,obj.le_device.name)
                obj.le_device.delta=deflection;
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        function obj=deflect_control_surface(obj,name,delta)
            if strcmp(obj.te_device.name,name)
                %% TODO: check if right format!
                obj.delta_te_cs=delta;
            elseif strcmp(obj.le_device.name,name)
                obj.delta_le_cs=delta;
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        
        function obj=compute_segment_coordinates(obj)
            
            if isempty(obj.real_sweep)
                %the following code does not work properly, e.g. input span
                %10m, sweep 50, dihed 20 -> real span = 15m
            % compute segment edge points
            p1=obj.pos+[-obj.c_r/4*cos(obj.Theta_r*pi/180), -obj.c_r/4*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180), +obj.c_r/4*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
            p2=obj.pos+[+obj.b*tan(obj.c4_sweep*pi/180), obj.b*cos(obj.dihed*pi/180), +obj.b*sin(obj.dihed*pi/180)]+[-obj.c_t/4*cos(obj.Theta_t*pi/180) -obj.c_t/4*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +obj.c_t/4*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
            p3=obj.pos+[+obj.b*tan(obj.c4_sweep*pi/180), obj.b*cos(obj.dihed*pi/180), +obj.b*sin(obj.dihed*pi/180)]-[-obj.c_t*3/4*cos(obj.Theta_t*pi/180) -obj.c_t*3/4*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +obj.c_t*3/4*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
            p4=obj.pos-[-obj.c_r*3/4*cos(obj.Theta_r*pi/180), -obj.c_r*3/4*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180), +obj.c_r*3/4*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
            
            obj.le_dihed=asind((p2(3)-p1(3))/obj.b);
            obj.xyz=[p1' p2' p3' p4'];
            else
                %the following code does not work and is commented out
%                 p1=[0 0 0];
%                 p2=[0 obj.b 0];
%                 p3=[obj.c_r obj.b 0];
%                 p4=[obj.c_r 0 0];
%                 a=obj.dihed;
%                 b=obj.Theta_r;
%                 
%                 c=obj.real_sweep;
%                 Lx=[1       0       0
%                     0   cosd(a)  sind(a)
%                     0   -sind(a) cosd(a)];
%                 
%                 Ly=[cosd(b) 0 sind(b)
%                     0      1    0
%                     -sind(b)  0   cosd(b)];
%                 
%                 Lz=[cosd(c) sind(c)   0
%                     -sind(c) cosd(c)  0
%                     0           0   1];
%                 
%                 T=Lz*Ly*Lx;
%                 
%                 p1=obj.pos+(T*p1')';
%                 p2=obj.pos+(T*p2')';
%                 p3=obj.pos+(T*p3')';
%                 p4=obj.pos+(T*p4')';
%                 obj.xyz=[p1' p2' p3' p4'];
                    disp('warning, this does not work properly')
            end
            
        end
        
        function obj=compute_controlsurface_coordinates(obj)                
            if obj.has_te_cs
                obj.delta_te_device=obj.te_device.delta;
            elseif obj.has_le_cs
                obj.delta_le_device=obj.le_device.delta;
            end
            
            p1=obj.xyz(:,1)';
            p2=obj.xyz(:,2)';
            p3=obj.xyz(:,3)';
            p4=obj.xyz(:,4)';
            
            if obj.c_te_device~=0
                if obj.te_device.is_tapered==1
                      t1=norm(p4-p1);
                      t2=norm(p3-p2);
                      p1=p1*(sum(obj.c_te_device)/t1)+p4*(1-sum(obj.c_te_device)/t1);
                      p2=p2*(sum(obj.c_te_device)*obj.c_t/obj.c_r/t2)+p3*(1-sum(obj.c_te_device)*obj.c_t/obj.c_r/t2);
%                     p1=p4+[-sum(obj.c_te_device)*cos(obj.Theta_r*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
%                     p2=p3+[-sum(obj.c_te_device)*cos(obj.Theta_t*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)]*obj.c_t/obj.c_r;
                else
                      t1=norm(p4-p1);
                      p1=p1*(sum(obj.c_te_device)/t1)+p4*(1-sum(obj.c_te_device)/t1);
                      p2=p2*(sum(obj.c_te_device)/t1)+p3*(1-sum(obj.c_te_device)/t1);
%                     p1=p4+[-sum(obj.c_te_device)*cos(obj.Theta_r*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
%                     p2=p3+[-sum(obj.c_te_device)*cos(obj.Theta_t*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
                end
                for i=1:length(obj.c_te_device)
                    %if obj.te_device.is_tapered==1
                    
                    hinge_vec=p2-p1;
                    
                    u=hinge_vec(1);
                    v=hinge_vec(2);
                    w=hinge_vec(3);
                    
                    x=p3(1);
                    y=p3(2);
                    z=p3(3);
                    
                    a=p1(1);
                    b=p1(2);
                    c=p1(3);
                    
                    Theta=obj.delta_te_device(1)*pi/180;
                    L=u^2+v^2+w^2;
                    
                    p3=[(a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*x*cos(Theta)+sqrt(L)*(-c*v+b*w-w*y+v*z)*sin(Theta);
                          (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*y*cos(Theta)+sqrt(L)*(c*u-a*w+w*x-u*z)*sin(Theta);
                          (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(Theta))+L*z*cos(Theta)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sin(Theta);]'/L;
                      
                    x=p4(1);
                    y=p4(2);
                    z=p4(3);  
                      
                      
                     p4=[(a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*x*cos(Theta)+sqrt(L)*(-c*v+b*w-w*y+v*z)*sin(Theta);
                          (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*y*cos(Theta)+sqrt(L)*(c*u-a*w+w*x-u*z)*sin(Theta);
                          (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(Theta))+L*z*cos(Theta)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sin(Theta);]'/L;
                      
                    obj.xyz_te_device(i,:,:)=[p1' p2' p3' p4'];
                    p1=p4;
                    p2=p3;
                end
            end
            
            if obj.c_le_device~=0
                
                obj.xyz_le_device=zeros(length(obj.c_le_device),3,4);
                
                p4=obj.xyz(:,1)'-[-sum(obj.c_le_device)*cos(obj.Theta_r*pi/180) -(sum(obj.c_le_device))*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) +(sum(obj.c_le_device))*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
                p3=obj.xyz(:,2)'-[-sum(obj.c_le_device)*cos(obj.Theta_t*pi/180) -(sum(obj.c_le_device))*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +(sum(obj.c_le_device))*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
                
                for i=length(obj.c_le_device):-1:1
                    dvec1=[(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*cos(obj.Theta_r*pi/180) (cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) -(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
                    dvec3=[(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*cos(obj.Theta_t*pi/180) (cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) -(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
                    dvec2=[obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*sin(obj.Theta_r*pi/180) -obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*sin(obj.dihed*pi/180) obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*cos(obj.Theta_r*pi/180)];
                    dvec4=[obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*sin(obj.Theta_t*pi/180) -obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*sin(obj.dihed*pi/180) obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*cos(obj.Theta_t*pi/180)];
                    
                    p2=p3-dvec3+dvec4;
                    p1=p4-dvec1+dvec2;
                    
                    obj.xyz_le_device(i,:,:)=[p1' p2' p3' p4'];
                    
                    p4=p1;
                    p3=p2;
                end
            end
        end
        
        function obj=compute_xyz_fixed(obj)
            if isempty(obj.c_te_device)&& isempty(obj.c_le_device)
                obj.xyz_fixed=obj.xyz;
            elseif isempty(obj.c_le_device)
                obj.xyz_fixed=[obj.xyz(:,1) obj.xyz(:,2) squeeze(obj.xyz_te_device(1,:,2))' squeeze(obj.xyz_te_device(1,:,1))'];
            elseif isempty(obj.c_te_device)
                obj.xyz_fixed=[squeeze(obj.xyz_le_device(end,:,4))' squeeze(obj.xyz_le_device(end,:,3))' obj.xyz(:,3) obj.xyz(:,4)];
            elseif ~isempty(obj.c_te_device)&&~isempty(obj.c_le_device)
                obj.xyz_fixed=[squeeze(obj.xyz_le_device(end,:,4))' squeeze(obj.xyz_le_device(end,:,3))'  squeeze(obj.xyz_te_device(1,:,2))' squeeze(obj.xyz_te_device(1,:,1))'];
            end
        end
        
        function obj = read_wingbox_coords(obj,pathNodeCoords,n,iNodes)
            span_grid=0:1/n:1;
            span_grid_aero=0:1/obj.n_span:1;
            front_coords=[];
            rear_coords=[];
            c4_coords=[];
            
            nodeCoords = importdata(pathNodeCoords);
            segmentNodeCoords = nodeCoords(iNodes:iNodes+n,:);
            
            for i=1:n+1
                 c4_coords(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            end
            
            obj.wingbox_coords=zeros(3,n+1,2);
            obj.wingbox_coords(:,:,1) = segmentNodeCoords';
            obj.wingbox_coords(:,:,2) = segmentNodeCoords';
            obj.wingbox_c4 = c4_coords;
%             
            for i=1:obj.n_span
                obj.c4_coords(:,i)=obj.xyz(:,4)+(obj.xyz(:,1)-obj.xyz(:,4))*0.75+span_grid_aero(i)*(obj.xyz(:,3)+(obj.xyz(:,2)-obj.xyz(:,3))*0.75-obj.xyz(:,4)-(obj.xyz(:,1)-obj.xyz(:,4))*0.75);
            end
        end
        
        function obj=compute_wingbox_coords(obj,frontspar,rearspar,n)
            %% TODO: generalize for n number of spars
            span_grid=0:1/n:1;
            span_grid_aero=0:1/obj.n_span:1;
            front_coords=[];
            rear_coords=[];
            c4_coords=[];
            
            front_sp=frontspar(1)*(1-span_grid)+frontspar(2)*span_grid;
            rear_sp=rearspar(1)*(1-span_grid)+rearspar(2)*span_grid;
            
            for i=1:n+1
                front_coords(:,i)=(obj.xyz(:,1)*(1-front_sp(i))+obj.xyz(:,4)*front_sp(i))+span_grid(i)*(obj.xyz(:,2)*(1-front_sp(i))+obj.xyz(:,3)*front_sp(i)-obj.xyz(:,1)*(1-front_sp(i))-obj.xyz(:,4)*front_sp(i));
                rear_coords(:,i)=(obj.xyz(:,1)*(1-rear_sp(i))+obj.xyz(:,4)*rear_sp(i))+span_grid(i)*(obj.xyz(:,2)*(1-rear_sp(i))+obj.xyz(:,3)*rear_sp(i)-obj.xyz(:,1)*(1-rear_sp(i))-obj.xyz(:,4)*rear_sp(i));
                c4_coords(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            end
            
            obj.wingbox_coords=zeros(3,n+1,2);
            obj.wingbox_coords(:,:,1)=front_coords(:,:);
            obj.wingbox_coords(:,:,2)=rear_coords(:,:);
            obj.wingbox_c4=c4_coords;
            
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            
            h_fr=zeros(1,length(front_sp));
            h_re=zeros(1,length(front_sp));
            
            for j=1:length(front_sp)
                profile_height_t=abs(interp1(coords_upper_t,profile_upper_t,front_sp(j),'lin','extrap')-interp1(coords_lower_t,profile_lower_t,front_sp(j),'lin','extrap'))*obj.c_t;
                profile_height_r=abs(interp1(coords_lower_r,profile_lower_r,front_sp(j),'lin','extrap')-interp1(coords_upper_r,profile_upper_r,front_sp(j),'lin','extrap'))*obj.c_r;
                h_fr(j)=profile_height_r*(1-span_grid(j))+profile_height_t*span_grid(j);
                
                profile_height_t=abs(interp1(coords_upper_t,profile_upper_t,rear_sp(j),'lin','extrap')-interp1(coords_lower_t,profile_lower_t,rear_sp(j),'lin','extrap'))*obj.c_t;
                profile_height_r=abs(interp1(coords_lower_r,profile_lower_r,rear_sp(j),'lin','extrap')-interp1(coords_upper_r,profile_upper_r,rear_sp(j),'lin','extrap'))*obj.c_r;
                h_re(j)=profile_height_r*(1-span_grid(j))+profile_height_t*span_grid(j);
            end
            obj.wingbox_height(:,1)=h_fr;
            obj.wingbox_height(:,2)=h_re;
            
            for i=1:obj.n_span
                obj.c4_coords(:,i)=obj.xyz(:,4)+(obj.xyz(:,1)-obj.xyz(:,4))*0.75+span_grid_aero(i)*(obj.xyz(:,3)+(obj.xyz(:,2)-obj.xyz(:,3))*0.75-obj.xyz(:,4)-(obj.xyz(:,1)-obj.xyz(:,4))*0.75);
            end
            
            
        end
        
        function obj=plot_segment(obj)
            
            if ~isempty(obj.te_device)
                for i=1:length(obj.c_te_device)
                    if ~obj.te_device.is_sym_defl
                        obj.te_device.delta=obj.te_device.delta_l_r(1);
                    end
                    obj=obj.compute_controlsurface_coordinates();
                    handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                    alpha(handle,0.4)
                    if obj.symmetric==1
                        if obj.te_device.is_sym_defl
                            handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',-squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                            alpha(handle,0.4)
                        else
%                             obj.te_device.delta
%                             obj.te_device.name
%                             zw=obj.te_device.delta(1);
                            obj.te_device.delta=obj.te_device.delta_l_r(2)*-1;
                            obj=obj.compute_controlsurface_coordinates();
                            handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',-squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                            alpha(handle,0.4)
                            obj.te_device.delta=obj.te_device.delta_l_r(1);
                            obj=obj.compute_controlsurface_coordinates();
                        end
                    end
                end
            end
            
            if ~isempty(obj.le_device)
                for i=1:length(obj.c_le_device)
                    handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                    alpha(handle,0.4)
                    if obj.symmetric==1
                        if obj.le_device.is_sym_defl
                            handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',-squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                            alpha(handle,0.4);
                        else
                            obj.te_device.delta=obj.te_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                            handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',-squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                            alpha(handle,0.4);
                            obj.te_device.delta=obj.te_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                        end
                    end
                end
            end
            
            handle=fill3(obj.xyz_fixed(1,:),obj.xyz_fixed(2,:),obj.xyz_fixed(3,:),'r');
            hold on
            alpha(handle,0.4)
            if obj.symmetric==1
                handle=fill3(obj.xyz_fixed(1,:),-obj.xyz_fixed(2,:),obj.xyz_fixed(3,:),'r');
                hold on
                alpha(handle,0.4)
            end
            axis equal
        end
        
        function pos=get_tip_c4point(obj)
            pos=obj.pos+[obj.b*tan(obj.c4_sweep*pi/180) obj.b*cos(obj.dihed*pi/180) +obj.b*sin(obj.dihed*pi/180)];
        end
        
        function c=get_tip_chord(obj)
            c=obj.c_t;
        end
        
        function obj=compute_forces(obj,panel_forces,panels,grid)
            n_pan=sum([obj.n_le_panels obj.n_chord  obj.n_te_panels]);
            obj.span_forces=zeros(3,obj.n_span);
            obj.span_moments=zeros(3,obj.n_span);
            span_grid=0.5/obj.n_span:1/obj.n_span:1-0.5/obj.n_span;
            for i=1:obj.n_span

                %zw=obj.c4_coords(:,i);
%                 hold on
%                                 plot3(zw(1,:),zw(2,:),zw(3,:),'o-');
%                 %                pause
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)-1+j;

                    obj.span_forces(:,i)=obj.span_forces(:,i)+panel_forces(:,idx)/(obj.b/obj.n_span);
                    r=grid(:,panels(4,idx))+(grid(:,panels(1,idx))-grid(:,panels(4,idx)))*0.75+0.5*(grid(:,panels(3,idx))+(grid(:,panels(2,idx))-grid(:,panels(3,idx)))*0.75-grid(:,panels(4,idx))-(grid(:,panels(1,idx))-grid(:,panels(4,idx)))*0.75);
%                     hold on
%                      fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'r');
                    obj.span_moments(:,i)=obj.span_moments(:,i)+cross(r-obj.c4_coords(:,i),panel_forces(:,idx))/(obj.b/obj.n_span);
                    % plot3(zw(1,:),zw(2,:),zw(3,:),'x-');
                end
            end
        end
        
        function grid=compute_deflected_grid_left(obj,panels,grid,deflections_structmesh,offset)
            span_grid=0:1/obj.n_span:1;
            n_pan=sum([obj.n_le_panels obj.n_chord  obj.n_te_panels]);
            
            le_ctr=1;
            le_loc_ctr=1;
            te_ctr=1;
            te_loc_ctr=1;
% %             
                deflections_structmesh(4,:)=-deflections_structmesh(4,:);
% %              % deflections_structmesh(5,:)=-deflections_structmesh(5,:);
          deflections_structmesh(6,:)=-deflections_structmesh(6,:);
    
%                            deflections_structmesh(4,:)=deflections_structmesh(4,:);
%              % deflections_structmesh(5,:)=-deflections_structmesh(5,:);
%                deflections_structmesh(6,:)=deflections_structmesh(6,:);   
%                
            for i=1:obj.n_span
                Theta=(obj.Theta_r*(1-span_grid(i))+obj.Theta_t*span_grid(i))*pi/180*0;
                
                a=-obj.dihed*pi/180;
                %a=obj.dihed*pi/180;
                c=obj.c4_sweep;
                
                Lx=[1       0       0
                    0   cos(a)  sin(a)
                    0   -sin(a) cos(a)];
                
                
                Lz=[cos(c) sin(c)   0
                    -sin(c) cos(c)  0
                    0           0   1];
                
                R=Lz*Lx*[0;Theta;0];

                c4_coords_edge(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
                c4_coords_edge(2,i)=-c4_coords_edge(2,i);
%                  hold on
%                  plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'x')
                c4_coords_edge(:,i)=c4_coords_edge(:,i)+deflections_structmesh(1:3,i);
                
            %     plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'ro')
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)+j-1+offset;
                    
                    grid(:,panels(2,idx))=grid(:,panels(2,idx))+deflections_structmesh(1:3,i);
                    
                    dist_x=norm(grid(1:3,panels(2,idx))-c4_coords_edge(1:3,i));
                    sgnx=sign(grid(1,panels(2,idx))-c4_coords_edge(1,i));
                    
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                    dx2=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    Theta=R(3);
                    dTheta=deflections_structmesh(6,i);
                    
                    dx3=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dy3=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    delta_twist_1=[dx2+dx3;-dy3;dz2];
                    
                    grid(:,panels(2,idx))=grid(:,panels(2,idx))+delta_twist_1;
                    
                   % plot3(grid(1,panels(2,idx)),grid(2,panels(2,idx)),grid(3,panels(2,idx)),'go')
                    if obj.has_le_cs
                        if le_ctr<=length(obj.n_le_panels)
                            le_loc_ctr=le_loc_ctr+1;
                            if le_loc_ctr>obj.n_le_panels(le_ctr)
                                le_loc_ctr=1;
                                le_ctr=le_ctr+1;
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                               % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                dTheta=deflections_structmesh(6,i);
                                Theta=R(3);
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_4=[dx2+dx3;-dy3;dz2];
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                            end
                            
                        end
                    end
                    
                    if obj.has_te_cs
                        if j>n_pan-sum(obj.n_te_panels)
                            if te_ctr<=length(obj.n_te_panels)
                                te_loc_ctr=te_loc_ctr+1;
                                if te_loc_ctr>obj.n_te_panels(te_ctr)
                                    te_loc_ctr=1;
                                    te_ctr=te_ctr+1;
                                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                                    dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                    sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                    dTheta=deflections_structmesh(5,i);
                                    Theta=R(2);
                                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    dTheta=deflections_structmesh(6,i);
                                    Theta=R(3);
                                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    delta_twist_4=[dx2+dx3;-dy3;dz2];
                                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                                end
                            end
                        elseif j==n_pan-sum(obj.n_te_panels)
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                           %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;-dy3;dz2];
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                        end
                    end
                end
                if  (isempty(obj.has_te_cs))||(obj.has_te_cs==0)
                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                    dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                    sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    dTheta=deflections_structmesh(6,i);
                    Theta=R(3);
                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    delta_twist_4=[dx2+dx3;-dy3;dz2];
                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                end
                le_ctr=1;
                te_ctr=1;
            end
            le_loc_ctr=1;
            te_loc_ctr=1;
            
            Theta=(obj.Theta_r*(1-span_grid(i+1))+obj.Theta_t*span_grid(i+1))*pi/180*0;
            
            a=-obj.dihed*pi/180;
            
            b=Theta;
            
            c=obj.c4_sweep;
            
            
            Lx=[1       0       0
                0   cos(a)  sin(a)
                0   -sin(a) cos(a)];
            
            
            
            Lz=[cos(c) sin(c)   0
                -sin(c) cos(c)  0
                0           0   1];
            
            R=Lz*Lx*[0;Theta;0];
            
            c4_coords_edge(:,i+1)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i+1)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            %c4_coords_edge(2,i)=-c4_coords_edge(2,i);
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'x')
            c4_coords_edge(:,i+1)=c4_coords_edge(:,i+1)+deflections_structmesh(1:3,i+1);
            
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'ro')
            for j=1:n_pan
                idx=obj.panel_start_idx+n_pan*(i-1)+j-1+offset;
                
                grid(:,panels(1,idx))=grid(:,panels(1,idx))+deflections_structmesh(1:3,i+1);
                
                dist_x=abs(grid(1,panels(1,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(1,idx))-c4_coords_edge(1,i+1));
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_2=[dx2+dx3;-dy3;dz2];
                
               % plot3(grid(1,panels(1,idx))-c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'gx')
                
                grid(:,panels(1,idx))=grid(:,panels(1,idx))+delta_twist_2;
                if obj.has_le_cs
                    if le_ctr<=length(obj.n_le_panels)
                        le_loc_ctr=le_loc_ctr+1;
                        if le_loc_ctr>obj.n_le_panels(le_ctr)
                            le_loc_ctr=1;
                            le_ctr=le_ctr+1;
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                         %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dTheta=deflections_structmesh(5,i+1);
                            Theta=R(2);
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            Theta=R(3);
                            dTheta=deflections_structmesh(6,i+1);
                            
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_3=[dx2+dx3;-dy3;dz2];
                            
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                        end
                    end
                end
                
                if obj.has_te_cs
                    if j>n_pan-sum(obj.n_te_panels)
                        if te_ctr<=length(obj.n_te_panels)
                            te_loc_ctr=te_loc_ctr+1;
                            if te_loc_ctr>obj.n_te_panels(te_ctr)
                                te_loc_ctr=1;
                                te_ctr=te_ctr+1;
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                            %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i+1);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                Theta=R(3);
                                dTheta=deflections_structmesh(6,i+1);
                                
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_3=[dx2+dx3;-dy3;dz2];
                                
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                            end
                        end
                    elseif j==n_pan-sum(obj.n_te_panels)
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;-dy3;dz2];
                        
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                    end
                end
            end
            if (isempty(obj.has_te_cs)) || (obj.has_te_cs==0)
                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_3=[dx2+dx3;-dy3;dz2];
                
                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
            end
            %grid(2,:)=-grid(2,:);
        end
        
        function grid=compute_deflected_grid(obj,panels,grid,deflections_structmesh)

            span_grid=0:1/obj.n_span:1;
            n_pan=sum([obj.n_le_panels obj.n_chord  obj.n_te_panels]);
            
            le_ctr=1;
            le_loc_ctr=1;
            te_ctr=1;
            te_loc_ctr=1;
            
            for i=1:obj.n_span
                %twist set to zerO?
                Theta=(obj.Theta_r*(1-span_grid(i))+obj.Theta_t*span_grid(i))*pi/180*0;
                
                a=obj.dihed*pi/180;
                c=obj.c4_sweep*pi/180;
                
                Lx=[1       0       0
                    0   cos(a)  sin(a)
                    0   -sin(a) cos(a)];
                
                
                Lz=[cos(c) sin(c)   0
                    -sin(c) cos(c)  0
                    0           0   1];
                
                %R is always zero beacuse twist is set to zero
                R=Lz*Lx*[0;Theta;0];
                
                c4_coords_edge(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
              %    hold on
              %    plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'x')
                c4_coords_edge(:,i)=c4_coords_edge(:,i)+deflections_structmesh(1:3,i);
                
         %        plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'ro')
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)+j-1;
                    
                    grid(:,panels(1,idx))=grid(:,panels(1,idx))+deflections_structmesh(1:3,i);
                    
                    dist_x=norm(grid(1:3,panels(1,idx))-c4_coords_edge(1:3,i));
                    sgnx=sign(grid(1,panels(1,idx))-c4_coords_edge(1,i));
                    
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                    dx2=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    Theta=R(3);
                    
                    dTheta=deflections_structmesh(6,i);
                    dx3=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dy3=sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    
                    
                    delta_twist_1=[dx2+dx3;dy3;dz2];
                    grid(:,panels(1,idx))=grid(:,panels(1,idx))+delta_twist_1;
                    
                  %    plot3(grid(1,panels(1,idx)),grid(2,panels(1,idx)),grid(3,panels(1,idx)),'go')
                    if obj.has_le_cs
                        if le_ctr<=length(obj.n_le_panels)
                            le_loc_ctr=le_loc_ctr+1;
                            if le_loc_ctr>obj.n_le_panels(le_ctr)
                                le_loc_ctr=1;
                                le_ctr=le_ctr+1;
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                dTheta=deflections_structmesh(6,i);
                                Theta=R(3);
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_4=[dx2+dx3;dy3;dz2];
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                            end
                            
                        end
                    end
                    
                    if obj.has_te_cs
                        if j>n_pan-sum(obj.n_te_panels)
                            if te_ctr<=length(obj.n_te_panels)
                                te_loc_ctr=te_loc_ctr+1;
                                if te_loc_ctr>obj.n_te_panels(te_ctr)
                                    te_loc_ctr=1;
                                    te_ctr=te_ctr+1;
                                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                                    dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                    sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                    dTheta=deflections_structmesh(5,i);
                                    Theta=R(2);
                                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    dTheta=deflections_structmesh(6,i);
                                    Theta=R(3);
                                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    delta_twist_4=[dx2+dx3;dy3;dz2];
                                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                                end
                            end
                        elseif j==n_pan-sum(obj.n_te_panels)
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                            % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;dy3;dz2];
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                        end
                    end
                end
                if  (isempty(obj.has_te_cs))||(obj.has_te_cs==0)
                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                    dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                    sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    dTheta=deflections_structmesh(6,i);
                    Theta=R(3);
                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    delta_twist_4=[dx2+dx3;dy3;dz2];
                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                end
                le_ctr=1;
                te_ctr=1;
            end
            % last spanwise row
            le_loc_ctr=1;
            te_loc_ctr=1;
            
            Theta=(obj.Theta_r*(1-span_grid(i+1))+obj.Theta_t*span_grid(i+1))*pi/180*0;
            
            a=obj.dihed*pi/180;
            
            b=Theta;
            
            c=obj.c4_sweep;
            
            
            Lx=[1       0       0
                0   cos(a)  sin(a)
                0   -sin(a) cos(a)];
            
            
            
            Lz=[cos(c) sin(c)   0
                -sin(c) cos(c)  0
                0           0   1];
            
            R=Lz*Lx*[0;Theta;0];
            
            c4_coords_edge(:,i+1)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i+1)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'x')
            c4_coords_edge(:,i+1)=c4_coords_edge(:,i+1)+deflections_structmesh(1:3,i+1);
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'ro')
            for j=1:n_pan
                idx=obj.panel_start_idx+n_pan*(i-1)+j-1;
                
                grid(:,panels(2,idx))=grid(:,panels(2,idx))+deflections_structmesh(1:3,i+1);
                
                dist_x=abs(grid(1,panels(2,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(2,idx))-c4_coords_edge(1,i+1));
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_2=[dx2+dx3;dy3;dz2];
                
                
                %plot3(grid(1,panels(2,idx))-c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'gx')
                
                grid(:,panels(2,idx))=grid(:,panels(2,idx))+delta_twist_2;
                if obj.has_le_cs
                    if le_ctr<=length(obj.n_le_panels)
                        le_loc_ctr=le_loc_ctr+1;
                        if le_loc_ctr>obj.n_le_panels(le_ctr)
                            le_loc_ctr=1;
                            le_ctr=le_ctr+1;
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                            dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                          %     plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dTheta=deflections_structmesh(5,i+1);
                            Theta=R(2);
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            Theta=R(3);
                            dTheta=deflections_structmesh(6,i+1);
                            
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_3=[dx2+dx3;dy3;dz2];
                            
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                        end
                    end
                end
                
                if obj.has_te_cs
                    if j>n_pan-sum(obj.n_te_panels)
                        if te_ctr<=length(obj.n_te_panels)
                            te_loc_ctr=te_loc_ctr+1;
                            if te_loc_ctr>obj.n_te_panels(te_ctr)
                                te_loc_ctr=1;
                                te_ctr=te_ctr+1;
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                             %    plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i+1);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                Theta=R(3);
                                dTheta=deflections_structmesh(6,i+1);
                                
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_3=[dx2+dx3;dy3;dz2];
                                
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                            end
                        end
                    elseif j==n_pan-sum(obj.n_te_panels)
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;dy3;dz2];
                        
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                    end
                end
            end
            if (isempty(obj.has_te_cs)) || (obj.has_te_cs==0)
                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                 % plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_3=[dx2+dx3;dy3;dz2];
                
                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
            end
        end
        
        function obj=mirror_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    
                    obj.te_device.delta=obj.te_device.delta*-1;
                end
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta=obj.le_device.delta*-1;
                end
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        function obj=right_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    obj.te_device.delta=obj.te_device.delta_l_r(2)*-1;
                end
                obj=obj.compute_controlsurface_coordinates();
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta=obj.le_device.delta_l_r(2)*-1;
                end
                obj=obj.compute_controlsurface_coordinates();
            end
            
        end
        
        function obj=left_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    obj.te_device.delta(1)=obj.te_device.delta_l_r(1);
                end
                obj=obj.compute_controlsurface_coordinates();
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta(1)=obj.le_device.delta_l_r(1);
                end
                obj=obj.compute_controlsurface_coordinates();
            end
        end
        
        function obj=compute_grid(obj,x_max,y_max,wake)
            
            grid3D=1;

            
            % number of spanwise panels
            obj.n_span=ceil(obj.b/y_max);
            % number of chordwise panels
            obj.n_chord=ceil(norm((obj.xyz_fixed(:,4)-obj.xyz_fixed(:,1))*0.5+(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,2))*0.5)/x_max);
            
            % compute span spacing vector
            %% TODO: optional use nonlinear cosine distribution of panels
            span_spacing=0:1/obj.n_span:1;
            
            % compute chordwise spacing vector
            chord_spacing_ctr=0:1/obj.n_chord:1;
            
            chord_spacing_le=[];
            chord_spacing_te=[];
            
            size_grid_le=0;
            size_grid_te=0;
            
            n_chordwise_points=obj.n_chord+1;
            n_chordwise_panels=obj.n_chord;
            
            obj.n_le_panels=[];
            obj.n_te_panels=[];
            
            if(obj.has_le_cs)
                n_le_devices=length(obj.xyz_le_device(:,1,1));
                for j=1:n_le_devices
                    obj.n_le_panels(j)=ceil(norm((obj.xyz_le_device(j,:,4)-obj.xyz_le_device(j,:,1))*0.5+(obj.xyz_le_device(j,:,3)-obj.xyz_le_device(j,:,2))*0.5)/x_max);
                end
                size_grid_le=length(span_spacing)*(sum(obj.n_le_panels)+length(obj.n_le_panels));
                n_chordwise_points=n_chordwise_points+sum(obj.n_le_panels)+length(obj.n_le_panels);
                n_chordwise_panels=n_chordwise_panels+sum(obj.n_le_panels);  
            end
            
            if(obj.has_te_cs)
                n_te_devices=length(obj.xyz_te_device(:,1,1));
                for j=1:n_te_devices
                    %% Problem: ceil is unstable 4.000000 will be 5
                    %% set proper limit
                    n=norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max;
                    if ceil(n)-n>0.9
                        obj.n_te_panels(j)=ceil(norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max)-1;
                    else
                        obj.n_te_panels(j)=ceil(norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max);
                    end
                end
                size_grid_te=length(span_spacing)*(sum(obj.n_te_panels)+length(obj.n_te_panels));
                n_chordwise_points=n_chordwise_points+sum(obj.n_te_panels)+length(obj.n_te_panels);
                n_chordwise_panels=n_chordwise_panels+sum(obj.n_te_panels);
            end
            
            grid_len=length(span_spacing)*length(chord_spacing_ctr)+size_grid_le+size_grid_te;
            grid=zeros(3,grid_len);
            grid_flat=zeros(3,grid_len);
            grid_upper=zeros(3,grid_len);
            grid_lower=zeros(3,grid_len);
            te_idx=zeros(1,grid_len);
            
            k=1;
            
            for i=1:length(span_spacing)
                % compute current chord
                c_c=obj.c_r*(1-span_spacing(i))+obj.c_t*span_spacing(i);
                % compute relative
                c_le_device_cs=[0 cumsum(obj.c_le_device)]/obj.c_r;
                
                relative_chord_pos_le=0;
                if(obj.has_le_cs)
                    for j=1:n_le_devices
                        % compute chord spacing
                        chord_spacing_le=0:1/ obj.n_le_panels(j):1;
                        r1=obj.xyz_le_device(j,:,1)+span_spacing(i)*(obj.xyz_le_device(j,:,2)-obj.xyz_le_device(j,:,1));
                        r2=obj.xyz_le_device(j,:,4)+span_spacing(i)*(obj.xyz_le_device(j,:,3)-obj.xyz_le_device(j,:,4));
                        for ii=1:length(chord_spacing_le)
                            if obj.le_device.is_tapered==1
                                relative_chord_pos_le=c_le_device_cs(j)+chord_spacing_le(ii)*obj.c_le_device(j)/obj.c_r;
                            else
                                relative_chord_pos_le=c_le_device_cs(j)*obj.c_r/c_c+chord_spacing_le(ii)*obj.c_le_device(j)/c_c;
                            end
                            skeleton_point=obj.compute_skeleton_point(relative_chord_pos_le,span_spacing(i));
                            grid(:,k)=r1+chord_spacing_le(ii)*(r2-r1);
                            grid_flat(:,k)=grid(:,k);
                            
                            if grid3D==1
                                [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos_le,span_spacing(i));
                                grid_upper(:,k)=grid(:,k);
                                grid_lower(:,k)=grid(:,k);
                                grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                                grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                                grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                                grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                            end
                            grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                            grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                            k=k+1;
                        end
                    end
                end
                
                r1=obj.xyz_fixed(:,1)+span_spacing(i)*(obj.xyz_fixed(:,2)-obj.xyz_fixed(:,1));
                r2=obj.xyz_fixed(:,4)+span_spacing(i)*(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,4));
                for j=1:length(chord_spacing_ctr)
                    
                    relative_chord_pos=sum(obj.c_le_device)/c_c+chord_spacing_ctr(j)*(c_c-sum(obj.c_te_device)-sum(obj.c_le_device))/c_c;
                    if obj.has_te_cs==1
                        if obj.te_device.is_tapered==1
                            %% careful only untapered te_device
                            relative_chord_pos=sum(obj.c_le_device)/c_c+chord_spacing_ctr(j)*(c_c-sum(obj.c_te_device)*c_c/obj.c_r-sum(obj.c_le_device))/c_c;
                        end
                    end
                    skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));
                    grid(:,k)=r1+chord_spacing_ctr(j)*(r2-r1);
                    grid_flat(:,k)=grid(:,k);
                    if grid3D==1
                       [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                       grid_upper(:,k)=grid(:,k);
                       grid_lower(:,k)=grid(:,k);
                       grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                       grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                       grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                       grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                    end
                    grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                    grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                    k=k+1;
                end
                
                if(obj.has_te_cs)
                    c_te_device_cs=[0 cumsum(obj.c_te_device)]/obj.c_r;
                    for j=1:n_te_devices
                        chord_spacing_te=0:1/obj.n_te_panels(j):1;
                        r1=obj.xyz_te_device(j,:,1)+span_spacing(i)*(obj.xyz_te_device(j,:,2)-obj.xyz_te_device(j,:,1));
                        r2=obj.xyz_te_device(j,:,4)+span_spacing(i)*(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,4));
                        for ii=1:length(chord_spacing_te)
                            if obj.te_device.is_tapered==1
                                relative_chord_pos=c_te_device_cs(j)+chord_spacing_te(ii)*sum(obj.c_te_device(j))/obj.c_r+sum(obj.c_le_device)/c_c+(c_c-sum(obj.c_te_device)*c_c/obj.c_r-sum(obj.c_le_device))/c_c;
                            else
                                relative_chord_pos=c_te_device_cs(j)*obj.c_r/c_c+chord_spacing_te(ii)*sum(obj.c_te_device(j))/c_c+sum(obj.c_le_device)/c_c+(c_c-sum(obj.c_te_device)-sum(obj.c_le_device))/c_c;
                            end
                            skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));  
                            grid(:,k)=r1+chord_spacing_te(ii)*(r2-r1);
                            grid_flat(:,k)=grid(:,k);
                            if grid3D==1
                                [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                                grid_upper(:,k)=grid(:,k);
                                grid_lower(:,k)=grid(:,k);
                                grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                                grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                                grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                                grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                            end
                            grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                            grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                            k=k+1;
                        end
                    end
                end
                te_idx(k-n_chordwise_points:k-1)=k-1;
                
                % if wake grid is desired save wake front line
                if wake==1
                    grid_wake(:,i)=grid(:,k-1);
                elseif wake==2
                    grid_wake(:,i)=grid(:,k-1);
                    %grid_wake(:,i)=grid(:,k-1)+0.20*231*2.7056e-04*(grid(:,k-1)-grid(:,k-2))/norm(grid(:,k-1)-grid(:,k-2));
                    %grid_wake(:,i)=grid(:,k-1)+0.25*(grid(:,k-1)-grid(:,k-2));
                end
            end
            
            x=[];
            
            % compute panel indices
            chordwise_pan=[obj.n_le_panels obj.n_chord  obj.n_te_panels];
            end_points=cumsum(chordwise_pan+1);
            offset=0;
            for ii=1:obj.n_span
                for jj=1:(length(end_points))
                    if jj==1
                        x=[x 1+offset:(end_points(jj)-1)+offset];
                    else
                        x=[x (end_points(jj-1)+1)+offset:(end_points(jj)-1)+offset];
                    end
                end
                offset=offset+sum(chordwise_pan+1);
            end
            
            panels=zeros(4,obj.n_span*n_chordwise_panels);
           
            % for potential flow solution
            obj.is_te=zeros(1,size(panels,2));
            
            for i=1:obj.n_span*n_chordwise_panels
                panels(:,i)=[x(i);x(i)+sum(chordwise_pan+1);x(i)+sum(chordwise_pan+1)+1;x(i)+1];
                if mod(i,n_chordwise_panels)==0
                    obj.is_te(i)=1;
                end
            end
            
            % if wake grid is desired
            if (wake==1)||(wake==2)
               len=length(span_spacing);
               for i=1:len-1
                   panels_wake(1,i)=i;
                   panels_wake(2,i)=i+1;
                   panels_wake(3,i)=i+1+len;
                   panels_wake(4,i)=i+len;
               end 
               obj.grid_wake=[grid_wake grid_wake];
               obj.panels_wake=panels_wake;
            end
            
            obj.grid=grid;
            obj.grid_vol_upper=grid_upper;
            obj.grid_vol_lower=grid_lower;
            obj.grid_flat=grid_flat;
            obj.panels=panels;
            obj.te_idx=te_idx;
        end
        
        function skeleton_point=compute_skeleton_point(obj,relative_chord_pos,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
           % skeleton_line_t=obj.c_t*0.5*(interp1(coords_lower_t,profile_lower_t,relative_chord_pos,'lin','extrap')+interp1(coords_upper_t,profile_upper_t,relative_chord_pos,'lin','extrap'));
            
            skeleton_line_t=obj.c_t*0.5*(nakeinterp1(coords_lower_t,profile_lower_t,relative_chord_pos)+nakeinterp1(coords_upper_t,profile_upper_t,relative_chord_pos));
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
           % skeleton_line_r=obj.c_r*0.5*(interp1(coords_lower_r,profile_lower_r,relative_chord_pos,'lin','extrap')+interp1(coords_upper_r,profile_upper_r,relative_chord_pos,'lin','extrap'));
            skeleton_line_r=obj.c_r*0.5*(nakeinterp1(coords_lower_r,profile_lower_r,relative_chord_pos)+nakeinterp1(coords_upper_r,profile_upper_r,relative_chord_pos));
            skeleton_point=skeleton_line_r*(1-relative_span_pos)+skeleton_line_t*relative_span_pos;
        end
        function skeletonAngle=compute_skeleton_angle(obj,relative_chord_pos,relative_chord_pos_0,relative_chord_pos_1,relative_span_pos)
            %relative_chord_pos_0: leading edge of the panel
            %relative_chord_pos_1: trailing edge of the panel
            
            % function to determine which panel of the discretized
            % camberline is crossed by the perpendicular segment of the chosen
            % panel (depending on approach) - detailed report has been
            % writen
            
            %tip initialization
            nPointsUpperTip=obj.profile_t(1,1);
            nPointsLowerTip=obj.profile_t(1,2);
            coordsUpperTip=obj.profile_t(2:1+nPointsUpperTip,1);
            profileUpperTip=obj.profile_t(2:1+nPointsUpperTip,2);
            coordsLowerTip=obj.profile_t(2+nPointsUpperTip:1+nPointsUpperTip+nPointsLowerTip,1);
            profileLowerTip=obj.profile_t(2+nPointsUpperTip:1+nPointsUpperTip+nPointsLowerTip,2);
            %root initialization
            nPointsUpperRoot=obj.profile_r(1,1);
            nPointsLowerRoot=obj.profile_r(1,2);
            coordsUpperRoot=obj.profile_r(2:1+nPointsUpperRoot,1);
            profileUpperRoot=obj.profile_r(2:1+nPointsUpperRoot,2);
            coordsLowerRoot=obj.profile_r(2+nPointsUpperRoot:1+nPointsUpperRoot+nPointsLowerRoot,1);
            profileLowerRoot=obj.profile_r(2+nPointsUpperRoot:1+nPointsUpperRoot+nPointsLowerRoot,2);
            %settings
            nSteps=100;
            coordsInterp=0:1/nSteps:1;
            %tip interpolation and skeleton line
            coordsUpperTipNorm=coordsUpperTip/(coordsUpperTip(end) - coordsUpperTip(1));
            profileUpperTipNorm=profileUpperTip/(coordsUpperTip(end) - coordsUpperTip(1));
            profileUpperTipInterp=interp1(coordsUpperTipNorm,profileUpperTipNorm,coordsInterp);
            
            coordsLowerTipNorm=coordsLowerTip/(coordsLowerTip(end) - coordsLowerTip(1));
            profileLowerTipNorm=profileLowerTip/(coordsLowerTip(end) - coordsLowerTip(1));
            profileLowerTipInterp=interp1(coordsLowerTipNorm,profileLowerTipNorm,coordsInterp);
                        
            skeletonLineTipInterp=0.5*(profileUpperTipInterp+profileLowerTipInterp);
            % skeletonAngleTip=atan(diff(skeletonLineTipInterp)/(1/nSteps));
            
            % Position of leading edge, trailing edge and 3/4 chord of the
            % panel
            y_relative_chord_pos_0_Tip=interp1(coordsInterp,skeletonLineTipInterp,relative_chord_pos_0);
            y_relative_chord_pos_1_Tip=interp1(coordsInterp,skeletonLineTipInterp,relative_chord_pos_1);
            y_relative_chord_pos=interp1([relative_chord_pos_0 relative_chord_pos_1],[y_relative_chord_pos_0_Tip y_relative_chord_pos_1_Tip],relative_chord_pos);
            Point_ref=[relative_chord_pos,y_relative_chord_pos];
            CamberLineTip=[coordsInterp',skeletonLineTipInterp']; % Coordinates of the camberline
            
            %Slope of the perpendicular segment (m_ref) on 3/4 chord of the panel
            if relative_chord_pos_0==relative_chord_pos_1
                m_ref=0;
            elseif y_relative_chord_pos_0_Tip==y_relative_chord_pos_1_Tip
                m_ref=1E10;
            else
                m=(y_relative_chord_pos_0_Tip - y_relative_chord_pos_1_Tip)/(relative_chord_pos_0 - relative_chord_pos_1);
                m_ref=(-1)/m;
            end
            
            % Find the crossed segment
            for iPointCamberLine=1:(size(CamberLineTip,1)-1)
                P_0=CamberLineTip(iPointCamberLine,:); % Point of the camberline
                P_1=CamberLineTip(iPointCamberLine+1,:); % Following point of the camberline
                
                yaux=y_relative_chord_pos+20; 
                
                Point_aux=[Point_ref(1)-((Point_ref(2)-yaux)/m_ref) yaux]; % Given point which has m_ref slope, and y coordinate higher than the highest y of the segment
                segment_1=Point_aux-Point_ref; % Segment perpendicular to the chosen panel
                segment_2=P_1-Point_ref; % Segment of the discretized camberline (from 3/4 chord to given point of the discretized camberline)
                segment_3=P_0-Point_ref; % Segment of the discretized camberline (from 3/4 chord to  a next given point of the discretized camberline)
                
                VecProd1=(segment_1(1,1)*segment_2(1,2))-(segment_1(1,2)*segment_2(1,1));
                VecProd2=(segment_1(1,1)*segment_3(1,2))-(segment_1(1,2)*segment_3(1,1));
                proof=VecProd1*VecProd2; % To find the panel crossed by the perpendicular segment, the third components of the cross product must have opposite directions (negative result) 
                
                if proof<=0
                    %                     if P_0(1,1)==P_1(1,1)
                    %                         relativeChordSkeletonAngleTip=0;
                    %                     elseif P_0(1,2)==P_1(1,2)
                    %                         relativeChordSkeletonAngleTip=pi/2;
                    %                     else
                    %                     m_Panel=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1));
                    %                     relativeChordSkeletonAngleTip=atan((-1) / m_Panel);
                    %                     end
                    
                    relativeChordSkeletonAngleTip=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1)); % Slope of the segment
                end
            end
            
            % The same implementation has been done to the root
            % camberline
            
            %root interpolation and skeleton line
            coordsUpperRootNorm=coordsUpperRoot/(coordsUpperRoot(end) - coordsUpperRoot(1));
            profileUpperRootNorm=profileUpperRoot/(coordsUpperRoot(end) - coordsUpperRoot(1));
            profileUpperRootInterp=interp1(coordsUpperRootNorm,profileUpperRootNorm,coordsInterp);
            
            coordsLowerRootNorm=coordsLowerRoot/(coordsLowerRoot(end) - coordsLowerRoot(1));
            profileLowerRootNorm=profileLowerRoot/(coordsLowerRoot(end) - coordsLowerRoot(1));
            profileLowerRootInterp=interp1(coordsLowerRootNorm,profileLowerRootNorm,coordsInterp);
                        
            skeletonLineRootInterp=0.5*(profileUpperRootInterp+profileLowerRootInterp);
            %             skeletonAngleRoot=atan(diff(skeletonLineRootInterp)/(1/nSteps));
            
            y_relative_chord_pos_0_Root=interp1(coordsInterp,skeletonLineRootInterp,relative_chord_pos_0);
            y_relative_chord_pos_1_Root=interp1(coordsInterp,skeletonLineRootInterp,relative_chord_pos_1);
            y_relative_chord_pos=interp1([relative_chord_pos_0 relative_chord_pos_1],[y_relative_chord_pos_0_Root y_relative_chord_pos_1_Root],relative_chord_pos);
            Point_ref=[relative_chord_pos,y_relative_chord_pos];
            CamberLineRoot=[coordsInterp',skeletonLineRootInterp'];
            
            if relative_chord_pos_0==relative_chord_pos_1
                m_ref=0;
            elseif y_relative_chord_pos_0_Root==y_relative_chord_pos_1_Root
                m_ref=1E10;
            else
                m=(y_relative_chord_pos_0_Root - y_relative_chord_pos_1_Root)/(relative_chord_pos_0 - relative_chord_pos_1);
                m_ref=(-1)/m;
            end
            
            for iPointCamberLine=1:(size(CamberLineRoot,1)-1)
                P_0=CamberLineRoot(iPointCamberLine,:);
                P_1=CamberLineRoot(iPointCamberLine+1,:);
                
                yaux=y_relative_chord_pos+20;
                
                Point_aux=[Point_ref(1)-((Point_ref(2)-yaux)/m_ref) yaux];
                segment_1=Point_aux-Point_ref;
                segment_2=P_1-Point_ref;
                segment_3=P_0-Point_ref;
                
                VecProd1=(segment_1(1,1)*segment_2(1,2))-(segment_1(1,2)*segment_2(1,1));
                VecProd2=(segment_1(1,1)*segment_3(1,2))-(segment_1(1,2)*segment_3(1,1));
                proof=VecProd1*VecProd2;
                
                if proof<=0
                    %                     if P_0(1)==P_1(1)
                    %                         relativeChordSkeletonAngleRoot=0;
                    %                     elseif P_0(2)==P_1(2)
                    %                         relativeChordSkeletonAngleRoot=pi/2;
                    %                     else
                    %                     m_Panel=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1));
                    %                     relativeChordSkeletonAngleRoot=atan((-1) / m_Panel);
                    %                     end
                    relativeChordSkeletonAngleRoot=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1)); % Slope of the segment
                end
            end
            
            %             relativeChordSkeletonAngleTip=interp1(coordsInterp(2:end),skeletonAngleTip,relative_chord_pos);
            %             relativeChordSkeletonAngleRoot=interp1(coordsInterp(2:end),skeletonAngleRoot,relative_chord_pos);
            
            skeletonAngle=relative_span_pos*relativeChordSkeletonAngleTip+(1-relative_span_pos)*relativeChordSkeletonAngleRoot;
            % plot profile:
%             figure
%             hold on
%             grid on
%             plot(coordsInterp,skeletonLineRootInterp)
%             plot(coordsLowerRootNorm,profileLowerRootNorm)
%             plot(coordsUpperRootNorm,profileUpperRootNorm)
        end
        
        function [lower_point,upper_point]=compute_thickness_point(obj,relative_chord_pos,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            lower_line_t=obj.c_t*nakeinterp1(coords_lower_t,profile_lower_t,relative_chord_pos);
            %lower_line_t=obj.c_t*interp1(coords_lower_t,profile_lower_t,relative_chord_pos,'lin','extrap');
            
            %upper_line_t=obj.c_t*interp1(coords_upper_t,profile_upper_t,relative_chord_pos,'lin','extrap');
            upper_line_t=obj.c_t*nakeinterp1(coords_upper_t,profile_upper_t,relative_chord_pos);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            %lower_line_r=obj.c_r*interp1(coords_lower_r,profile_lower_r,relative_chord_pos,'lin','extrap');
            lower_line_r=obj.c_r*nakeinterp1(coords_lower_r,profile_lower_r,relative_chord_pos);
            
            %upper_line_r=obj.c_r*interp1(coords_upper_r,profile_upper_r,relative_chord_pos,'lin','extrap');
            upper_line_r=obj.c_r*nakeinterp1(coords_upper_r,profile_upper_r,relative_chord_pos);
            
            lower_point=lower_line_r*(1-relative_span_pos)+lower_line_t*relative_span_pos;
            upper_point=upper_line_r*(1-relative_span_pos)+upper_line_t*relative_span_pos; 
        end
        
        function profile=get_profile(obj,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            %% TODO: profiles don#t always have same size
            profile_upper=profile_upper_r*(1-relative_span_pos)+profile_upper_t*relative_span_pos;
            profile_lower=profile_lower_r*(1-relative_span_pos)+profile_lower_t*relative_span_pos;
            
            coords_upper=coords_upper_r*(1-relative_span_pos)+coords_upper_t*relative_span_pos;
            coords_lower=coords_lower_r*(1-relative_span_pos)+coords_lower_t*relative_span_pos;
            
            profile=[[length(profile_upper) coords_upper' coords_lower']; [length(profile_lower) profile_upper' profile_lower']]';
        end
        
        function obj=write_tecplot_grid(obj,fileID,name,id)
            nxt=0;
            fprintf(fileID,'zone t="%sSEGMENT%i_UPPER", i=%i, j=1, k=1, f=point \n',name,id,obj.te_idx(1));
            if id==1
                i=1;
            else
                i=obj.te_idx(1)+1;
            end
            while i<=size(obj.grid_vol_upper,2)
                if nxt==1
                  i=obj.te_idx(end)-obj.te_idx(1)+1;  
                  fprintf(fileID,'zone t="%sSEGMENT%i_UPPER", i=%i, j=1, k=1, f=point \n',name,id+i,obj.te_idx(i)-i+1); 
                  nxt=0;
                end
                if (i==obj.te_idx(i))
                  nxt=1;
                end
              fprintf(fileID,'%f %f %f\n',obj.grid_vol_upper(1,i),obj.grid_vol_upper(2,i),obj.grid_vol_upper(3,i));         
              i=i+1;
            end
            nxt=0;
            fprintf(fileID,'zone t="%sSEGMENT%i_LOWER", i=%i, j=1, k=1, f=point \n',name,id,obj.te_idx(1)); 
            if id==1
                i=1;
            else
                i=obj.te_idx(1)+1;
            end
            while i<=length(obj.grid_vol_lower)
              if nxt==1
                i=obj.te_idx(end)-obj.te_idx(1)+1; 
                fprintf(fileID,'zone t="%sSEGMENT%i_LOWER", i=%i, j=1, k=1, f=point \n',name,id+i,obj.te_idx(i)-i+1); 
                nxt=0;
              end   
              if (i==obj.te_idx(i))
                  nxt=1;
              end
              fprintf(fileID,'%f %f %f\n',obj.grid_vol_lower(1,i),obj.grid_vol_lower(2,i),obj.grid_vol_lower(3,i));         
              i=i+1;
            end
        end
        
        function obj=read_xml_definition(obj,xmlstruct)
            if strcmp(xmlstruct.tag,'SEGMENT')
                obj.pos(1)=str2double(xmlstruct.child(1).child(1).value);
                obj.pos(2)=str2double(xmlstruct.child(1).child(2).value);
                obj.pos(3)=str2double(xmlstruct.child(1).child(3).value);
                obj.dihed=str2double(xmlstruct.child(2).value);
                if strcmp(xmlstruct.child(3).tag,'REAL_SWEEP')
                    obj.le_sweep=[];
                    obj.c4_sweep=[];
                    obj.real_sweep=str2double(xmlstruct.child(3).value);
                end
                obj.le_sweep=str2double(xmlstruct.child(3).value);
                obj.b=str2double(xmlstruct.child(4).value);
                obj.c_r=str2double(xmlstruct.child(5).value);
                obj.c_t=str2double(xmlstruct.child(6).value);
                obj.Theta_r=str2double(xmlstruct.child(7).value);
                obj.Theta_t=str2double(xmlstruct.child(8).value);
                obj.profile_name_r=xmlstruct.child(9).value;
                obj.profile_name_t=xmlstruct.child(10).value;
                
                obj.has_te_cs=0;
                obj.has_le_cs=0;
                if ~isempty(obj.le_sweep)
                obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
                end
                obj=obj.complete_params_from_stdinit();
                obj=obj.load_airfoils();
                obj=obj.compute_segment_coordinates();
                
            else
                fprintf('Unknown Data Format');
            end
        end
    end
end
