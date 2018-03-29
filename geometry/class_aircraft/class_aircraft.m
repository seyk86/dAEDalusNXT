classdef class_aircraft
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> aircraft name    
        name;
        
        %> fuselage object(s)
        fuselages;
        fuselage_structure_properties;
        %> wing object(s)
        wings;
        wings_structural_properties;
        
        nacelles;
        pylons;
        
        %
        engines;
        
        aero_surfaces;
        
        control_surfaces;
        control_deflections;

        S_wet_fus;
        S_wet_wing;
        S_wet_stabilizer;
        S_wet_engine;
        
        W_TO;       % Gross Take-off weight
        OWE;
        W_E;        % Empty weight
        W_F;        % Mission fuel weight
        P_TO;       % Take-off Power
        
        CL_max_L;
        CL_max_TO;
        CL_max;
        
        CD;
        CD_f;
        grid;
        grid_flat;
        grid_deflected;
        
        panels;
        te_idx;
        
        % volume grid
        grid_3D;  %flag if 3D grid required
        grid_vol;
        panels_vol;
        is_te;
        is_te_vol;
        te_idx_vol;
        opposite_te_vol;
        
        grid_wake;
        panels_wake;
        
        reference;

        weights;
        panel_to_beam_element;
        grid_settings
        J;
        
        control;
        acc;
        boundary_conditions;
    end
    
    methods
        %> function to read xml file
        obj=read_xml_definition(obj,filename);
        
        function obj = class_aircraft(wing,varargin)
            if nargin==1
                obj.CD_f=0;
                aero_surface(1)=wing;
                obj.wings=aero_surface;
                obj.name='blank';
            else
                obj=read_xml_definition(obj,wing);
            end
            obj.grid_settings=class_grid_settings;
            
           mkdir('results',obj.name);
        end
        
        function obj = add_aerosurface(obj,aero_surface)
            obj.aero_surfaces=[obj.aero_surfaces aero_surface];
        end
        
        function obj = add_wing(obj,wing)
            obj.wings=[obj.wings wing];
        end
        
        function obj = add_fuselage(obj,fuselage)
            obj.fuselages=[obj.fuselages fuselage];
        end
        
        function obj = add_nacelle(obj,nacelle)
            obj.nacelles=[obj.nacelles nacelle];
        end
        
        function obj=plot(obj)
            hold on
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i).plot(); 
            end
            for i=1:length(obj.wings)
                obj.wings(i).plot();
            end

            for i=1:length(obj.fuselages)
                obj.fuselages(i).plot_fuselage_surface();
            end
            
            for i=1:length(obj.nacelles)
                obj.nacelles(i).plot_fuselage_surface();
            end
        end
        
        function obj = compute_CD_f(obj,state,S_ref)
            obj.CD_f=0;
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i)=obj.aero_surfaces(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.aero_surfaces(i).CD_f*(1+obj.wings(i).symmetric);
            end
            
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).compute_friction_drag(state,S_ref);
                if obj.wings(i).symmetric==1
                    obj.CD_f=obj.CD_f+obj.wings(i).CD_f*2;
                else
                    obj.CD_f=obj.CD_f+obj.wings(i).CD_f;
                end
            end
            
            for i=1:length(obj.fuselages)
                obj.fuselages(i)=obj.fuselages(i).compute_wetted_area();
                obj.fuselages(i)=obj.fuselages(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.fuselages(i).CD_f;
            end
            
            for i=1:length(obj.nacelles)
                obj.nacelles(i)=obj.nacelles(i).compute_wetted_area();
                obj.nacelles(i)=obj.nacelles(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.nacelles(i).CD_f;
            end
        end
        
        function obj=compute_force_interpolation_matrix(obj,aircraft_structure)
            obj.panel_to_beam_element=zeros(size(obj.panels,2),15);
            for i=1:length(obj.wings)
                obj.panel_to_beam_element=obj.wings(i).compute_force_interpolation_matrix(obj.panel_to_beam_element,obj.panels,obj.grid_flat);
            end
            
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).T_matrix(obj.panel_to_beam_element,obj.grid,obj.panels,aircraft_structure.beam(i));
            end
        end

        function obj=compute_beam_forces(obj,F_body, M_body ,aircraft_structure,varargin)
            if nargin==5
                nw=varargin{1};
            else
                nw=length(obj.wings);
            end
            for i=1:nw
                obj.wings(i)=obj.wings(i).compute_beam_forces(obj.panel_to_beam_element,F_body, M_body ,aircraft_structure.beam(i));
            end
        end
        
        function obj=compute_forces(obj,F_aero,panels,grid) 
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).compute_c4_forces(F_aero,panels,grid);
            end
            
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i)=obj.aero_surfaces(i).compute_c4_forces(F_aero,panels,grid);
            end
        end
        
        function obj=compute_deflected_grid(obj,deflections_structmesh)
            x_max=obj.grid_settings.y_max_grid_size;
            obj.grid_deflected=obj.grid;
            for i=1:length(obj.wings)
                obj.grid_deflected=obj.wings(i).compute_deflected_grid(obj.panels,obj.grid_deflected,deflections_structmesh(i).def);
                 %obj.grid_deflected=obj.wings(i).grid_deflected;
            end
            if obj.grid_settings.aerodynamic_fuselage==1
                for j=i+1:i+length(obj.fuselages)
                    obj.grid_deflected=obj.fuselages(j-i).compute_deflected_grid_flat(obj.panels,obj.grid_deflected,deflections_structmesh(j).def,obj.fuselages(j-i).center_coords(3,1),x_max);
                end
                for k=j+1:j+length(obj.nacelles)
                    obj.grid_deflected=obj.nacelles(k-j).compute_deflected_grid_flat(obj.panels,obj.grid_deflected,deflections_structmesh(j+length(obj.nacelles)-(k-j)+1).def,obj.nacelles(k-j).center_coords(3,1),x_max);
                end
            end
        end
        
        function obj=compute_wingbox_coords(obj)
            if ~isempty(obj.wings_structural_properties)
                for i=1:length(obj.wings)
                    obj.wings(i)=obj.wings(i).compute_wingbox_coords(obj.grid_settings.dy_max_struct_grid,obj.wings_structural_properties(i).frontspar,obj.wings_structural_properties(i).rearspar);
                end
            else
                for i=1:length(obj.wings)
                    obj.wings(i)=obj.wings(i).compute_wingbox_coords(obj.grid_settings.dy_max_struct_grid);
                end
            end
        end
        
        function obj=compute_shell_coords(obj)
            for i=1:length(obj.fuselages)
                obj.fuselages(i)=obj.fuselages(i).compute_shell_coords(obj.grid_settings.y_max_grid_size);
            end
        end
        
        function obj=compute_acceleration(obj,wingaero,wingstructure,state)
            
            CX=wingaero.CX;
            CZ=wingaero.CZ;
            CY=wingaero.CY;
            CL=wingaero.CL;
            CM=wingaero.CM;
            CN=wingaero.CN;
            rho=wingaero.state.rho_air;
            V=norm(wingaero.Uinf);
            Sref=wingaero.reference.S_ref;
            bref=wingaero.reference.b_ref;
            cref=wingaero.reference.c_ref;
            
            m=wingstructure.f_compute_totalMass;
            
            J=wingstructure.f_compute_moment_of_inertia(wingaero.reference.p_ref);
            J1=J(1,1);%state.aircraft_state.I_xyz(1,1);
            J2=J(2,2);%state.aircraft_state.I_xyz(2,2);
            J3=J(3,3);%state.aircraft_state.I_xyz(3,3);
            
            obj.J=[J1 J2 J3];
            ax=1/m*(1/2*rho*V^2*CX*Sref);
            ay=1/m*(1/2*rho*V^2*CY*Sref);
            az=1/m*(1/2*rho*V^2*CZ*Sref);
            ap=1/J1*(1/2*rho*V^2*CL*Sref*bref);
            aq=1/J2*(1/2*rho*V^2*CM*Sref*cref);
            ar=1/J3*(1/2*rho*V^2*CN*Sref*bref);
            
            obj.acc=[ax ay az ap aq ar]; 
        end

        function obj=compute_grid(obj)
            % read desired grid size from grid_settings structure
            x_max=obj.grid_settings.x_max_grid_size;
            y_max=obj.grid_settings.y_max_grid_size;
            wake=obj.grid_settings.wake;

            % initialize variables for 2D grid
            grid=[];
            grid_flat=[];
            te_idx=[];
            panels=[];
            
            % initialize variables for 3D grid
            obj.grid_vol=[];
            obj.panels_vol=[];
            obj.is_te_vol=[];
            obj.is_te=[];
            obj.te_idx_vol=[];
            obj.opposite_te_vol=[];
            
            obj.grid_wake=[];
            obj.panels_wake=[]; 
            
            % assemble grid from all wings
            for i=1:length(obj.wings)
                % first compute all sub grids
                obj.wings(i)=obj.wings(i).compute_grid(x_max,y_max,wake);
                
                % all operations for 2D grid
                grid_len_b4=length(grid);
                grid=[grid obj.wings(i).grid];
                grid_flat=[grid_flat obj.wings(i).grid_flat];
                panel_len_b4=length(panels);
                panels=[panels,obj.wings(i).panels+grid_len_b4];
                te_idx=[te_idx,obj.wings(i).te_idx+grid_len_b4];
                obj.wings(i).grid_start_idx=grid_len_b4+1;
                obj.wings(i).panel_start_idx=panel_len_b4+1;
                
                obj.wings(i)=obj.wings(i).update_idx();

                obj.is_te=[obj.is_te obj.wings(i).is_te_vol(1:end/2)];
                % all operations for 3D grid
                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.wings(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.wings(i).panels_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.wings(i).is_te_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.wings(i).te_idx_vol+grid_len_b4_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.wings(i).opposite_te_vol+panel_len_b4_vol];
                obj.wings(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.wings(i).panel_start_idx_vol=panel_len_b4_vol+1;
                
                if (wake==1)||(wake==2)
                    grid_len_wake_b4_vol=length(obj.grid_wake);
                    obj.grid_wake=[obj.grid_wake obj.wings(i).grid_wake];
                    obj.panels_wake=[obj.panels_wake obj.wings(i).panels_wake+grid_len_wake_b4_vol];
                end
            end
            
            obj.grid=grid;
            obj.grid_flat=grid_flat;
            obj.panels=panels;
            obj.te_idx=te_idx;
            if obj.grid_settings.aerodynamic_fuselage==1
                obj=compute_fuselage_grid(obj);
            end
        end    
        
                
        function obj=f_update_grid(obj)
            
        end

        function obj=compute_fuselage_grid(obj)
            
            % read desired grid size from grid_settings structure
            x_max=obj.grid_settings.y_max_grid_size;
            y_max=obj.grid_settings.x_max_grid_size;
            wake=obj.grid_settings.wake;
            
            % initialize variables for 2D grid
           % grid=obj.grid_vol;
            grid_flat=obj.grid;
            te_idx=obj.te_idx_vol;
           % panels=obj.panels_vol;
            te_idx_flat=obj.te_idx;
            panels_flat=obj.panels;
            len=[];
            
            for i=1:length(obj.fuselages)
                
                obj.fuselages(i)=obj.fuselages(i).compute_shell_coords(x_max);
                obj.fuselages(i)=obj.fuselages(i).compute_grid();
                obj.fuselages(i)=obj.fuselages(i).compute_grid_flat(obj.fuselages(i).center_coords(3,1),x_max);
                
                % all operations for 2D grid
                grid_len_b4=length(grid_flat);
                panel_len_b4=length(panels_flat);
                grid_flat=[grid_flat obj.fuselages(i).grid_flat];
                panels_flat=[panels_flat,obj.fuselages(i).panels_flat+grid_len_b4];
                te_idx_flat=[te_idx_flat,obj.fuselages(i).te_idx_flat+grid_len_b4];
                obj.fuselages(i).grid_start_idx_flat=grid_len_b4+1;
                obj.fuselages(i).panel_start_idx_flat=panel_len_b4+1; 
                
                
                obj.fuselages(i)=obj.fuselages(i).update_idx();
%                 grid_len_b4=length(grid);
%                 grid=[grid obj.fuselages(i).grid_flat];
%                 panel_len_b4=length(panels);
%                 panels=[panels,obj.fuselages(i).panels+grid_len_b4];
%                 te_idx=[te_idx,obj.fuselages(i).te_idx+grid_len_b4];
%                 obj.fuselages(i).grid_start_idx_vol=grid_len_b4+1;
%                 obj.fuselages(i).panel_start_idx=panel_len_b4+1; 
                %obj.wings(i)=obj.wings(i).update_idx();

                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.fuselages(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.fuselages(i).panels_vol+grid_len_b4_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.fuselages(i).te_idx_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.fuselages(i).is_te_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.fuselages(i).opposite_te_vol+panel_len_b4_vol];
                obj.fuselages(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.fuselages(i).panel_start_idx_vol=panel_len_b4_vol+1;
            end
             
            for i=1:length(obj.nacelles)
                obj.nacelles(i)=obj.nacelles(i).compute_shell_coords(x_max);
                obj.nacelles(i)=obj.nacelles(i).compute_grid();
                obj.nacelles(i)=obj.nacelles(i).compute_grid_flat(obj.nacelles(i).center_coords(3,1),x_max);
                
                grid_len_b4=length(grid_flat);
                grid_flat=[grid_flat obj.nacelles(i).grid_flat];
                panels_flat=[panels_flat,obj.nacelles(i).panels_flat+grid_len_b4];
                te_idx_flat=[te_idx_flat,obj.nacelles(i).te_idx_flat+grid_len_b4];
                obj.nacelles(i).grid_start_idx_flat=grid_len_b4+1;
                obj.nacelles(i).panel_start_idx_flat=panel_len_b4+1; 
                
                obj.nacelles(i)=obj.nacelles(i).update_idx();

                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.nacelles(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.nacelles(i).panels_vol+grid_len_b4_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.nacelles(i).te_idx_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.nacelles(i).is_te_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.nacelles(i).opposite_te_vol+panel_len_b4_vol]; 
                obj.nacelles(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.nacelles(i).panel_start_idx_vol=panel_len_b4_vol+1;
                
            end 
            obj.grid=grid_flat;
            obj.grid_flat=grid_flat;
            obj.panels=panels_flat;
            obj.te_idx=te_idx_flat;
        end
        
        function obj=write_modes_in_tecplot(obj,aircraft_structure,nmodes,exaturation_factor)
            for i=1:length(obj.control_surfaces)
                obj=obj.f_set_control_surface(obj.control_surfaces{i},0);
            end
            obj=obj.compute_grid();
            for mod=1:nmodes
                for omega=[0 pi/2 3*pi/2]
                    aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,mod)*exaturation_factor*sin(omega);
                    aircraft_structure=aircraft_structure.f_postprocess();
                    obj=obj.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    obj.write_grid_deflected(['results/' obj.name '/mode_' num2str(mod) '_' num2str(omega)],[0 0 0]',[0 0 0]',[0 0 0]')
                    %obj.plot_grid_deflected(['results/' obj.name '/mode_' num2str(mod) '_' num2str(omega)]);
                end
            end
        end
        
        function obj=write_modes_in_paraview(obj, aircraft_structure,Modes,exagFactor,folder, varargin)
            if ~isempty(varargin)
                shape=varargin{1};
            else
                shape=[];
            end
                mkdir(folder);
            for i=Modes
                aircraft_structure.nodal_deflections=shape+exagFactor*aircraft_structure.modeshapes(:,i);
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=obj.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
                aircraft.write_grid_deflected([folder, '\modeshape_', sprintf('%.3i',i)],[0 0 0]',[0 0 0]',0);
            end                
        end
        
        function obj=write_mode_animation_in_paraview(obj,aircraft_structure,Modes,exagFactor,folder,frames,varargin)
            if ~isempty(varargin)
                shape=varargin{1};
            else
                shape=aircraft_structure.modeshapes(:,1)*0;
            end
            inputExagFactor=exagFactor;
            exagFactor=ones(1,max(Modes));
            if length(inputExagFactor)==length(Modes)
                exagFactor(Modes)=inputExagFactor;
            end
            if length(inputExagFactor)==1
                exagFactor=ones(1,max(Modes))*inputExagFactor;
            end
            mkdir(folder);
            for i=Modes
                for j=1:frames
                    aircraft_structure.nodal_deflections=shape+cos((j/frames)*2*pi)*exagFactor(i)*aircraft_structure.modeshapes(:,i);
                    aircraft_structure=aircraft_structure.f_postprocess();
                    aircraft=obj.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    aircraft.write_grid_deflected([folder, '\modeshape_', sprintf('%.3i',i),'_frame_', int2str(j)],[0 0 0]',[0 0 0]',0);
                end
            end
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
        
        function obj=write_tecplot_volgrid(obj,filename)
            mode='W';
            fileID = fopen([filename '.tp'],mode);
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid_vol),length(obj.panels_vol));
            
            for i =1:length(obj.grid_vol)
                fprintf(fileID,'%f %f %f %f\n',obj.grid_vol(1,i),obj.grid_vol(2,i),obj.grid_vol(3,i),1);
            end
            
            for i =1:length(obj.panels_vol)
                fprintf(fileID,'%i %i %i %i \n',obj.panels_vol(1,i),obj.panels_vol(2,i),obj.panels_vol(3,i),obj.panels_vol(4,i));
            end

            fclose(fileID);
        end
        
        function obj=write_tecplot_wingbox(obj,filename)

            mode='W';

            fileID = fopen([filename '.tp'],mode);
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            for nwing=1:length(obj.wings)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,1),obj.wings(nwing).wingbox_coords(2,i,1),obj.wings(nwing).wingbox_coords(3,i,1)+obj.wings(nwing).wingbox_height(i,1)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,1),obj.wings(nwing).wingbox_coords(2,i,1),obj.wings(nwing).wingbox_coords(3,i,1)-obj.wings(nwing).wingbox_height(i,1)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,2),obj.wings(nwing).wingbox_coords(2,i,2),obj.wings(nwing).wingbox_coords(3,i,2)+obj.wings(nwing).wingbox_height(i,2)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,2),obj.wings(nwing).wingbox_coords(2,i,2),obj.wings(nwing).wingbox_coords(3,i,2)-obj.wings(nwing).wingbox_height(i,2)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                
            end
            fclose(fileID);
        end
        
        
        function obj=write_tecplot_grid(obj,filename)
            fileID = fopen(filename,'W');
            fprintf(fileID,'TITLE = "Aircraft"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z"\n');
            
            for i=1:length(obj.wings)
                obj.wings(i).write_tecplot_grid(fileID);
            end
            fclose(fileID);
        end
        
        function obj=plot_grid_vol(obj)
            hold on
            for i=1:length(obj.wings)
                obj.wings(i).plot_grid_vol();
            end
            for i=1:length(obj.fuselages)
                obj.fuselages(i).plot_grid_vol();
            end
            for i=1:length(obj.nacelles)
                obj.nacelles(i).plot_grid_vol();
            end
            axis equal
            axis tight
            grid on
        end
        
        function obj=plot_wingbox_coordinates(obj)
           figure 
           hold on
           for i=1:length(obj.wings)
              plot3(obj.wings(i).wingbox_coords(1,:,1),obj.wings(i).wingbox_coords(2,:,1) ,obj.wings(i).wingbox_coords(3,:,1),'-rx')
              plot3(obj.wings(i).wingbox_coords(1,:,2),obj.wings(i).wingbox_coords(2,:,2) ,obj.wings(i).wingbox_coords(3,:,2),'-bx')
           end
        end
               
        function obj=plot_grid_deflected(obj,filename)
            hold on
            nel_fus=0;
            
%             if ~isempty(obj.fuselages)
%                 for i=1:length(obj.fuselages)
%                     nel_fus=nel_fus+length(obj.fuselages(i).panels);
%                 end
%             end
            for i=1:length(obj.panels)
                handle= fill3(obj.grid_deflected(1,obj.panels(:,i)), obj.grid_deflected(2,obj.panels(:,i)),obj.grid_deflected(3,obj.panels(:,i)),'b');
                alpha(handle,0.4);
            end
            axis equal
            axis tight
            grid on
    
        end
        
        function obj=write_grid_deflected(obj,filename,pos,Euler,ref)
             mode='W';
             pos(1)=0;
                             Lx=[1       0       0
                    0   cos(Euler(1))  -sin(Euler(1))
                    0   sin(Euler(1)) cos(Euler(1))];
                
                Ly=[cos(Euler(2)) 0 -sin(Euler(2))
                    0      1    0
                    sin(Euler(2))  0   cos(Euler(2))];
                
                Lz=[cos(Euler(3)) sin(Euler(3))   0
                    -sin(Euler(3)) cos(Euler(3))  0
                    0           0   1];
                
                
                M_BI=Lz*Ly*Lx;
             
            fileID = fopen([filename '.tp'],mode);
            append=0;
            if append~=1
                fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
                fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            for i =1:length(obj.grid)
                ri=obj.grid_deflected(:,i)-ref;
                pt=M_BI'*ri+pos(1:3);
                fprintf(fileID,'%f %f %f %f\n',pt(1),pt(2),pt(3),0);
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end
            fclose(fileID);
        end

        function obj=f_set_control_surface(obj,name,deflection) 
            for i=1:length(obj.wings)
                for j=1:length(obj.wings(i).wing_segments)
                    if  ~isempty(obj.wings(i).wing_segments(j).te_device)
                        if strcmp(obj.wings(i).wing_segments(j).te_device.name,name)
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(name,deflection);
                        end
                        
                        if strcmp([obj.wings(i).wing_segments(j).te_device.name '_left'],name)
                            obj.wings(i).wing_segments(j).te_device.delta(1)=deflection;
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(obj.wings(i).wing_segments(j).te_device.name,-deflection,'right');
                        end
                        
                        if strcmp([obj.wings(i).wing_segments(j).te_device.name '_right'],name)
                            obj.wings(i).wing_segments(j).te_device.delta(2)=deflection;
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(obj.wings(i).wing_segments(j).te_device.name,-deflection,'left');
                        end
                        
                    elseif ~isempty(obj.wings(i).wing_segments(j).le_device)
                        if strcmp(obj.wings(i).wing_segments(j).le_device.name,name)
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(name,deflection);
                        end
                        
                        if strcmp([obj.wings(i).wing_segments(j).le_device.name '_left'],name)
                            
                            obj.wings(i).wing_segments(j).le_device.delta(1)=deflection;
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(obj.wings(i).wing_segments(j).le_device.name,deflection);
                        end
                        
                        if strcmp([obj.wings(i).wing_segments(j).le_device.name 'right'],name)
                            obj.wings(i).wing_segments(j).le_device.delta(2)=deflection;
                            obj.wings(i).wing_segments(j)=obj.wings(i).wing_segments(j).f_deflect_control_surface(obj.wings(i).wing_segments(j).le_device.name,deflection);
                        end
                        
                    end
                end
            end
            for i=1:length(obj.control_surfaces)
                if(strcmp(obj.control_surfaces{i},name))
                   obj.control_deflections{i}=deflection; 
                end
            end
        end
        
        function [] = check_panel_to_beam_element(obj, aircraft_structure,iPanel)
            %% panel to beam element checks
            figure
            hold on
            grid on
            axis equal
            % which wing
            for i=1:size(obj.wings,2) 
                if iPanel<=obj.wings(i).panel_start_idx+size(obj.wings(i).panels,2)
                    wingId=i;
                    break;
                end
            end
            
            %aircraft_ow.plot_grid;
            if iPanel-30<obj.wings(wingId).panel_start_idx;
                startPanel=obj.wings(wingId).panel_start_idx;
            else
                startPanel=iPanel-30;
            end
            if iPanel+30>obj.wings(wingId).panel_start_idx+size(obj.wings(i).panels,2)
                endPanel=obj.wings(wingId).panel_start_idx+size(obj.wings(i).panels,2)-1;
            else
                endPanel=iPanel+30;
            end
            for i=startPanel:endPanel
                h=fill3(obj.grid(1,obj.panels(:,i)),obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                set(h,'facealpha',.25)
            end
            h=fill3(obj.grid(1,obj.panels(:,iPanel)),obj.grid(2,obj.panels(:,iPanel)),obj.grid(3,obj.panels(:,iPanel)),'r');

            if ~obj.panel_to_beam_element(iPanel,1)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,1),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,1)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,2),3))
            end

            if ~obj.panel_to_beam_element(iPanel,6)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,6),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,6)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,7),3))
            end

            if ~obj.panel_to_beam_element(iPanel,11)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,11),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,11)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,12),3))
            end
        end
        function [] = check_panel_to_beam_element2(obj, aircraft_structure,iBeam,iBe)

            %determine panel range for  wingidx for iBe
            panelRange=[obj.wings(iBeam).panel_start_idx:obj.wings(iBeam).panel_start_idx+size(obj.wings(iBeam).panels,2)-1]';

            %create matrix with panelIds belonging to iBe and Arearatio of Panel
            ptba=[panelRange obj.panel_to_beam_element(panelRange,1:2); panelRange obj.panel_to_beam_element(panelRange,6:7);panelRange obj.panel_to_beam_element(panelRange,11:12)];
            ptba2=ptba(ptba(:,2)==iBe,1:3);

            %plot beamelement
            figure; 
            hold on;
            for i=iBe-5:iBe+5
                 if ~(i>aircraft_structure.beam(iBeam).nel) && i>0
                     P1=aircraft_structure.beam(iBeam).node_coords(i,:);
                     P2=aircraft_structure.beam(iBeam).node_coords(i+1,:);
                     scatter3(P1(1),P1(2),P1(3),'bs');
                     scatter3(P2(1),P2(2),P2(3),'bs');
                     Pm=(P1+P2)/2;
                     text(Pm(1),Pm(2),Pm(3), num2str(i))
                 end
            end
            %check if panels are missing
            missing=setdiff(min(ptba2(:,1)):max(ptba2(:,1)), ptba2(:,1))';
            ptba2=[ptba2; [missing zeros(length(missing),1),  zeros(length(missing),1)]];
            %plotting also couple of panels before minNo and after maxNo
            n=10;
            ptba2=[ptba2; [max(ptba2(:,1))+1:max(ptba2(:,1))+n ; zeros(1,n); zeros(1,n)]'];
            ptba2=[ptba2; [min(ptba2(:,1))-n:min(ptba2(:,1)-1) ; zeros(1,n); zeros(1,n)]'];

            for i=1:size(ptba2,1)
                 if any(panelRange==ptba2(i,1))
                     %plot panel and fill according to area ratio
                     h=fill3(obj.grid(1,obj.panels(:,ptba2(i,1))),obj.grid(2,obj.panels(:,ptba2(i,1))),obj.grid(3,obj.panels(:,ptba2(i,1))),'g');
                     set(h,'facealpha',ptba2(i,3))
                 end
            end
            axis equal;
        end

    end
end
