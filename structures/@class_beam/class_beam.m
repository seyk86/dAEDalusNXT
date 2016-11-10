%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%

%> @file class_beam.m
%> @brief File containing the class for a finite element beam
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
%> @brief base class for a finite element beam
%> This class is the base class for a finite element beam, the follwing
%> coordinate system is used
%> Coordinate Definition (according to Airbus Wing Geometry Definitions ):
%>   y= along starboard wing
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
% ======================================================================

classdef class_beam < matlab.mixin.Heterogeneous
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class Attributes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        %> beam is being self-designed (=0) or given from external FEM (=1)
        isExternalFEM=0;
        %> external beam FEM mass matrix and stiffness matrix
        externalFEM = struct('Mext',[],'Kext',[]);
        %> name of the beam object
        beam;
        %> name of the beam object
        identifier;                     
        %> number of elements   
        nel;   
        %> array of beam elements
        beamelement;    
        %> DOF's for Finite Element (6 DOF implemented)
        el_ndof=6;           
        %> number of system DOF's
        ndof;
        %> number of free system DOF's
        free_ndof;                      
        %> number of boundary conditions
        n_bc=0;     
        %> array of boundary conditions
        boundary_condition;             
        %> coordinates of first node in global coordinate system
        r_Ref=[0 0 0];                  
        %> nodal coordinates in global coordinate system, original state
        node_coords;  
        %> nodal coordinates in global coordinate system in deflected state
        node_coords_deflected;          
        %> concentrated external nodal forces in x,y,z axis
        nodal_forces;       
        %> concentrated external nodal moments around x,y,z axis
        nodal_moments;   
        %> distributed external nodal loads fx,fy,fz,mx,my,mz
        nodal_loads;               
        %> concentrated external nodal loads fx,fy,fz,mx,my,mz
        nodal_loads_def;                   
        %> concentrated eccentric nodal masses m,x,y,z,I1,I2,I3 
        nodal_masses;                   
        %> flag if element load vector needs to be updated before solving
        update_Q=1;         
        %> flag if element stiffness matrix to be updated before solving
        update_K=1;                     
        %> flag if elemet mass matrix is to be updated before soving
        update_M=1;
        %> system stiffness matrix
        K;   
        %> system Mass matrix
        M;
        %> system Lumped mass matrix
        M_lumped;
        %> system load vector 
        Q; 
        %> system deflection vector
        w;      
        %> structural system matrix
        Kff;
        %> mass system matrix 
        Mff; 
        
        Mff_lumped;
        
        Mfp; %coupling matrix between free and prescribed DOF's
       
        Mpp; %Mass matrix for prescribed DOF's
     
        %> coupling matrix between free and prescribed DOF's
        Kfp;     
        %> stiffness matrix for prescribed DOF's
        Kpp;                            
        %> free deformations
        wf; 
        %> prescribed deformations
        wp;                             
        Ff;
        Fp;
        %> remaining forces
        Frem;                 
        %> vector for sorting free and prescribed DOF's in the matrices
        sort_vec;                       
        %> internal force
        P;              
        %> internal forces at free nodes
        Pf;                            
        %> error for non linear loop
        err;                            
        Kt;
        
        %> store information about what position in Kff equals which DOF and node
        Kff_node_dof_info;               
        sort_vec_f;
        
        %> structural stiffness matrix
        Ks;
        %> structural load vector 
        Qs;         
        %> vector of nodal_deflections in all DOF's in global coordinates
        nodal_deflections; 
        %> vector of nodal_deflections in all DOF's in local coordinates
        nodal_deflections_loc;          
        %> reaction forces at clamped positions
        reaction_forces;                
        
        %> flag identifying symmetric structure
        is_sym;
        
        %% general values
        
        %> structural mass
        m_total=0.0;
        %> internal forces in global coordinates
        node_loadings=0.0;
        %> internal forces in local coordinates
        node_loadings_loc=0.0;
        
        %% environmental factors
        %> graviational acceleration
        g=9.81;  
        %> load factor
        load_factor=1;   
        %> pitch angle
        gamma=0;
        
        %>  fuel density,put into beam derivate wing later
        fuel_density=807.5;   
        %> to identify the current loadcase
        loadcase_index; 
        
        
        %% dynamic solution results
        modeshapes;
        modefrequencies;
        
        
        %% solver settings class
        settings
        
        
        T_sort;
        db_sort;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Class Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=public)
        % =================================================================
        %> @brief Class constructor
        %>
        %> Initializes the finite element beam
        %>
        %> @param nel number of beam elements
        %> @param crosssection string with type of crosssection of beam
        %> @param varargin string containing beam name
        %>
        %> @return instance of the class_beam
        % =================================================================
         function obj = class_beam(nel,crosssection,varargin)
             %set identifier if passed to function
             if ~isempty(varargin)
                obj.identifier=cell2mat(varargin{1});
             end
             
             %initalize variables
             obj.nel=nel;
             obj.ndof=obj.el_ndof*(obj.nel+1);
             
             obj.K=zeros(obj.ndof);
             obj.M=zeros(obj.ndof);
             obj.M_lumped=zeros(obj.ndof);
             obj.Q=zeros(obj.ndof,1);
             obj.P=zeros(obj.ndof,1);
                          
             
             obj.nodal_deflections=zeros(nel*obj.el_ndof+obj.el_ndof,1);
             obj.nodal_deflections_loc=zeros(nel*obj.el_ndof+obj.el_ndof,1);

             obj.node_loadings=zeros(nel*obj.el_ndof+obj.el_ndof,1);
             obj.nodal_loads=zeros(nel*obj.el_ndof+obj.el_ndof,1);
             obj.nodal_loads_def=zeros(nel*obj.el_ndof+obj.el_ndof,1);
             obj.nodal_masses=zeros((nel+1),7);
             % only required for non-linear elements... possibly
             % differentiate later

             % unpretty and performance killing workaround due to lack of object oriented
             % functionality of MATALB, or my incapability
                for i=1:nel
                    if strcmp(crosssection,'none')
                        beamelement(i)=class_beamelement(crosssection,obj);
                    else
                        beamelement(i)=class_beamelement(eval(crosssection),obj);
                    end
                    beamelement(i).nodal_deflections_loc=zeros(2*obj.el_ndof,1); 
                end
             obj.beamelement=beamelement;
         end         
        % =================================================================
        %> @brief add a boundary condition to the beam
        %>
        %> Adds a structural boundary condition to the beam
        %>
        %> @param bc boundary condition of type class_boundary_condition
        %>
        %> @return instance of the class_beam
        % =================================================================
         function obj=f_add_boundary_condition(obj,bc)
             obj.n_bc=obj.n_bc+1;
             if obj.n_bc>1
                obj.boundary_condition(obj.n_bc)=bc;
             else
                obj.boundary_condition=bc; 
             end
         end
        % =================================================================
        %> @brief set the calculation state of the beam model
        %>
        %> Sets the state of the beam used for computation
        %>
        %> @param state state of type class_aircraft_state
        %>
        %> @return instance of the class_beam
        % =================================================================
         function obj=f_set_state(obj,state)
             obj.load_factor=state.load_factor;
             obj.loadcase_index=state.loadcase_index;
         end
        % =================================================================
        %> @brief reset the calculation state to standard values
        %> (loadfactor=1)
        %>
        %> Resets the state of the beam back to standard values
        %>
        %> @return instance of the class_beam
        % =================================================================
         function obj=f_reset_state(obj)
            obj.load_factor=1;
            obj.loadcase_index=0;
         end
        % =================================================================
        %> @brief makes a deep copy of the class
        %>
        %> Performs a deep copy of the class
        %>
        %> @return deep copy of the class_beam
        % =================================================================
         function new = f_copy(this)
            % Instantiate new object of the same class.
             new = eval([class(this) '(' num2str(this.nel) ',''' class(this.beamelement(1).crosssection) ''')']);

            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
            new.(p{i}) = this.(p{i});
            end
         end
        % =================================================================
        %> @brief returns the name of the beam as a string
        %>
        %> @return name of beam
        % =================================================================
         function name=f_get_name(obj)
             if ~isempty(obj.identifier)
                name=obj.identifier;
             else
                 name='identifier not set';
             end
         end
         
         
         
         function obj=f_update_node_coords(obj)
             obj.node_coords(:,1)=obj.node_coords(:,1)+obj.r_Ref(1);
             obj.node_coords(:,2)=obj.node_coords(:,2)+obj.r_Ref(2);
             obj.node_coords(:,3)=obj.node_coords(:,3)+obj.r_Ref(3);
         end
         
        % =================================================================
        %> @brief calls the solver
        %>
        %> @param eval_nonlin flag for linear or nonlinear evaluation
        %>  (0=linear, 1=nonlinear)
        %> @param add_eigenmass consider the eigenmass of the beam for the
        %>  analysis
        %> @param add_fuelmass consider the fuelmass of the beam for the
        %>  analysis
        %> @param add_engineforces consider the engineforces of the beam for the
        %>  analysis
        %> @param add_gearforces consider the gearforces of the beam for the
        %>  analysis
        %>
        %> @return solved instance of the class_beam
        % =================================================================
         function obj =f_solve(obj) 
             
             eval_nonlin=obj.settings.nonlinear;
             add_eigenmass=obj.settings.gravity;
             add_fuelmass=obj.settings.fuel_mass;
             add_engineforces=obj.settings.engine;
             add_gearforces=obj.settings.landing_gear;
             % f_assemble(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces)            
             if(eval_nonlin)
                disp(['               solving structures (nonlinear) for: ' obj.identifier]); 
                %% assemble system matrices and vectors
                obj.err=100;
                ctr=1;
                %% reset deflections before calculation
                obj=obj.f_reset();
              %   ls=1;
                for ls=0.2:0.01:1
               % disp(['                 nonlinear loadstep: ' num2str(ls)]);
                    obj.err=100;
                    while obj.err> 2e-1
                        obj=obj.f_nonlinassemble(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces,ls); 
                        %% apply boundary conditions
                        obj=obj.f_set_BC();
                        %% call nonlinear system solver
                        obj=obj.f_nonlinsolve();
                        %% postprocess results
                        obj=obj.f_postprocess_nonlin();
                    
                        disp(['                    nonlinear iteration: ' num2str(ctr) '    error: ' num2str(obj.err)]);
                        %% increase counter
                        ctr=ctr+1;
                    end
               end
             else
                disp(['               solving structures (linear) for: ' obj.identifier]); 
                %% assemble system matrices and vectors
                obj=obj.f_assemble(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces); 
                %% set boundary conditions
                %obj=obj.f_applyBC();
                obj=obj.f_set_BC();
                %% call linear system solver
                obj=obj.f_linsolve();
                %% postprocess results
                obj=obj.f_postprocess();
             end
         end
        % =================================================================
        %> @brief calls the eigenvalue and eigenmode solver
        %>
        %>
        %> @return solved instance of the class_beam
        % =================================================================         
         function obj=f_solve_modes(obj)
            K=obj.Kff;
            M=obj.Mff;
            A=K\M; % inverse iteration method (http://en.wikipedia.org/wiki/Modal_analysis_using_FEM )
            [obj.modeshapes, omega2]= eig(A);
            
            j=1;
            for i = 1:size(obj.Mff,2)
                [ EMA(j), idx(j)] = max(abs(obj.modeshapes(:,i))); % returns the maximum of each column vector from the V matrix
                r(j) = mod(idx(j),6); % output 0 means the 6th DOF
                j=j+1;
            end
            obj.modefrequencies = (1./sqrt(diag(omega2)))/(2*pi); % Since inverse iteration is used in determining matrix A

         end
        function obj=f_solve_modes_simon(obj,side)
             %linker fl�gel constrained an wurzel
             if strcmp(side,'left')
             K=obj.Kff(1:size(obj.beamelement,2)*6/2,1:size(obj.beamelement,2)*6/2);
             M=obj.Mff(1:size(obj.beamelement,2)*6/2,1:size(obj.beamelement,2)*6/2);
             %rechter fl�gel constrained an wurzel
             elseif strcmp(side,'right')
             K=obj.Kff(size(obj.beamelement,2)*6/2+7:end,size(obj.beamelement,2)*6/2+7:end);
             M=obj.Mff(size(obj.beamelement,2)*6/2+7:end,size(obj.beamelement,2)*6/2+7:end);
             %beide seiten frei frei
             else
             K=obj.Kff;
             M=obj.Mff;
             end
            A=K\M; % inverse iteration method (http://en.wikipedia.org/wiki/Modal_analysis_using_FEM )
            [obj.modeshapes, omega2]= eig(A);
            
            j=1;
%             for i = 1:size(obj.Mff,2)
%                 [ EMA(j), idx(j)] = max(abs(obj.modeshapes(:,i))); % returns the maximum of each column vector from the V matrix
%                 r(j) = mod(idx(j),6); % output 0 means the 6th DOF
%                 j=j+1;
%             end
            obj.modefrequencies = (1./sqrt(diag(omega2)))/(2*pi); % Since inverse iteration is used in determining matrix A

         end
        % =================================================================
        %> @brief performs self structural design based on the input loads
        %>
        %> @param eval_nonlin flag for linear or nonlinear evaluation
        %>  (0=linear, 1=nonlinear)
        %> @param add_eigenmass consider the eigenmass of the beam for the
        %>  analysis
        %> @param add_fuelmass consider the fuelmass of the beam for the
        %>  analysis
        %> @param add_engineforces consider the engineforces of the beam for the
        %>  analysis
        %> @param add_gearforces consider the gearforces of the beam for the
        %>  analysis
        %> @param weights weights structure
        %>  analysis
        %>
        %> @return solved instance of the class_beam
        % =================================================================
         function obj=f_load_based_self_design(obj,weights,overwrite)
             
            % solve for initial internal forces, always with linear
            % theory!!!
            disp('          performing load based self design') 
            disp('            self design initial calculation');
            mem=obj.settings.nonlinear;
            obj.settings.nonlinear=0;
            obj=obj.f_solve();
            obj.settings.nonlinear=mem;
            
            %convergence criterion is the bending moment at the wingroot,
            %think of something better later...
            prev_val=obj.node_loadings(4);
            err=1000;
            its=1;
            %% Structural Layout Iteration
            while(err>5)
                %% perform structural layout (determine wing box size from loading)
                % first perform initial layout
                obj=obj.f_structural_layout(overwrite);
                % calculate mass
                obj=obj.f_calc_mass(weights);
                % resolve system with new mass
                obj=obj.f_solve();
                
                err=100*abs((prev_val-obj.m_total)/abs(obj.m_total));
                
                disp(['            self design iteration: ' num2str(its) '      error: ' num2str(err) '%']);
              
                prev_val=obj.m_total;
                its=its+1;
            end
         end
         
         %%
         % Calculates and returns the moment of inertia matrix.
         % Uses the mass matrix obj.Mff_lumped and needs the reference point as input.
         % All calculations are done in the body fixed system.
         function  Inertia = f_compute_moment_of_inertia(obj, referencePoint)
            k=1;
            for j=1:1:size((obj.node_coords),1)
                s_inertial(k:k+5)=[obj.node_coords(j,1:3)' [0 0 0]'];
                k=k+6;
            end
            n_dof=length(s_inertial);
            I_Hat=[];
            for i=1:n_dof/6
                I_Hat=[I_Hat; eye(3,3);zeros(3,3)];
            end
            s_inertial=(s_inertial'-I_Hat*referencePoint')';
            k=1;
            for j=1:1:size((obj.node_coords),1)
                delta_inertial(k:k+2)=obj.nodal_deflections((j-1)*6+1:(j-1)*6+3);
                delta_inertial(k+3:k+5)=obj.nodal_deflections((j-1)*6+4:j*6);
                k=k+6;
            end
            b=s_inertial'+delta_inertial';
            b_Hat_Skew=zeros(n_dof,3);
            for bi=1:6:length(b)-1
                b_Hat_Skew(bi:bi+2,1:3)=[0 b(bi+2) -b(bi+1); -b(bi+2) 0 b(bi); b(bi+1) -b(bi) 0];       % for negative cross product
            end
            for bi=4:6:length(b)-1
                b_Hat_Skew(bi:bi+2,1:3)=eye(3,3);   % 3x3 identity matrix
            end
            Inertia=b_Hat_Skew'*obj.Mff_lumped*b_Hat_Skew;
         
         end
         
        %%
        % Calculates and returns the total mass of the beam (only the
        % masses that are in the mass matrix are considered).
        % Uses the mass matrix Mff_lumped
        function mTotal = f_compute_totalMass(obj)
            % delete all rotational DOFs in the mass matrix
            MReduced = obj.Mff_lumped;
            for i=1:length(obj.Mff_lumped)/6
                j = length(obj.Mff_lumped)/6+1-i;
                MReduced(j*6-2:j*6,:) = [];
                MReduced(:,j*6-2:j*6) = [];
            end
            % sum of all masses
            mTotal = sum(sum(MReduced))/3;
        end
         
        function obj=f_set_solver_settings(obj,solver_settings)
            obj.settings=solver_settings;

        end

         obj=f_init_stdgeometry(obj,le,w,h,nu,phi,t_sk_up,t_sk_lo,t_sp_fr,t_sp_re);
        
         obj = f_structural_layout(obj,overwrite);
         
         obj = f_assemble(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces); 
         obj=f_init_stdLoad(obj,qx,qy,qz,dqx,dqy,dqz,mt);
         
         
         obj = f_init_std_Ixyz(obj,nu,phi,Ix,Iy,Iz,J,A,le);
         
         obj = f_init_stdBeam(obj,E,G,m);

         obj = f_calc_mass(obj,weights );
        

         % =================================================================
         %> @brief plot external forces acting on beam
         %>
         %> @param add_eigenmass plot loads due to eigenmass additionally
         %> @param add_fuelmass plot loads due to fuelmass additionally
         %> @param engine plot loads due to engine additionally
         %> @param gear plot loads due to gear additionally
         % ================================================================
          function plot_externalforces(wing,add_eigenmass,add_fuelmass,engine,gear)
            hold on
           
            midp=zeros(3,1);
            
            x_dist=zeros(length(wing.beamelement)+1,1);
            y_dist=zeros(length(wing.beamelement)+1,1);
            z_dist=zeros(length(wing.beamelement)+1,1);
            
            for i=1:1:length(wing.beamelement)
                for j=1:3
                    midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
                end
                h=wing.beamelement(i).crosssection.h;
                w=wing.beamelement(i).crosssection.w;
                
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                twist=wing.beamelement(i).epsilon;
                z_dist(i)=h/2;
                x_dist(i)=w/2;
                z_dist(i+1)=h/2;
                x_dist(i+1)=w/2;
                plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu,twist,0],[1 1 1 1 1 1 1 1],0.3,1);
            end
            
            axis equal
            grid on
            
            qx=cell2mat({wing.beamelement(:).qx});
            qy=cell2mat({wing.beamelement(:).qy});
            qz=cell2mat({wing.beamelement(:).qz});
            halfspan=sum(cell2mat({wing.beamelement(:).le}));
            
            max_q=max([abs(qx),abs(qy),abs(qz)]);
            
            scale=2*max_q/halfspan;
            
            qmrot=zeros(3,length(wing.beamelement));
            qfrot=zeros(3,length(wing.beamelement));
            
            for i=1:length(wing.beamelement)    
                if add_eigenmass
                    T=wing.beamelement(i).f_rotVec(0,wing.beamelement(i).gamma,0);  % transformations matrix from NED to aerodynamic system
                    qmrot(1:3,i)=T*[0,0,wing.beamelement(i).qm]';%wing.beamelement(i).qm];           % distributed loading due to eigenmass in x,y,z coordinates
                end
                
                if add_fuelmass
                    T=wing.beamelement(i).f_rotVec(0,wing.beamelement(i).gamma,0);  % transformations matrix from NED to aerodynamic system
                    qfrot(:,i)=T*[0;0;wing.beamelement(i).qf];           % distributed loading due to eigenmass in x,y,z coordinates
                end
            end
            
            %plotforce(handle,qx,scale,wing.node_coords,'z',0.5);
            %plotforce(handle,qy,scale,wing.node_coords,'z',0.5);
%             for i=1:1:length(wing.node_coords)
%                wing.node_coords(i,:)=wing.node_coords(i,:)+wing.r_Ref(:)';
%             end
    
            xxx=subplot(1,1,1);
            %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_co
            %ords(:,3),'-k','LineWidth',2)
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

            %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_coords(:,3),'-x','LineWidth',2)
            %plot3(engine.cgpos(1),engine.cgpos(2),engine.cgpos(3),'c+');
            %plot3(gear.pos(1),gear.pos(2),gear.pos(3),'go');
         end
         
         % ================================================================
         %> @brief plot geometry of beam
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
         function plot_geometry(wing,varargin)
            hold on
            
               col=1;
               tr=1;
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
                w=wing.beamelement(i).crosssection.w;
                le=wing.beamelement(i).le;
                nu=wing.beamelement(i).nu;
                twist=wing.beamelement(i).epsilon;
                sweep=wing.beamelement(i).phi;
                plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu,twist,sweep],[col col col col col col col col],tr,1);
            end
            axis equal
            grid on
         end  
         
         % ================================================================
         %> @brief plot resulting structure including thicknesses
         %>
         % ================================================================
%          function plot_structure(wing,varargin)
%             mycolormap=jet(256);
%             hold on
%             midp=zeros(3,1);
%             
% 
%                 for i=1:1:length(wing.beamelement) 
%                     t_sk_up(i)=cell2mat({wing.beamelement(i).wingsection.t_sk_up});
%                     t_sk_lo(i)=cell2mat({wing.beamelement(i).wingsection.t_sk_lo});
%                     t_sp_fr(i)=cell2mat({wing.beamelement(i).wingsection.t_sp_fr});
%                     t_sp_re(i)=cell2mat({wing.beamelement(i).wingsection.t_sp_re});
%                 end
%             
%             optargin = size(varargin,2);
%             if optargin==2
%                 t_min=varargin{2};
%                 t_max=varargin{1};
%             else
%                  t_max=0;
%                 t_min=100;
%                 tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re];
%                 t_max=max([tk t_max]);
%                 t_min=min([tk t_min]);
%             end
%             
%             for i=1:1:length(wing.beamelement)
%                 
%                 for j=1:3
%                     midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
%                 end
%                 
%                 h=wing.beamelement(i).wingsection.h;
%                 w=wing.beamelement(i).wingsection.w;
%                 le=wing.beamelement(i).le;
%                 nu=wing.beamelement(i).nu;
%                 twist=wing.beamelement(i).epsilon;
% 
%                 col_i=min(round((t_sk_lo(i)-t_min)*255/(t_max-t_min)+1),255);
%                 colidx(i,1)=col_i;
%                 col=ones(8,1)*mycolormap(col_i,:);
%                 
%                a=-wing.beamelement(i).nu;
%                 b=-wing.beamelement(i).epsilon;
%                 %lower skin
%                 plotcube([midp(1),midp(2)-h/2*sin(a),midp(3)-h/2*cos(a)],[w,le,t_sk_lo(i)],[twist,0,nu],ones(8,1)*col_i/255,1,1);
%                 
%                 
%                 % front spar
%                 col_i=min(round((t_sp_fr(i)-t_min)*255/(t_max-t_min)+1),255);
%                 % colidx(i,2)=col_i;
%                 col=ones(8,1)*mycolormap(col_i,:);
%                 
%                 % rot=wing.nodal_deflections((i-1)*wing.el_ndof+4:(i-1)*wing.el_ndof+6);
%                 
% 
%                 c=0;
%               
%                 Lx=[1       0       0
%                 0   cos(a)  sin(a)
%                 0   -sin(a) cos(a)];
% 
%                 Ly=[cos(b) 0 -sin(b)
%                 0      1    0
%                 sin(b)  0   cos(b)];
% 
% 
%                 M=Ly*Lx;
%                 
%                 diff=[-w/2 0 0];
%                 
%                 xxx=M*diff';
% 
%                 plotcube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_fr(i),le,h],[twist,0,nu],ones(8,1)*col_i/255,1,1);
%                 %rear spar
%                 
%                 diff=[w/2 0 0];
%                 xxx=M*diff';
%                 col_i=min(round((t_sp_re(i)-t_min)*255/(t_max-t_min)+1),255);
%                 plotcube([midp(1)+xxx(1),midp(2)+xxx(2),midp(3)+xxx(3)],[t_sp_re(i),le,h],[twist,0,nu],ones(8,1)*col_i/255,1,1);
%                 %upper skin
%                 col_i=min(round((t_sk_up(i)-t_min)*255/(t_max-t_min)+1),255);
%                 %col=ones(8,1)*mycolormap(col_i,:);
%                 plotcube([midp(1),midp(2)+h/2*sin(a),midp(3)+h/2*cos(a)],[w,le,t_sk_up(i)],[twist,0,nu],ones(8,1)*col_i/255,1,1);
%             end
%             axis equal
%             grid on
%             
%          end
         
         % ================================================================
         %> @brief plot deformations of structure
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

                mid_str=0;%norm_def(2+6*(i-1))-norm_def(6+2+6*(i-1));
                
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
         %> @brief plot internal forces
         %>
         % ================================================================
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
            
%             aQx=wing.node_loadings(1:6:end);
%             aQy=wing.node_loadings(2:6:end);
%             aQz=wing.node_loadings(3:6:end);
%             aMx=wing.node_loadings(4:6:end);
%             aMy=wing.node_loadings(5:6:end);
%             aMz=wing.node_loadings(6:6:end);
%              
            optargin = size(varargin,2);
            if optargin==9
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
            end
 
            % plot Mx 
            handle1=subplot(3,2,4);
            plotforce(handle1,aMx,scale_bendingmoment,wing.node_coords,'z',1,'r',1);
            title('Mx')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            axis equal
            
            % plot Qz 
            handle2=subplot(3,2,1);
            plotforce(handle2,aQz,scale_force,wing.node_coords,'z',1,'r',1);
            title('Qz')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            axis equal
            
            handle4=subplot(3,2,2);
            plotforce(handle4,aQy,scale_force,wing.node_coords,'z',1,'r',1);
            title('Qy')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            axis equal
                
            handle5=subplot(3,2,5);
            plotforce(handle5,aMy,scale_torsionmoment,wing.node_coords,'z',1,'r',1);
            title('My')
            view(-90,0)
            ylim([med_y-0.6*delta,med_y+0.6*delta]);
            zlim([med_z-0.6*delta,med_z+0.6*delta]);
            axis equal
            
            handle3=subplot(3,2,3);
            plotforce(handle3,aMz,scale_bendingmoment,wing.node_coords,'x',1,'r',1);
            title('Mz')
            view(-90,90)
            ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
            xlim([med_x-0.6*deltaxy,med_x+0.6*deltaxy]);
            axis equal
            
            handle6=subplot(3,2,6);
            plotforce(handle6,aQx,scale_force,wing.node_coords,'x',1,'r',1);
            title('Qx')
            view(-90,90)
            ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
            xlim([med_x-0.6*deltaxy,med_x+0.6*deltaxy]);
            axis equal
         end
         
         % ================================================================
         %> @brief compute rotation matrix for rotation of 3D vector
         %>
         %> @param a rotation angle (rad) about x-axis
         %> @param a rotation angle (rad) about y-axis
         %> @param a rotation angle (rad) about z-axis
         % ================================================================
         function T=f_rot_vec3(obj,a,b,c)
                   Lx=[1       0       0
            0   cos(a)  sin(a)
            0   -sin(a) cos(a)];

            Ly=[cos(b) 0 -sin(b)
                0      1    0
               sin(b)  0   cos(b)];

            Lz=[cos(c) sin(c)   0
            -sin(c) cos(c)  0
            0           0   1];

            T=Lz*Ly*Lx;

         end
         function obj=f_reset(obj)
                obj.wf=obj.wf*0;
                obj.nodal_deflections=obj.nodal_deflections.*0;
                for i=1:obj.nel
                    obj.beamelement(i).nodal_deflections_loc=zeros(2*obj.el_ndof,1); 
                end  
         end
         
    end  
    methods (Access=private) 
        % =================================================================
        %> @brief reset deflections (function only required for non-linear
        %>  beam)
        %>
        %> Adds a structural boundary condition to the beam
        %>
        %> @param bc boundary condition of type class_boundary_condition
        %>
        %> @return instance of the class_beam
        % =================================================================
   

     end
end

