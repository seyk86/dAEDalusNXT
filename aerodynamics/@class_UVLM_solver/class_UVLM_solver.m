%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_UVLM_solver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> name of the computation case
        case_name;
        %> aerodynamic state used for solution (type class_aerodynamic_state)
        state;
        %> settings for UVLM computation
        settings;
        %> reference values
        reference;
        %> grid point coordinates
        grid;
        %> panel edge indices of grid
        panels;
        %> panel area
        area;
        %> is the current panel a trailing edge panel?
        is_te;
        %> grid
        grid_vring;
        %> grid point coordinates for wake
        grid_wake;
        
        grid_wake_old;
        
        %> required for postprocessing
        gridpoint_to_panel;
        %> panel edge incices of wake grid
        panels_wake;
        
        t_step;
        n_step;
        %> array containing the force vector application point
        fvap;
        %> collocation point coordinates
        colloc;
        %> normal vector on collocation point
        colloc_nvec;
        
        %% Linearized UVLM
        %> grid point coordinates at reference condition
        grid_init;
        %> collocation point coordinates at reference condition
        colloc_init;
        %> normal vector on collocation points at reference condition
        colloc_nvec_init;
        %> derivatives at the reference condition of the normal vector on collocation points with respect to grid point coordinates
        d_colloc_nvec;
        
        %> aerodynamic influence coefficient matrix wing on wing
        Abb;
        
        Abb_ff;
        %> aerodynamic influence coefficient matrix wake on wing
        Abw;
        
        Abw_ff;
        %> aerodynamic influence coefficient matrix wake on wake
        Aww;
        
        Aww_ff;
        %> aerodynamic influence coefficient matrix wing on wake
        Awb;
        
        Awb_ff;
        Abb_x;
        Aww_x;
        Awb_x;
        Abb_y;
        Aww_y;
        Awb_y;
        Abb_z;
        Aww_z;
        Awb_z;
        
        Abw_x;
        Abw_y;
        Abw_z;
        
        %> time propagation matrix for wing
        Cbb;
        %> time propagation matrix for wake
        Cbw;
        
        %> no penetration boundary condition vector
        b;
        
        b_check;
        %> vorticity vector
        Gamma;
        %> vorticity of previous timestep
        Gamma_prv;
        
        %> vorticity of the wake panels
        Gamma_wake;
        alpha_s;
        F_body;
        % Moment on every panel due to unsteady force application at midpoint
        M_body;
        
        F_body_unsteady;
        F_body2;
        F_aero;
        F_drag;
        
        t_vec;
        row_length;
        % constant free stream velocity over entire aircraft
        Uinf;
        
        % use if varying free stream velocity over aircraft x-axis
        Uinf_x;
        % respective coordinate of Uinf_x
        x_kin;
        
        cl;
        cdi;
        cp;
        cp_st;
        cp_un;
        cz;
        cy;
        cx;
        cl_n;
        rho;
        qinf;
        r;
        r_mid;
        r_dwn;
        dvortex;
        
        Ma_corr;
        
        gust_shape;
        gust_flag;
        gust_t_start;
        uind;
        
        C_modes;
        %%%% RESULTS
        %% aerodynamic coefficients
        %> Coefficients in body axis
        CX;
        CY;
        CZ;
        Cl;
        Cy;
        Cdi;
        Cdi2;
        CL;
        CM;
        CN;
        %> Alpha Derivatives in body axis
        CXa;
        CYa;
        CZa;
        Cla;
        Cya;
        Cdia;
        CLa;
        CMa;
        CNa;
        %> Beta Derivatives in body axis
        CXb;
        CYb;
        CZb;
        Clb;
        Cyb;
        Cdib;
        CLb;
        CMb;
        CNb;
        %% damping derivatives
        % pitch
        CXq;
        CYq;
        CZq;
        Clq;
        Cyq;
        Cdiq;
        CLq;
        CMq;
        CNq;
        % roll
        CXp;
        CYp;
        CZp;
        Clp;
        Cyp;
        Cdip;
        CLp;
        CMp;
        CNp;
        % yaw
        CXr;
        CYr;
        CZr;
        Clr;
        Cyr;
        Cdir;
        CLr;
        CMr;
        CNr;
        
        
        CX_complex;
        CZ_complex;
        Cl_complex;
        CM_complex;
        CY_complex;
        CL_complex;
        CN_complex;
        
        C_MODES_complex;
        
        Q_modes;
        Q_modes_complex;
        
        
        flag;
        
        crossp;
        grid_wake_init;
        panels_wake_init;
        panel_len;
        alpha_vec;
        
        u_test;
        acc;
        vel;
        pos;
        
        v_ind;
        %state space model valid for one velocity
        sys;
        %lpv state space model valid for all speeds
        lpvSys;
        max_gust_loads;
    end
    
    methods
        
        obj=solve_unsteady_aeroelastic_gust_response(obj,aircraft,aircraft_structure,Uds,H,t_start,t_end,init_def);
        obj=solve_unsteady_aeroelastic_gust_response_linear(obj,aircraft,aircraft_structure,Uds,H,t_start,t_end,init_def,flag_wake);
        obj=solve_unsteady_aeroelastic_gust_response_rom(obj,aircraft,aircraft_structure,Uds,H,t_start,t_end,init_def,n_modes,Cdi);
        obj=solve_free_flying_aeroelastic_gust_response(obj,aircraft,aircraft_structure,ac_state,Uds,H,t_start,t_end,init_def,Cdi);
        obj=solve_free_flying_rigid_gust_response(obj,ref_state,Uds,H,t_start,t_end,aircraft,Cdi);
        obj=solve_unsteady_aeroelastic_theta_response(obj,aircraft,aircraft_structure,Uds,H,t_start,t_end,init_def);
        obj=solve_free_flying_aeroelastic_simulation_fixed_axis(obj,aircraft,aircraft_structure,ref_state,Uds,H,t_start,t_end,init_def,Cdi,control_signals);
        obj=solve_free_flying_aeroelastic_simulation_ATTAS(obj,aircraft,aircraft_structure,ref_state,Uds,H,t_start,t_end,init_def,Cdi,control_signals);
        obj=initialize_time_domain_solution(obj,t_step);
        obj=solve_time_domain_aerodynamics(obj,aircraft,x_body,V,alpha,beta,pqr,rho_air,itr,approach);

        function obj=class_UVLM_solver(name,grid,is_te,panels,state,grid_wake,panels_wake,reference,settings)
            
            obj.case_name=name;
            
            obj.settings=settings;
            % initialization procedure
            obj.state=state;
            obj.rho=state.rho_air;
            obj.Ma_corr=sqrt(1-state.Ma^2);
            obj.Ma_corr=1;
            
            %% TODO: Ma correction
            % obj.Ma_corr=1;%sqrt(1-state.Ma^2);
            obj.state.p_ref=reference.p_ref;
            obj.grid=grid;
            obj.panels=panels;
            obj.is_te=is_te;
            obj.reference=reference;
            obj.r=zeros(3,length(obj.panels));
            
            %obj.grid_wake_init=grid_wake;
            obj.panels_wake_init=panels_wake;
            
            % compute collocation point coordinates
            for i=1:length(obj.panels)
                %from p_ref to 1/4
                obj.r(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.75+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.25-obj.state.p_ref';
                %from p_ref to 3/4
                obj.r_dwn(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.25+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.75-obj.state.p_ref';
                %from p_ref to 1/2
                obj.r_mid(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.5+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.5-obj.state.p_ref';
            end
            
            if obj.Ma_corr<1
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,obj.Ma_corr);
                obj.grid(1,:)=grid(1,:)/obj.Ma_corr;
                obj.grid(2,:)=grid(2,:);
                obj.grid(3,:)=grid(3,:);
                
                obj.grid_wake(1,:)=grid_wake(1,:)/obj.Ma_corr;
                obj.grid_wake(2,:)=grid_wake(2,:);
                obj.grid_wake(3,:)=grid_wake(3,:);
            else
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,1);
            end
            
            obj.gridpoint_to_panel=zeros(length(obj.grid),4);
            for i=1:length(obj.panels)
                for j=1:4
                    idx=obj.panels(j,i);
                    k=1;
                    while obj.gridpoint_to_panel(idx,k)~=0
                        k=k+1;
                    end
                    if k<=4
                        obj.gridpoint_to_panel(idx,k)=i;
                    end
                end
            end
            
            obj=obj.initialize_vring_grid(0.002);
            obj=obj.initialize_wake_grid(10,0.002);
            
            obj.qinf=1/2*state.rho_air*norm(obj.Uinf)^2;
            obj=obj.compute_colloc_points();
            obj=obj.compute_panel_area();
        end
        
        function obj=initialize_vring_grid(obj,t_step)
            grid_wake=[];
            
            wake_start=0.25;
            % compute actual grid (vring grid)
            is_done=zeros(1,length(obj.grid));
            for i=1:1:length(obj.panels)
                if is_done(obj.panels(1,i))==0
                    obj.grid_vring(:,obj.panels(1,i))=0.75*obj.grid(:,obj.panels(1,i))+0.25*obj.grid(:,obj.panels(4,i));
                    is_done(obj.panels(1,i))=1;
                end
                if is_done(obj.panels(2,i))==0
                    obj.grid_vring(:,obj.panels(2,i))=0.75*obj.grid(:,obj.panels(2,i))+0.25*obj.grid(:,obj.panels(3,i));
                    is_done(obj.panels(2,i))=1;
                end
                if is_done(obj.panels(4,i))==0
                    if obj.is_te(i)==0
                        obj.grid_vring(:,obj.panels(4,i))=obj.grid(:,obj.panels(4,i))+0.25*(obj.grid(:,obj.panels(4,i+1))-obj.grid(:,obj.panels(1,i+1)));
                        is_done(obj.panels(4,i))=1;
                    else
                        %obj.grid_vring(:,obj.panels(4,i))=obj.grid(:,obj.panels(4,i))+wake_start*norm(obj.Uinf)*t_step*(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i)))/norm((obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i))));
                        obj.grid_vring(:,obj.panels(4,i))=obj.grid(:,obj.panels(4,i))+wake_start*(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i)));
                        is_done(obj.panels(4,i))=1;
                        %grid_wake=[grid_wake obj.grid_vring(:,obj.panels(4,i))];
                    end
                end
                if is_done(obj.panels(3,i))==0
                    if obj.is_te(i)==0
                        obj.grid_vring(:,obj.panels(3,i))=obj.grid(:,obj.panels(3,i))+0.25*(obj.grid(:,obj.panels(3,i+1))-obj.grid(:,obj.panels(2,i+1)));
                        is_done(obj.panels(3,i))=1;
                    else
                        %obj.grid_vring(:,obj.panels(3,i))=obj.grid(:,obj.panels(3,i))+wake_start*norm(obj.Uinf)*t_step*(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i)))/norm((obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i))));
                        obj.grid_vring(:,obj.panels(3,i))=obj.grid(:,obj.panels(3,i))+wake_start*(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i)));
                        is_done(obj.panels(3,i))=1;
                        %grid_wake=[grid_wake obj.grid_vring(:,obj.panels(3,i))];
                    end
                end
            end
        end
        
        function obj=initialize_wake_grid(obj,n_step,t_step)
            
            %grid_wake=obj.grid_wake_init;
            obj.grid_wake=[];%obj.grid_wake_init;
            panels_wake=obj.panels_wake_init;
            obj.panels_wake=[];%obj.panels_wake_init;
            
            % set wake
            obj.n_step=n_step;
            obj.t_step=t_step;
            
            idx_1_prv=1;
            idx_2_prv=0;
            grid_te_line=[];
            prv_jump_idx=0;
            for i=1:1:length(obj.panels)
                if obj.is_te(i)
                    if i==1
                        idx_1_prv=obj.panels(3,i);
                        idx_2_prv=obj.panels(4,i);
                        if obj.panels(3,i)>obj.panels(4,i)
                            grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(4,i)) obj.grid_vring(:,obj.panels(3,i))];
                        else
                            grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(3,i)) obj.grid_vring(:,obj.panels(4,i))];
                        end
                    else
                        if obj.panels(3,i)>obj.panels(4,i)
                            if obj.panels(4,i)==idx_1_prv
                                grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(3,i))];
                            else
                                if prv_jump_idx==0
                                    grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(4,i))  obj.grid_vring(:,obj.panels(3,i))];
                                else
                                    grid_te_line=[grid_te_line  grid_te_line(:,prv_jump_idx:end)  obj.grid_vring(:,obj.panels(4,i))  obj.grid_vring(:,obj.panels(3,i))];
                                end
                                prv_jump_idx=size(grid_te_line,2)-1;
                            end
                        else
                            if obj.panels(3,i)==idx_2_prv
                                grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(4,i))];
                            else
                                if prv_jump_idx==0
                                    grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(3,i)) obj.grid_vring(:,obj.panels(4,i))];
                                else
                                    grid_te_line=[grid_te_line grid_te_line(:,prv_jump_idx:end)  obj.grid_vring(:,obj.panels(3,i)) obj.grid_vring(:,obj.panels(4,i))];
                                end
                                prv_jump_idx=size(grid_te_line,2)-1;
                            end
                            
                        end
                        idx_1_prv=obj.panels(3,i);
                        idx_2_prv=obj.panels(4,i);
                    end
                end
            end
            grid_te_line=[grid_te_line grid_te_line(:,prv_jump_idx:end)];
            
            obj.grid_wake=grid_te_line;
            
            obj.panels_wake=panels_wake;
            
            is_done=zeros(1,length(obj.grid_wake));
            
            for i=1:length(obj.panels_wake)
                if is_done(obj.panels_wake(3,i))==0
                    %obj.grid_wake(:,obj.panels_wake(3,i))=obj.grid_wake(:,obj.panels_wake(3,i))+[obj.Uinf(1)/obj.Ma_corr 0 0]'*obj.t_step*n_step;
                    obj.grid_wake(:,obj.panels_wake(3,i))=obj.grid_wake(:,obj.panels_wake(3,i))+[obj.Uinf(1) 0 0]'*obj.t_step*n_step;
                    is_done(obj.panels_wake(3,i))=1;
                end
                if is_done(obj.panels_wake(4,i))==0
                    %obj.grid_wake(:,obj.panels_wake(4,i))=obj.grid_wake(:,obj.panels_wake(4,i))+[obj.Uinf(1)/obj.Ma_corr 0 0]'*obj.t_step*n_step;
                    obj.grid_wake(:,obj.panels_wake(4,i))=obj.grid_wake(:,obj.panels_wake(4,i))+[obj.Uinf(1) 0 0]'*obj.t_step*n_step;
                    is_done(obj.panels_wake(4,i))=1;
                end
            end
            
            jump_idx_grid=0;
            jump_idx_panel=[];
            
            inv_wake_panel=zeros(length(panels_wake),1);
            
            for i=1:length(panels_wake)
                if panels_wake(3,i)>panels_wake(4,i)
                    inv_wake_panel(i)=0;
                else
                    inv_wake_panel(i)=1;
                end
            end
            
            for i=1:length(panels_wake)-1
                if abs(panels_wake(1,i)-panels_wake(1,i+1))>1.5
                    jump_idx_panel=[jump_idx_panel i];
                    if panels_wake(3,i)>panels_wake(4,i)
                        jump_idx_grid=[jump_idx_grid panels_wake(3,i)];
                    else
                        jump_idx_grid=[jump_idx_grid panels_wake(4,i)];
                    end
                end
            end
            
            jump_idx_grid=[jump_idx_grid length(obj.grid_wake)];
            jump_idx_panel=[jump_idx_panel length(obj.panels_wake)];
            
            grid_wake=obj.grid_wake;
            row_len=size(grid_wake,2)/2;
            obj.row_length=row_len;
            panel_len=size(panels_wake,2);
            
            test_wake=[];
            for i=0:1:n_step
                for j=1:length(jump_idx_grid)-1
                    next_row=grid_wake(:,(jump_idx_grid(j)+1):(jump_idx_grid(j)+(jump_idx_grid(j+1)-jump_idx_grid(j))/2))*(1-i*obj.t_step/(obj.t_step*n_step))+grid_wake(:,(jump_idx_grid(j)+(jump_idx_grid(j+1)-jump_idx_grid(j))/2+1):(jump_idx_grid(j)+(jump_idx_grid(j+1)-jump_idx_grid(j))))*(i*obj.t_step/(obj.t_step*n_step));
                    %obj.grid_wake=[obj.grid_wake next_row];
                    %next_row(1,:)=next_row(1,:)/obj.Ma_corr;
                    test_wake=[test_wake next_row];
                end
            end
            
            %             if obj.Ma_corr<1
            %                 test_wake(1,:)=test_wake(1,:)/obj.Ma_corr;
            %             end
            
            obj.grid_wake=test_wake;
            
            obj.grid_wake_old=test_wake;
            panel_row=obj.panels_wake;
            first_row=obj.panels_wake;
            
            first_row(3:4,:)=first_row(3:4,:)+row_len;
            
            last_row=obj.panels_wake;
            last_row(1:2,:)=first_row(1:2,:)+row_len*(1+n_step);
            
            obj.panels_wake=first_row;
            k=0;
            j=1;
            n=0;
            
            for i=1:n_step*length(panel_row)
                %if j~=length(jump_idx_panel)
                if inv_wake_panel((i-n*(length(panel_row))))==0
                    %if inv_wake_panel((i-n*(length(panel_row))))==0
                    test_panels(1,i)=i+k;
                    test_panels(2,i)=i+1+k;
                    test_panels(3,i)=i+1+length(grid_wake)/2+k;
                    test_panels(4,i)=i+length(grid_wake)/2+k;
                else
                    test_panels(1,i)=i+k+1;
                    test_panels(2,i)=i+k;
                    test_panels(3,i)=i+length(grid_wake)/2+k;
                    test_panels(4,i)=i+length(grid_wake)/2+k+1;
                end
                if (i-n*(length(panel_row)))==jump_idx_panel(j)
                    if j<length(jump_idx_panel)
                        j=j+1;
                        k=k+1;
                    else
                        j=1;
                        k=k+1;
                        n=n+1;
                    end
                end
            end
            
            obj.panels_wake=test_panels;%[obj.panels_wake last_row];
            % initialize time marching matrices
            obj.Cbb=zeros(size(obj.panels_wake,2),size(obj.panels,2));
            k=1;
            for i=1:size(obj.panels,2)
                if obj.is_te(i)==1
                    %len_b=norm(obj.grid(:,obj.panels(2,i))-obj.grid(:,obj.panels(1,i)))+norm(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i)))+norm(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(3,i)))+norm(obj.grid(:,obj.panels(1,i))-obj.grid(:,obj.panels(4,i)));
                    %len_w=norm(obj.grid_wake(:,obj.panels_wake(2,k))-obj.grid_wake(:,obj.panels_wake(1,k)))+norm(obj.grid_wake(:,obj.panels_wake(3,k))-obj.grid_wake(:,obj.panels_wake(2,k)))+norm(obj.grid_wake(:,obj.panels_wake(4,k))-obj.grid_wake(:,obj.panels_wake(3,k)))+norm(obj.grid_wake(:,obj.panels_wake(1,k))-obj.grid_wake(:,obj.panels_wake(4,k)));
                    obj.Cbb(k,i)=1;
                    k=k+1;
                end
            end
%            obj.Cbw=zeros(length(obj.panels_wake),length(obj.panels_wake));
%            for i=1:length(obj.panels_wake)-panel_len
%                obj.Cbw(i+panel_len,i)=1;
%            end
            obj.panel_len=panel_len;
        end
        
        function obj=f_set_state(obj,state)
            obj.Ma_corr=sqrt(1-state.Ma^2);
            obj.Uinf=state.V_inf;
            obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,obj.Ma_corr);
            obj.qinf=1/2*state.rho_air*norm(obj.Uinf)^2;
            obj.rho=state.rho_air;
            obj.state=state;
        end
        
        function obj=set_grid(obj,grid,panels)
            obj.grid(1,:)=grid(1,:)/obj.Ma_corr;
            obj.grid(2,:)=grid(2,:);
            obj.grid(3,:)=grid(3,:);
            obj.panels=panels;
            obj=obj.compute_colloc_points();
        end
        
        function obj=update_p_ref(obj,p_ref)
            obj.state.p_ref=p_ref;
            obj.reference.p_ref=p_ref;
            obj.r=zeros(3,length(obj.panels));
            % careful (PG correction)
            % compute collocation point coordinates
            obj.grid(1,:)=obj.grid(1,:)*obj.Ma_corr;
            for i=1:length(obj.panels)
                obj.r(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.75+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.25-obj.state.p_ref';
                obj.r_dwn(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.25+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.75-obj.state.p_ref';
            end
            obj.grid(1,:)=obj.grid(1,:)/obj.Ma_corr;
        end
        
        function obj=update_grid(obj)
            is_done=zeros(1,length(obj.grid));
            for i=1:1:length(obj.panels)
                if is_done(obj.panels(1,i))==0
                    obj.grid_vring(:,obj.panels(1,i))=0.75*obj.grid(:,obj.panels(1,i))+0.25*obj.grid(:,obj.panels(4,i));
                    is_done(obj.panels(1,i))=1;
                end
                if is_done(obj.panels(2,i))==0
                    obj.grid_vring(:,obj.panels(2,i))=0.75*obj.grid(:,obj.panels(2,i))+0.25*obj.grid(:,obj.panels(3,i));
                    is_done(obj.panels(2,i))=1;
                end
                if is_done(obj.panels(3,i))==0
                    if obj.is_te(i)==0
                        obj.grid_vring(:,obj.panels(3,i))=obj.grid(:,obj.panels(3,i))+0.25*(obj.grid(:,obj.panels(3,i+1))-obj.grid(:,obj.panels(2,i+1)));
                        is_done(obj.panels(3,i))=1;
                    else
                        obj.grid_vring(:,obj.panels(3,i))=obj.grid(:,obj.panels(3,i))+0.25*(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i)));
                        is_done(obj.panels(3,i))=1;
                    end
                end
                if is_done(obj.panels(4,i))==0
                    if obj.is_te(i)==0
                        obj.grid_vring(:,obj.panels(4,i))=obj.grid(:,obj.panels(4,i))+0.25*(obj.grid(:,obj.panels(4,i+1))-obj.grid(:,obj.panels(1,i+1)));
                        is_done(obj.panels(4,i))=1;
                    else
                        obj.grid_vring(:,obj.panels(4,i))=obj.grid(:,obj.panels(4,i))+0.25*(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i)));
                        is_done(obj.panels(4,i))=1;
                    end
                end
            end
            
            obj.grid(1,:)=obj.grid(1,:)*obj.Ma_corr;
            for i=1:length(obj.panels)
                obj.r(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.75+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.25-obj.state.p_ref';
                obj.r_dwn(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.25+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.75-obj.state.p_ref';
            end
            obj.grid(1,:)=obj.grid(1,:)/obj.Ma_corr;
            obj=obj.update_trailing_edge_line();
            obj=obj.compute_colloc_points();
        end
        
        function obj=update_trailing_edge_line(obj)
            idx_1_prv=0;
            idx_2_prv=0;
            grid_te_line=[];
            for i=1:1:length(obj.panels)
                idx_1=obj.panels(3,i);
                idx_2=obj.panels(4,i);
                if obj.is_te(i)
                    if i==1
                        if idx_1>idx_2
                            grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(4,i))  obj.grid_vring(:,obj.panels(3,i))];
                        else
                            grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(4,i))];
                        end
                    else
                        if idx_1>idx_2
                            if obj.panels(4,i)==idx_1_prv
                                grid_te_line=[grid_te_line obj.grid_vring(:,obj.panels(3,i))];
                            else
                                grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(4,i))  obj.grid_vring(:,obj.panels(3,i))];
                            end
                        else
                            if obj.panels(3,i)==idx_2_prv
                                grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(4,i))];
                            else
                                grid_te_line=[grid_te_line  obj.grid_vring(:,obj.panels(3,i)) obj.grid_vring(:,obj.panels(4,i))];
                            end
                            
                        end
                        idx_1_prv=obj.panels(3,i);
                        idx_2_prv=obj.panels(4,i);
                    end
                end
            end
            obj.grid_wake(:,1:obj.row_length)=grid_te_line;
        end
%         
%         function new = f_copy(this)
%             % Instantiate new object of the same class.
%             class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,wingaero.state,aircraft.reference);
%             new = eval([class(this) '(' this.grid ',' this.te_idx ',' this.panels, ',' this.state ',' this.referencen ')']);
%             % Copy all non-hidden properties.
%             p = properties(this);
%             for i = 1:length(p)
%                 new.(p{i}) = this.(p{i});
%             end
%         end
%         
        function obj=compute_colloc_points(obj)
            obj.colloc=0.125*obj.grid(:,obj.panels(1,:))+0.125*obj.grid(:,obj.panels(2,:))+0.375*obj.grid(:,obj.panels(3,:))+0.375*obj.grid(:,obj.panels(4,:));
            obj.fvap=0.375*obj.grid(:,obj.panels(1,:))+0.375*obj.grid(:,obj.panels(2,:))+0.125*obj.grid(:,obj.panels(3,:))+0.125*obj.grid(:,obj.panels(4,:));
            for i=1:1:length(obj.panels)
                % obj.colloc(:,i)=0.25*obj.grid_vring(:,obj.panels(1,i))+0.25*obj.grid_vring(:,obj.panels(2,i))+0.25*obj.grid_vring(:,obj.panels(3,i))+0.25*obj.grid_vring(:,obj.panels(4,i));
                %obj.colloc(:,i)=0.125*obj.grid(:,obj.panels(1,i))+0.125*obj.grid(:,obj.panels(2,i))+0.375*obj.grid(:,obj.panels(3,i))+0.375*obj.grid(:,obj.panels(4,i));
                %obj.fvap(:,i)=0.375*obj.grid(:,obj.panels(1,i))+0.375*obj.grid(:,obj.panels(2,i))+0.125*obj.grid(:,obj.panels(3,i))+0.125*obj.grid(:,obj.panels(4,i));
                %n_vec=cross(obj.grid_vring(:,obj.panels(3,i))-obj.grid_vring(:,obj.panels(1,i)),obj.grid_vring(:,obj.panels(4,i))-obj.grid_vring(:,obj.panels(2,i)));
                %% cross is too slow
                %n_vec=cross(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(1,i)),obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(2,i)));
                vec1=obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(1,i));
                vec2=obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(2,i));
                n_vec(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2);
                n_vec(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3);
                n_vec(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1);
                obj.colloc_nvec(:,i)=n_vec/sqrt(n_vec(1)*n_vec(1)+n_vec(2)*n_vec(2)+n_vec(3)*n_vec(3));
            end
        end

        function obj=compute_d_colloc_nvec(obj)
            
            colloc_nvec_init=zeros(3,length(obj.panels));
            d_colloc_nvec=zeros(3,12,length(obj.panels));
            
            for i=1:1:length(obj.panels)
                
                vec1=obj.grid_init(:,obj.panels(3,i))-obj.grid_init(:,obj.panels(1,i));
                vec2=obj.grid_init(:,obj.panels(4,i))-obj.grid_init(:,obj.panels(2,i));
                
                n_vec(1,1)=vec1(2)*vec2(3)-vec1(3)*vec2(2);
                n_vec(2,1)=vec1(3)*vec2(1)-vec1(1)*vec2(3);
                n_vec(3,1)=vec1(1)*vec2(2)-vec1(2)*vec2(1);
                
                norm_n_vec=sqrt(n_vec(1)*n_vec(1)+n_vec(2)*n_vec(2)+n_vec(3)*n_vec(3));
                
                colloc_nvec_init(:,i)=n_vec/norm_n_vec;
                
                d_nvec=[0 obj.grid_init(3,obj.panels(2,i))-obj.grid_init(3,obj.panels(4,i)) obj.grid_init(2,obj.panels(4,i))-obj.grid_init(2,obj.panels(2,i))...
                    0 obj.grid_init(3,obj.panels(3,i))-obj.grid_init(3,obj.panels(1,i)) obj.grid_init(2,obj.panels(1,i))-obj.grid_init(2,obj.panels(3,i))...
                    0 obj.grid_init(3,obj.panels(4,i))-obj.grid_init(3,obj.panels(2,i)) obj.grid_init(2,obj.panels(2,i))-obj.grid_init(2,obj.panels(4,i))...
                    0 obj.grid_init(3,obj.panels(1,i))-obj.grid_init(3,obj.panels(3,i)) obj.grid_init(2,obj.panels(3,i))-obj.grid_init(2,obj.panels(1,i));...
                    obj.grid_init(3,obj.panels(4,i))-obj.grid_init(3,obj.panels(2,i)) 0 obj.grid_init(1,obj.panels(2,i))-obj.grid_init(1,obj.panels(4,i))...
                    obj.grid_init(3,obj.panels(1,i))-obj.grid_init(3,obj.panels(3,i)) 0 obj.grid_init(1,obj.panels(3,i))-obj.grid_init(1,obj.panels(1,i))...
                    obj.grid_init(3,obj.panels(2,i))-obj.grid_init(3,obj.panels(4,i)) 0 obj.grid_init(1,obj.panels(4,i))-obj.grid_init(1,obj.panels(2,i))...
                    obj.grid_init(3,obj.panels(3,i))-obj.grid_init(3,obj.panels(1,i)) 0 obj.grid_init(1,obj.panels(1,i))-obj.grid_init(1,obj.panels(3,i));...
                    obj.grid_init(2,obj.panels(2,i))-obj.grid_init(2,obj.panels(4,i)) obj.grid_init(1,obj.panels(4,i))-obj.grid_init(1,obj.panels(2,i)) 0 ...
                    obj.grid_init(2,obj.panels(3,i))-obj.grid_init(2,obj.panels(1,i)) obj.grid_init(1,obj.panels(1,i))-obj.grid_init(1,obj.panels(3,i)) 0 ...
                    obj.grid_init(2,obj.panels(4,i))-obj.grid_init(2,obj.panels(2,i)) obj.grid_init(1,obj.panels(2,i))-obj.grid_init(1,obj.panels(4,i)) 0 ...
                    obj.grid_init(2,obj.panels(1,i))-obj.grid_init(2,obj.panels(3,i)) obj.grid_init(1,obj.panels(2,i))-obj.grid_init(1,obj.panels(4,i)) 0];
                
                for j=1:4
                    d_norm_nvec(1,1+(j-1)*3:j*3)=sum(diag(n_vec)*d_nvec(:,1+(j-1)*3:j*3)/norm_n_vec);
                end
                
                d_colloc_nvec(:,:,i)=(d_nvec*norm_n_vec-n_vec*d_norm_nvec)/norm_n_vec^2;
                
                obj.colloc_nvec_init(:,i)=colloc_nvec_init(:,i);
                obj.d_colloc_nvec(:,:,i)=d_colloc_nvec(:,:,i);
                
                obj.colloc_init(:,i)=0.125*obj.grid_init(:,obj.panels(1,i))+0.125*obj.grid_init(:,obj.panels(2,i))+0.375*obj.grid_init(:,obj.panels(3,i))+0.375*obj.grid_init(:,obj.panels(4,i));
                obj.colloc(:,i)=0.125*obj.grid(:,obj.panels(1,i))+0.125*obj.grid(:,obj.panels(2,i))+0.375*obj.grid(:,obj.panels(3,i))+0.375*obj.grid(:,obj.panels(4,i));
                
            end
        end
        
        function obj=compute_panel_area(obj)
            for i=1:1:length(obj.panels)
                r1=obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(1,i));
                r2=obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(2,i));
                r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                obj.area(i)=0.5*norm(r1xr2);
            end
        end
        
        function wij=compute_downwash_wing_wing(obj,panel_idx,colloc_len)
            vortex=obj.grid_vring(:,obj.panels(1:4,panel_idx)) ;
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                for i=1:length(vortex)
                    if i==4
                        r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                        r2=obj.colloc(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                        r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB/(4*pi);
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,colloc_idx)');
            end
        end
        
        function [wij_x,wij_y,wij_z]=compute_downwash_wing_wing_sw(obj,panel_idx,colloc_len)
            vortex=obj.grid_vring(:,obj.panels(1:4,panel_idx)) ;
            skip3=0;
            for colloc_idx=1:1:colloc_len
%                 if colloc_idx==panel_idx
%                     start=1;
%                 else
%                     start=2;
%                 end
%                 if panel_idx>1
%                     if obj.is_te(panel_idx-1)
%                         skip3=1;
%                     else
%                         skip3=0;
%                     end
%                 end
                w=[0 0 0];
                for i=2:2:length(vortex)
                    if i==4
                        r1=obj.fvap(:,colloc_idx)-vortex(:,i);
                        r2=obj.fvap(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.fvap(:,colloc_idx)-vortex(:,i);
                        r2=obj.fvap(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    
                    if i==3 && skip3==1
                        wAB=[0 0 0];
                    else
                    
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    
                    if (sum(conj(r0).*r1dr2))<-10^5
                    (sum(conj(r0).*r1dr2))
                    end
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    end

                    if ~isnan(wAB(1))
                         if (norm(r1xr2)>0.01)
                            w=w+wAB/(4*pi);
                         else
                            xxx=0;
                         end
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wij_x(1,colloc_idx)=w(1);
                wij_y(1,colloc_idx)=w(2);
                wij_z(1,colloc_idx)=w(3);
            end
        end
        
        function wij=compute_downwash_wing_wing_ff(obj,panel_idx,colloc_len)
            vortex=obj.grid_vring(:,obj.panels(1:4,panel_idx)) ;
            for colloc_idx=1:1:colloc_len
                % w=[0 0 0];
                w=0;
                dist=sqrt(sum((0.25*sum(vortex,2)-obj.colloc(:,colloc_idx)).^2));
                a=0.5*norm((vortex(:,3)-vortex(:,1)))+0.5*norm((vortex(:,4)-vortex(:,2)));
                % check if far field
                if (dist/a)<2
                    a=0;
                    for i=1:length(vortex)
                        if i==4
                            r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                            r2=obj.colloc(:,colloc_idx)-vortex(:,1);
                            r0=vortex(:,1)-vortex(:,i);
                        else
                            r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                            r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                            r0=vortex(:,i+1)-vortex(:,i);
                        end
                        %r1xr2=cross(r1,r2);
                        r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        r1dr2=(r1/norm(r1)-r2/norm(r2));
                        wAB=0;%-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                        %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                        if ~isnan(wAB(1))
                            w=w+wAB/(4*pi);
                        end
                    end
                else
                    % far field
                    disp('farfield');
                    r1=vortex(:,3)-vortex(:,1);
                    r2=vortex(:,4)-vortex(:,2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r0=0.25*sum(vortex,2);
                    r=obj.colloc(:,colloc_idx);
                    r0=r-r0;
                    A=0.5*norm(r1xr2);
                    dotp=r1xr2(1)*r0(1)+r1xr2(2)*r0(2)+r1xr2(3)*r0(3);
                    w=-A;%/(4*pi*norm(r0)^3)*[3*dotp*r0(1)/norm(r0)^2-r1xr2(1);3*dotp*r0(2)/norm(r0)^2-r1xr2(2);3*dotp*r0(3)/norm(r0)^2-r1xr2(3)]';
                    w=dist;
                    %w=-A/(4*pi*norm(r0)^3)*[3*r0(1)*r0(1)/norm(r0)^2-1;3*r0(1)*r0(2)/norm(r0)^2;3*r0(1)*r0(3)/norm(r0)^2]';
                    %den=(4*pi)*((r0(1))^2+(r0(2))^2+r0(3)^2)^(5/2);
                    % w=-[3*A*r0(1)*r0(3);3*A*r0(2)*r0(3);-A*(r0(1))^2+(r0(2))^2-2*r0(3)^2]'/den;
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wij(1,colloc_idx)=w;%sum(conj(w).*obj.colloc_nvec(:,colloc_idx)');
            end
        end
        
        
        function wij=compute_downwash_wing_wing_mex(obj)
            [wij]=compute_influence_vring(obj.grid_vring,obj.panels,0,obj.colloc,obj.colloc_nvec,0);
            wij=wij/(4*pi);
        end
        
        function wij=compute_influence_vring_ff_in(obj,grid,panels,colloc,colloc_nvec)
            vortex=grid(:,panels(1:4,1)) ;
            w=[0 0 0];
            colloc_idx=1;
            dist=sqrt(sum((0.25*sum(vortex,2)).^2-colloc(:,1).^2));
            a=0.5*norm((vortex(:,3)-vortex(:,1)))+0.5*norm((vortex(:,4)-vortex(:,2)));
            % check if far field
            if dist/a<8
                for i=1:length(vortex)
                    if i==4
                        r1=colloc(:,colloc_idx)-vortex(:,i);
                        r2=colloc(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                        r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB/(4*pi);
                    end
                end
            else
                % far field
                r1=vortex(:,3)-vortex(:,1);
                r2=vortex(:,4)-vortex(:,2);
                r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                r0=0.25*sum(vortex,2);
                r=obj.colloc(:,colloc_idx);
                r0=r-r0;
                A=0.5*norm(r1xr2);
                dotp=r1xr2(1)*r0(1)+r1xr2(2)*r0(2)+r1xr2(3)*r0(3);
                w=-A/(4*pi*norm(r0)^3)*[3*dotp*r0(1)/norm(r0)^2-r1xr2(1);3*dotp*r0(2)/norm(r0)^2-r1xr2(2);3*dotp*r0(3)/norm(r0)^2-r1xr2(3)]';
                %w=-A/(4*pi*norm(r0)^3)*[3*r0(1)*r0(1)/norm(r0)^2-1;3*r0(1)*r0(2)/norm(r0)^2;3*r0(1)*r0(3)/norm(r0)^2]';
                %den=(4*pi)*((r0(1))^2+(r0(2))^2+r0(3)^2)^(5/2);
                % w=-[3*A*r0(1)*r0(3);3*A*r0(2)*r0(3);-A*(r0(1))^2+(r0(2))^2-2*r0(3)^2]'/den;
            end
            %wij=dot(w,obj.colloc_nvec(:,panel_idx));
            wij(1,colloc_idx)=sum(conj(w).*colloc_nvec(:,colloc_idx)');
        end
        %
        function wij=compute_downwash_wing_wake(obj,panel_wake_idx,colloc_len)
            vortex=obj.grid_wake(:,obj.panels_wake(1:4,panel_wake_idx));
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                for i=1:length(vortex)
                    if i==4
                        r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                        r2=obj.colloc(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                        r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB/(4*pi);
                    end
                end
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,colloc_idx)');
            end
        end
        
        function [wij_x,wij_y,wij_z]=compute_downwash_wing_wake3D(obj,panel_wake_idx,colloc_len)
            vortex=obj.grid_wake(:,obj.panels_wake(1:4,panel_wake_idx));
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                for i=2:2:length(vortex)
                    if i==4
                        r1=obj.fvap(:,colloc_idx)-vortex(:,i);
                        r2=obj.fvap(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.fvap(:,colloc_idx)-vortex(:,i);
                        r2=obj.fvap(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB/(4*pi);
                    end
                end
                wij_x(1,colloc_idx)=w(1);
                wij_y(1,colloc_idx)=w(2);
                wij_z(1,colloc_idx)=w(3);
            end
        end
        
        function wij=compute_downwash_wing_wake_mex(obj)
            [wij]=compute_influence_vring(obj.grid_wake,obj.panels_wake,0,obj.colloc,obj.colloc_nvec,0);
            wij=wij/(4*pi);
        end
        
        function wij=compute_downwash_wake_wing(obj,panel_idx,colloc_len)
            vortex=obj.grid_vring(:,obj.panels(1:4,panel_idx));
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                r=sqrt((0.25*sum(vortex,2)).^2-obj.grid_wake(:,colloc_idx).^2);
                a=0.5*sqrt((vortex(:,3)-vortex(:,1)).^2)+0.5*sqrt((vortex(:,3)-vortex(:,1)).^2);
                % check if far field
                if norm(r)/norm(a)<3
                    for i=1:length(vortex)
                        if i==4
                            r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                            r2=obj.grid_wake(:,colloc_idx)-vortex(:,1);
                            r0=vortex(:,1)-vortex(:,i);
                        else
                            r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                            r2=obj.grid_wake(:,colloc_idx)-vortex(:,i+1);
                            r0=vortex(:,i+1)-vortex(:,i);
                        end
                        %r1xr2=cross(r1,r2);
                        r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        r1dr2=(r1/norm(r1)-r2/norm(r2));
                        wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                        %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                        
                        r1xr0=[r1(2)*r0(3)-r1(3)*r0(2),r1(3)*r0(1)-r1(1)*r0(3),r1(1)*r0(2)-r1(2)*r0(1)];
                        coredist=norm((norm(r1xr0))/(norm(r0)));
                        
                        if ~isnan(wAB(1)) && (coredist>0.1)
                            w=w+wAB/(4*pi);
                        end
                    end
                else
                    % far field
                    r1=vortex(:,3)-vortex(:,1);
                    r2=vortex(:,4)-vortex(:,2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r0=0.25*sum(vortex,2);
                    r=obj.grid_wake(:,colloc_idx);
                    A=0.5*r1xr2;
                    den=(4*pi)*((r(1)-r0(1))^2+(r(2)-r0(2))^2+r(3)^2)^(5/2);
                    w=-[3*A*(r(1)-r0(1))*r(3);3*A*(r(2)-r0(2))*r(3);-A*((r(1)-r0(1))^2+(r(2)-r0(2))^2-2*r(3)^2)]/den;
                end
                wij(1,colloc_idx)=w(3);%sum(conj(w).*obj.colloc_wake_nvec(:,colloc_idx)');
            end
        end
        
        function wij=compute_downwash_wake_wake(obj,panel_wake_idx,colloc_len)
            vortex=obj.grid_wake(:,obj.panels_wake(1:4,panel_wake_idx));
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                for i=1:length(vortex)
                    if i==4
                        r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                        r2=obj.grid_wake(:,colloc_idx)-vortex(:,1);
                        r0=vortex(:,1)-vortex(:,i);
                    else
                        r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                        r2=obj.grid_wake(:,colloc_idx)-vortex(:,i+1);
                        r0=vortex(:,i+1)-vortex(:,i);
                    end
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    r1xr0=[r1(2)*r0(3)-r1(3)*r0(2),r1(3)*r0(1)-r1(1)*r0(3),r1(1)*r0(2)-r1(2)*r0(1)];
                    coredist=norm((norm(r1xr0))/(norm(r0)));
                    
                    if ~isnan(wAB(1)) && (coredist>0.2)
                        w=w+wAB/(4*pi);
                    end
                end
                wij(1,colloc_idx)=w(3);%sum(conj(w).*obj.colloc_wake_nvec(:,colloc_idx)');
            end
        end
        
        function wij=compute_downwash_wake_wake_ff(obj,panel_wake_idx,colloc_len)
            vortex=obj.grid_wake(:,obj.panels_wake(1:4,panel_wake_idx));
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                dist=sqrt((0.25*sum(vortex,2)).^2-obj.grid_wake(:,colloc_idx).^2);
                a=0.5*sqrt((vortex(:,3)-vortex(:,1)).^2)+0.5*sqrt((vortex(:,4)-vortex(:,2)).^2);
                % check if far field
                if norm(dist)/norm(a)<12
                    
                    % near field
                    for i=1:length(vortex)
                        if i==4
                            r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                            r2=obj.grid_wake(:,colloc_idx)-vortex(:,1);
                            r0=vortex(:,1)-vortex(:,i);
                        else
                            r1=obj.grid_wake(:,colloc_idx)-vortex(:,i);
                            r2=obj.grid_wake(:,colloc_idx)-vortex(:,i+1);
                            r0=vortex(:,i+1)-vortex(:,i);
                        end
                        %r1xr2=cross(r1,r2);
                        r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                        r1dr2=(r1/norm(r1)-r2/norm(r2));
                        wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                        %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                        
                        r1xr0=[r1(2)*r0(3)-r1(3)*r0(2),r1(3)*r0(1)-r1(1)*r0(3),r1(1)*r0(2)-r1(2)*r0(1)];
                        coredist=norm((norm(r1xr0))/(norm(r0)));
                        
                        if ~isnan(wAB(1)) && (coredist>0.2)
                            w=w+wAB/(4*pi);
                        end
                    end
                else
                    % far field
                    r1=vortex(:,3)-vortex(:,1);
                    r2=vortex(:,4)-vortex(:,2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r0=0.25*sum(vortex,2);
                    r=obj.grid_wake(:,colloc_idx);
                    r0=r-r0;
                    A=0.5*norm(r1xr2);
                    dotp=r1xr2(1)*r0(1)+r1xr2(2)*r0(2)+r1xr2(3)*r0(3);
                    w=-A/(4*pi*norm(r0)^3)*[3*dotp*r0(1)/norm(r0)^2-r1xr2(1);3*dotp*r0(2)/norm(r0)^2-r1xr2(2);3*dotp*r0(3)/norm(r0)^2-r1xr2(3)];
                    % den=(4*pi)*((r(1)-r0(1))^2+(r(2)-r0(2))^2+r(3)^2)^(5/2);
                    %w=[3*A*(r(1)-r0(1))*r(3);3*A*(r(2)-r0(2))*r(3);-A*((r(1)-r0(1))^2+(r(2)-r0(2))^2-2*r(3)^2)]'/den;
                end
                wij(1,colloc_idx)=w(3);%sum(conj(w).*obj.colloc_wake_nvec(:,colloc_idx)');
            end
        end
        
        function wij=compute_downwash_mex(obj,panel_idx)
            
        end
        function obj=compute_influence_coeff_matrix_wake(obj)
            n=length(obj.colloc);
            nw=size(obj.panels_wake,2);
            ng=size(obj.grid_wake,2);
            Cbw=zeros(n,nw);
            
            %% TODO: check problem
            %% influence coefficent from wake grid to bound collication points
            Cbw=compute_influence_vring(obj.grid_wake,obj.panels_wake,0,obj.colloc,obj.colloc_nvec,0);
            %Cbw=compute_influence_vring_ff(obj.grid_wake,obj.panels_wake,0,obj.colloc,obj.colloc_nvec,0);
            obj.Abw=Cbw/(4*pi);
            
     
            
            Cbb_x2=zeros(n,n);
            Cbb_y2=zeros(n,n);
            Cbb_z2=zeros(n,n);
            
            [Cbb_x2,Cbb_y2,Cbb_z2]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.fvap,0,0);
            obj.Abw_x=Cbb_x2/(4*pi);
            obj.Abw_y=Cbb_y2/(4*pi);
            obj.Abw_z=Cbb_z2/(4*pi);
            
%             [Cbb_x2,Cbb_y2,Cbb_z2]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.grid_wake,0,0);
%             obj.Aww_x=Cbb_x2/(4*pi);
%             obj.Aww_y=Cbb_y2/(4*pi);
%             obj.Aww_z=Cbb_z2/(4*pi);
%             
%             [Cbb_x,Cbb_y,Cbb_z]=compute_influence_vring_3D(obj.grid_vring,obj.panels,0,obj.grid_wake,0,0);
%             obj.Awb_x=Cbb_x/(4*pi);
%             obj.Awb_y=Cbb_y/(4*pi);
%             obj.Awb_z=Cbb_z/(4*pi);
% %             
%             Cbb_x=zeros(n,n);
%             Cbb_y=zeros(n,n);
%             Cbb_z=zeros(n,n);
%             for i=1:n
%                 [Cbb_x(1:n,i),Cbb_y(1:n,i),Cbb_z(1:n,i)]=obj.compute_downwash_wing_wing_sw(i,n);
%             end
%             obj.Abb_x=Cbb_x;
%             obj.Abb_y=Cbb_y;
%             obj.Abb_z=Cbb_z;

%                             [Aww_x,Aww_y,Aww_z]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.colloc,0,0);
%                             obj.Aww_x=Aww_x/(4*pi);
%                             obj.Aww_y=Aww_y/(4*pi);
%                             obj.Aww_z=Aww_z/(4*pi);
          end
        
        function obj=compute_influence_coeff_matrix(obj)
            n=length(obj.colloc);
            nw=size(obj.panels_wake,2);
            ng=size(obj.grid_wake,2);
            
            C=zeros(n,n);
            %% influence coefficent from bound to bound
            C=compute_influence_vring(obj.grid_vring,obj.panels,0,obj.colloc,obj.colloc_nvec,0);
            %C=compute_influence_vring_ff(obj.grid_vring,obj.panels,0,obj.colloc,obj.colloc_nvec,0);
            obj.Abb=C/(4*pi);
            
            %            for i=1:n
            %               [C(1:n,i)]=obj.compute_downwash_wing_wing_ff(i,n);
            %            end
            %            obj.Abb=C;
            %
            Cbw=zeros(n,nw);
            
            %% TODO: check problem
            %% influence coefficent from wake grid to bound collication points
            Cbw=compute_influence_vring(obj.grid_wake,obj.panels_wake,0,obj.colloc,obj.colloc_nvec,0);
            %Cbw=compute_influence_vring_ff(obj.grid_wake,obj.panels_wake,0,obj.colloc,obj.colloc_nvec,0);
            obj.Abw=Cbw/(4*pi);
            
            %             Cbw=zeros(n,nw);
            %             for i=1:nw
            %                 [Cbw(1:n,i)]=obj.compute_downwash_wing_wake(i,n);
            %             end
            %             obj.Abw=Cbw;
            
            Cwb=zeros(ng,n);
            Cwb=compute_influence_vring_wake(obj.grid_vring,obj.panels,0,obj.grid_wake,0,0);
            %Cwb=compute_influence_vring_wake_ff(obj.grid_vring,obj.panels,0,obj.grid_wake,0,0);
            obj.Awb=Cwb/(4*pi);
            
            %             for i=1:n
            %                             [Cwb(1:ng,i)]=obj.compute_downwash_wake_wing(i,ng);
            %             end
            %             obj.Awb=Cwb;
            %
            Cww=zeros(ng,nw);
            Cww=compute_influence_vring_wake(obj.grid_wake,obj.panels_wake,0,obj.grid_wake,0,0);
            %Cww=compute_influence_vring_wake_ff(obj.grid_wake,obj.panels_wake,0,obj.grid_wake,0,0);
            % time_normal=toc
            obj.Aww=Cww/(4*pi);
            %
            %                         for i=1:nw
            %                             [Cww(1:ng,i)]=obj.compute_downwash_wake_wake(i,ng);
            %                         end
            %
            %                         obj.Aww=Cww;
            %             xxx=1;
            %             if xxx==1
            
% 
%             Cbb_x=zeros(n,n);
%             Cbb_y=zeros(n,n);
%             Cbb_z=zeros(n,n);
%             for i=1:nw
%                 [Cbb_x(1:n,i),Cbb_y(1:n,i),Cbb_z(1:n,i)]=obj.compute_downwash_wing_wake3D(i,n);
%             end
%             obj.Abw_x=Cbb_x;
%             obj.Abw_y=Cbb_y;
%             obj.Abw_z=Cbb_z;
            
            Cbb_x2=zeros(n,n);
            Cbb_y2=zeros(n,n);
            Cbb_z2=zeros(n,n);
            
            [Cbb_x2,Cbb_y2,Cbb_z2]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.fvap,0,0);
            obj.Abw_x=Cbb_x2/(4*pi);
            obj.Abw_y=Cbb_y2/(4*pi);
            obj.Abw_z=Cbb_z2/(4*pi);
            
            
            Cbb_x=zeros(n,n);
            Cbb_y=zeros(n,n);
            Cbb_z=zeros(n,n);
            [Cbb_x,Cbb_y,Cbb_z]=compute_influence_vring_3D(obj.grid_vring,obj.panels,0,obj.fvap,0,0);
            obj.Abb_x=Cbb_x/(4*pi);
            obj.Abb_y=Cbb_y/(4*pi);
            obj.Abb_z=Cbb_z/(4*pi);
            
%             [Cbb_x2,Cbb_y2,Cbb_z2]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.grid_wake,0,0);
%             obj.Aww_x=Cbb_x2/(4*pi);
%             obj.Aww_y=Cbb_y2/(4*pi);
%             obj.Aww_z=Cbb_z2/(4*pi);
%             
%             [Cbb_x,Cbb_y,Cbb_z]=compute_influence_vring_3D(obj.grid_vring,obj.panels,0,obj.grid_wake,0,0);
%             obj.Awb_x=Cbb_x/(4*pi);
%             obj.Awb_y=Cbb_y/(4*pi);
%             obj.Awb_z=Cbb_z/(4*pi);
% %             
%             Cbb_x=zeros(n,n);
%             Cbb_y=zeros(n,n);
%             Cbb_z=zeros(n,n);
%             for i=1:n
%                 [Cbb_x(1:n,i),Cbb_y(1:n,i),Cbb_z(1:n,i)]=obj.compute_downwash_wing_wing_sw(i,n);
%             end
%             obj.Abb_x=Cbb_x;
%             obj.Abb_y=Cbb_y;
%             obj.Abb_z=Cbb_z;

%                             [Aww_x,Aww_y,Aww_z]=compute_influence_vring_wake3D(obj.grid_wake,obj.panels_wake,0,obj.colloc,0,0);
%                             obj.Aww_x=Aww_x/(4*pi);
%                             obj.Aww_y=Aww_y/(4*pi);
%                             obj.Aww_z=Aww_z/(4*pi);
        end
        
        function obj=determine_boundary_conditions(obj)
            
%             for i=1:1:length(obj.colloc)
%                 %obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf);% was slow
%                 obj.b(i)=-(obj.colloc_nvec(1,i)*obj.Uinf(1)+obj.colloc_nvec(2,i)*obj.Uinf(2)+obj.colloc_nvec(3,i)*obj.Uinf(3));
%             end
            % 245 times faster !!!!
            obj.b=-(obj.colloc_nvec'*obj.Uinf')';
        end
        
        function obj=determine_boundary_conditions_gust_deflection_linear(obj,prv_colloc,deflection_induced_speed_init,timestep)
            
            for i=1:1:length(obj.colloc)
                
                obj.colloc(:,i)=0.125*obj.grid(:,obj.panels(1,i))+0.125*obj.grid(:,obj.panels(2,i))+0.375*obj.grid(:,obj.panels(3,i))+0.375*obj.grid(:,obj.panels(4,i));
                
                deflection_induced_speed(:,i)=(-obj.colloc(:,i)+prv_colloc(:,i))/timestep;
                deflection_induced_speed(3,i)=deflection_induced_speed(3,i)*obj.Ma_corr;
                
                x_colloc=obj.colloc_init(1,i);
                
                Uinf_x=nakeinterp1(obj.x_kin',obj.Uinf_x(1,:)',x_colloc);
                Uinf_y=nakeinterp1(obj.x_kin',obj.Uinf_x(2,:)',x_colloc);
                Uinf_z=nakeinterp1(obj.x_kin',obj.Uinf_x(3,:)',x_colloc);
                
                delta_disp=[obj.grid(:,obj.panels(1,i)); obj.grid(:,obj.panels(2,i)); obj.grid(:,obj.panels(3,i)); obj.grid(:,obj.panels(4,i))];
                
                obj.colloc_nvec(:,i)=obj.colloc_nvec_init(:,i)+obj.d_colloc_nvec(:,:,i)*delta_disp;
                
                obj.b(i)=-(dot(obj.colloc_nvec_init(:,i),([Uinf_x; Uinf_y; Uinf_z]+deflection_induced_speed(:,i)))+dot(obj.d_colloc_nvec(:,:,i)*delta_disp,([Uinf_x; Uinf_y; Uinf_z]+deflection_induced_speed_init(:,i))));
                
            end
        end
        
        function obj=determine_boundary_conditions_gust_deflection_rotation(obj,p,q,r,prv_colloc,timestep)
            deflection_induced_speed=(-obj.colloc+prv_colloc)/timestep;
            dq=q;%*obj.Ma_corr;
            deflection_induced_speed(3,:)=deflection_induced_speed(3,:);%*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                x_colloc=obj.colloc(1,i);
                Uinf_x=nakeinterp1(obj.x_kin',obj.Uinf_x(1,:)',x_colloc);
                Uinf_y=nakeinterp1(obj.x_kin',obj.Uinf_x(2,:)',x_colloc);
                Uinf_z=nakeinterp1(obj.x_kin',obj.Uinf_x(3,:)',x_colloc);
                dV=cross(obj.r_dwn(:,i),[p dq r]); %check cross product
                obj.b(i)=-dot(obj.colloc_nvec(:,i),[Uinf_x Uinf_y Uinf_z]+deflection_induced_speed(:,i)'+dV);
            end
        end
        
        function obj=determine_boundary_conditions_deflection(obj,prv_colloc,timestep)
            deflection_induced_speed=(-obj.colloc+prv_colloc)/timestep;
            deflection_induced_speed(3,:)=deflection_induced_speed(3,:);%*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                % dot product very slow in matlab
                %  obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+deflection_induced_speed(:,i)');
                obj.b(i)=-(obj.colloc_nvec(1,i)*(obj.Uinf(1)+deflection_induced_speed(1,i)))-(obj.colloc_nvec(2,i)*(obj.Uinf(2)+deflection_induced_speed(2,i)))-(obj.colloc_nvec(3,i)*(obj.Uinf(3)+deflection_induced_speed(3,i)));
            end
        end
        
        function obj=determine_boundary_conditions_gust_deflection(obj,prv_colloc,timestep)
            deflection_induced_speed=(-obj.colloc+prv_colloc)/timestep;
            deflection_induced_speed(3,:)=deflection_induced_speed(3,:);%*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                x_colloc=obj.colloc(1,i);
                %                Uinf_x=lininterp1(obj.x_kin,obj.Uinf_x(1,:),x_colloc);
                %                Uinf_y=lininterp1(obj.x_kin,obj.Uinf_x(2,:),x_colloc);
                %                Uinf_z=lininterp1(obj.x_kin,obj.Uinf_x(3,:),x_colloc);
                Uinf_x=nakeinterp1(obj.x_kin',obj.Uinf_x(1,:)',x_colloc);
                Uinf_y=nakeinterp1(obj.x_kin',obj.Uinf_x(2,:)',x_colloc);
                Uinf_z=nakeinterp1(obj.x_kin',obj.Uinf_x(3,:)',x_colloc);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),[Uinf_x Uinf_y Uinf_z]+deflection_induced_speed(:,i)');
            end
        end
        
        
        
        function obj=determine_boundary_conditions_gust(obj)
            for i=1:1:length(obj.colloc)
                x_colloc=obj.colloc(1,i);
                %                 Uinf_x=lininterp1(obj.x_kin,obj.Uinf_x(1,:),x_colloc);
                %                 Uinf_y=lininterp1(obj.x_kin,obj.Uinf_x(2,:),x_colloc);
                %                 Uinf_z=lininterp1(obj.x_kin,obj.Uinf_x(3,:),x_colloc);
                Uinf_x=nakeinterp1(obj.x_kin',obj.Uinf_x(1,:)',x_colloc);
                Uinf_y=nakeinterp1(obj.x_kin',obj.Uinf_x(2,:)',x_colloc);
                Uinf_z=nakeinterp1(obj.x_kin',obj.Uinf_x(3,:)',x_colloc);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),[Uinf_x Uinf_y Uinf_z]);
            end
        end
        
        function obj=determine_boundary_conditions_q(obj,dtheta)
            dq=dtheta*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[0 dq 0]); %check cross product
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV);
            end
        end
        
        function obj=determine_boundary_conditions_p(obj,dtheta)
            dp=dtheta*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[dp 0 0]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV);
            end
        end
        
        function obj=determine_boundary_conditions_r(obj,dtheta)
            dr=dtheta*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[0 0 dr]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV);
            end
        end
        
        function obj=determine_boundary_conditions_pqr(obj,p,q,r)
            dq=q*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[p dq r]); %check cross product
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV);
            end
        end
        
        function obj=perform_wake_rollup(obj)
            
            gamma_scaling_factor=zeros(length(obj.panels_wake),1);
            circum_prev=zeros(1,length(obj.panels_wake));
            circum_nxt=zeros(1,length(obj.panels_wake));
            for i=1:length(obj.panels_wake)
                p1=obj.grid_wake(:,obj.panels_wake(1,i));
                p2=obj.grid_wake(:,obj.panels_wake(2,i));
                p3=obj.grid_wake(:,obj.panels_wake(3,i));
                p4=obj.grid_wake(:,obj.panels_wake(4,i));
                circum_prev(i)=norm(p2-p1)+norm(p3-p2)+norm(p4-p3)+norm(p1-p4);
            end
            %             obj.uind=-obj.Awb_x*obj.Gamma-obj.Aww_x*obj.Gamma_wake;
            %             obj.grid_wake(1,obj.row_length+1:end)=obj.grid_wake(1,obj.row_length+1:end)+obj.uind(obj.row_length+1:end)'*obj.t_step;
            %             obj.uind=-obj.Awb_y*obj.Gamma-obj.Aww_y*obj.Gamma_wake;
            %             obj.grid_wake(2,obj.row_length+1:end)=obj.grid_wake(2,obj.row_length+1:end)+obj.uind(obj.row_length+1:end)'*obj.t_step;
            obj.uind=-obj.Awb_z*obj.Gamma-obj.Aww_z*obj.Gamma_wake;%%+
            obj.grid_wake(3,obj.row_length+1:end)=obj.grid_wake(3,obj.row_length+1:end)+obj.uind(obj.row_length+1:end)'*obj.t_step;
            
            for i=1:length(obj.panels_wake)
                p1=obj.grid_wake(:,obj.panels_wake(1,i));
                p2=obj.grid_wake(:,obj.panels_wake(2,i));
                p3=obj.grid_wake(:,obj.panels_wake(3,i));
                p4=obj.grid_wake(:,obj.panels_wake(4,i));
                circum_nxt(i)=norm(p2-p1)+norm(p3-p2)+norm(p4-p3)+norm(p1-p4);
                gamma_scaling_factor(i)=circum_prev(i)/circum_nxt(i);
            end
            obj.Gamma_wake=obj.Gamma_wake.*gamma_scaling_factor;
        end
        
        function obj=shed_wake(obj)
            for i=obj.n_step:-1:1
                obj.grid_wake(1,1+obj.row_length*i:obj.row_length*(i+1))=obj.grid_wake(1,1+obj.row_length*(i-1):obj.row_length*i)+obj.Uinf(1)*obj.t_step;%/obj.Ma_corr;
                obj.grid_wake(2,1+obj.row_length*i:obj.row_length*(i+1))=obj.grid_wake(2,1+obj.row_length*(i-1):obj.row_length*i)+obj.Uinf(2)*obj.t_step;
                obj.grid_wake(3,1+obj.row_length*i:obj.row_length*(i+1))=obj.grid_wake(3,1+obj.row_length*(i-1):obj.row_length*i)+obj.Uinf(3)*obj.t_step;
            end
        end
        
        function obj=shed_wake_gust(obj)
            for i=obj.n_step:-1:1
                for j=1:obj.row_length
                    x_colloc=obj.grid_wake(1,1+obj.row_length*(i-1)-1+j);
                    %  Uinf_x=lininterp1(obj.x_kin,obj.Uinf_x(1,:),x_colloc);
                    Uinf_x=obj.Uinf_x(1,1);
                    % Uinf_y=lininterp1(obj.x_kin,obj.Uinf_x(2,:),x_colloc);
                    Uinf_y=obj.Uinf_x(2,1);
                    %Uinf_z=lininterp1(obj.x_kin,obj.Uinf_x(3,:),x_colloc);
                    Uinf_z=nakeinterp1(obj.x_kin',obj.Uinf_x(3,:)',x_colloc);
                    obj.grid_wake(1,1+obj.row_length*i+j-1)=obj.grid_wake(1,1+obj.row_length*(i-1)+j-1)+Uinf_x*obj.t_step;%/obj.Ma_corr;
                    obj.grid_wake(2,1+obj.row_length*i+j-1)=obj.grid_wake(2,1+obj.row_length*(i-1)+j-1)+Uinf_y*obj.t_step;
                    obj.grid_wake(3,1+obj.row_length*i+j-1)=obj.grid_wake(3,1+obj.row_length*(i-1)+j-1)+Uinf_z*obj.t_step;
                end
            end
        end
        
        function obj=initialize_unsteady_computation_settings(obj,n_args,args)
            obj.settings.modal_data=0;
 
            if obj.Ma_corr<1
                obj.Uinf=a2bf(obj.state.V_A,obj.state.alpha,obj.state.beta,obj.Ma_corr);
            else
                obj.Uinf=a2bf(obj.state.V_A,obj.state.alpha,obj.state.beta,1);
            end
            % check if structural mode coupling coefficients are desired
            if n_args>=6
                obj.settings.modal_data=1;
            end
            
            if obj.settings.modal_data==1
                obj.settings.n_mode=args{1};
            end
        end
        
        function obj=initialize_unsteady_computation(obj)
            
            % initialize bound and wake grid
            obj=obj.initialize_vring_grid(obj.t_step);
            
            % compute required number of streamwise wake panels
            n_wake=obj.settings.wakelength_factor*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
            
            obj=obj.initialize_wake_grid(ceil(n_wake),obj.t_step);
            
            % initialize vorticity vector for bound and wake grid
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            obj.Gamma=zeros(length(obj.panels),1);
            
            % compute initial influence coefficient matrix
            if isempty(obj.Abb)
                obj=obj.compute_influence_coeff_matrix();
            elseif(size(obj.Aww,2)~=length(obj.panels_wake))
                obj=obj.compute_influence_coeff_matrix();
            end
            % initialize variables
            
            obj.CX=[];
            obj.CY=[];
            obj.CZ=[];
            obj.Cl=[];
            obj.CL=[];
            obj.CM=[];
            obj.CN=[];
            obj.Cdi=[];
            obj.Cdi2=[];
            obj.alpha_s=[];
            % initialization in case modal coefficients are desired
            if obj.settings.modal_data==1
                obj.C_modes=zeros(obj.settings.n_mode,length(obj.t_vec));
            end
        end
        
        function obj=initialize_wake_circulation_pqr(obj, dp,dq,dr)
            % here the wake circulation is initialized until a steady value
            % is reached
            obj=obj.determine_boundary_conditions_pqr(dp,dq,dr);
            obj.Gamma=linsolve(obj.Abb,obj.b');
            obj.Gamma_wake=obj.Cbb*obj.Gamma;
            i=0;
            while abs((sum(abs(obj.Gamma_prv))-sum(abs(obj.Gamma)))/sum(abs(obj.Gamma)))>10^(-10)
                obj=obj.determine_boundary_conditions_pqr(dp,dq,dr);
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                i=i+1;
            end
        end
        
        function obj=complex_coeffs_from_timeseries(obj,omega)
            %TODO: Reduce code here, multiple executions -> function?
            % compute complex aerodynamic coefficents by curvefit and
            % transformation into complex plane
            t_vec=obj.t_vec;
            %options = optimoptions('lsqcurvefit','TolX',1E-14,'TolFun',1E-14);
            options = optimset('TolX',1E-20,'TolFun',1E-20);
            amp_lim=1e-8;
            LB=[0 -pi -Inf];
            UB=[Inf pi Inf];
            
            %% CX
            Init_Amplitude=(abs(max(obj.CX(length(t_vec)/4:end)))+abs(min(obj.CX(length(t_vec)/4:end))))/2;
            Init_Phase=-1.5;
            Init_Offset=0.5*max(obj.CX(length(t_vec)/4:end))+0.5*min(obj.CX(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CX(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            %x=x.*(abs(x)>obj.settings.coeff_eps);
            obj.CX_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CX over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CX,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %% CY
            Init_Amplitude=(abs(max(obj.CY(length(t_vec)/4:end)))+abs(min(obj.CY(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.CY(length(t_vec)/4:end))+0.5*min(obj.CY(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CY(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.CY_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CY,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %% CZ
            Init_Amplitude=(abs(max(obj.CZ(length(t_vec)/4:end)))+abs(min(obj.CZ(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.CZ(length(t_vec)/4:end))+0.5*min(obj.CZ(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CZ(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.CZ_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CZ over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CZ,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            %% Cl

            Init_Amplitude=(abs(max(obj.Cl(length(t_vec)/4:end)))+abs(min(obj.Cl(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.Cl(length(t_vec)/4:end))+0.5*min(obj.Cl(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.Cl(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.Cl_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CZ over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.Cl,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %% CL
            Init_Amplitude=(abs(max(obj.CL(length(t_vec)/4:end)))+abs(min(obj.CL(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.CL(length(t_vec)/4:end))+0.5*min(obj.CL(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CL(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.CL_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CM over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CL,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %% CM
            Init_Amplitude=(abs(max(obj.CM(length(t_vec)/4:end)))+abs(min(obj.CM(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.CM(length(t_vec)/4:end))+0.5*min(obj.CM(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CM(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.CM_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CM over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CM,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %% CN
            Init_Amplitude=(abs(max(obj.CN(length(t_vec)/4:end)))+abs(min(obj.CN(length(t_vec)/4:end))))/2;
            Init_Phase=0;
            Init_Offset=0.5*max(obj.CN(length(t_vec)/4:end))+0.5*min(obj.CN(length(t_vec)/4:end));
            x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.CN(length(t_vec)/2:end),LB,UB,options);
            if x(1)<amp_lim
                x(1)=0;
                x(2)=0;
            end
            obj.CN_complex=[real(x(1)*exp(1i*(x(2)-pi/2))) imag(x(1)*exp(1i*(x(2)-pi/2))) x(3)];
            % if obj.settings.debug mode is selected plot CM over time
            if obj.settings.debug==1
                figure
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.CN,'r');
                hold on
                plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
            end
            
            %%
            if obj.settings.modal_data==1
                % TODO check LSQcurvefitting here. does not work in some
                % cases, check also initial phase estimation
                for mod=1:obj.settings.n_mode
                    Init_Amplitude=(abs(max(obj.C_modes(mod,length(t_vec)/4:end)))+abs(min(obj.C_modes(mod,length(t_vec)/4:end))))/2;
                    Init_Phase=0;
                    Init_Offset=0.5*max(obj.C_modes(mod,length(t_vec)/4:end))+0.5*min(obj.C_modes(mod,length(t_vec)/4:end));
                    x = lsqcurvefit(@(x,t_vec) x(1)*cos(omega*t_vec+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],t_vec(length(t_vec)/2:end),obj.C_modes(mod,(length(t_vec))/2:end),LB,UB,options);
                    if x(1)<amp_lim
                        x(1)=0;
                        x(2)=0;
                    end
                    %x=x.*(abs(x)>obj.settings.coeff_eps);
                    %x(2)*180/pi[x(1)*cos(x(2)) x(1)*sin(x(2)) x(3)];
                    obj.C_MODES_complex(mod,1)=real(x(1)*exp(1i*(x(2)-pi/2)));
                    obj.C_MODES_complex(mod,2)=imag(x(1)*exp(1i*(x(2)-pi/2)));
                    obj.C_MODES_complex(mod,3)=[x(3)];
                    %amplitude=x(1)
                    %phase=x(2)*180/pi
                    % if obj.settings.debug mode is selected plot C_MODES over time
                    if obj.settings.debug==1
                        figure
                        plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,obj.C_modes(mod,:),'r');
                        hold on
                        plot(norm(obj.Uinf)*t_vec/obj.reference.c_ref,x(1)*cos(omega*t_vec+x(2))+x(3),'b');
                    end
                end
            end
        end
        
        function obj=solve(obj)
            
            obj.t_step=obj.reference.c_ref/(norm(obj.Uinf)*16);
            %obj.t_step=obj.reference.c_ref/(norm(obj.Uinf)*4);
            n_wake=2*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
            
            obj=obj.initialize_vring_grid(obj.t_step);
            obj=obj.initialize_wake_grid(ceil(n_wake),obj.t_step);
            
            obj=obj.compute_influence_coeff_matrix();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.Abb,obj.b');
            err=10E10;
            prv_err=0;
            obj.Gamma_wake=obj.Cbb*obj.Gamma;
            obj.t_vec=[];
            %         obj.plot_gamma
            obj.CZ=[];
            while abs(prv_err-err)>0.1
                prv_err=err;
                err=sum(abs(obj.Gamma_wake));
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                
                % recompute influence coefficient matrices
                
                obj=obj.f_postprocess;
                obj=obj.shed_wake();
                
                if ~isempty(obj.t_vec)
                    obj.t_vec=[obj.t_vec obj.t_vec(end)+obj.t_step];
                else
                    obj.t_vec=[0];
                end
                %               close all
                % obj.plot_gamma
                %               drawnow
                % recompute influence coefficient matrices
                prv_err-err
            end
            %             for i =1:10
            %                 obj=obj.compute_influence_coeff_matrix();
            %                 obj=obj.perform_wake_rollup;
            %                 figure
            %                 obj.plot_gamma
            %                 drawnow
            %             end
            %
            %             obj=obj.perform_wake_rollup;
            %b_check=obj.Abb*obj.Gamma+obj.Abw*obj.Gamma_wake-obj.b;
            % plot(obj.CZ)
        end
        
        function obj=solve_unsteady(obj)
            %   obj=obj.compute_influence_coeff_matrix();
            
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            t_vec=0:obj.t_step:0.5;
            
            alpha_vec=5+0*sin(2*pi*10*t_vec);
            for i=1:length(t_vec)
                obj.Uinf=[norm(obj.Uinf)*cosd(alpha_vec(i)) 0 norm(obj.Uinf)*sind(alpha_vec(i))];
                obj=obj.determine_boundary_conditions();
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj.Gamma_wake=obj.Cbb*obj.Gamma+obj.Cbw*obj.Gamma_wake;
                obj=obj.shed_wake();
                obj=obj.perform_wake_rollup;
                % recompute influence coefficient matrices
                % obj=obj.compute_influence_coeff_matrix();
                t_vec(i)
                close all
                obj.plot_cp_wake
            end
        end
        
        function obj=solve_unsteady_gust(obj)
            filename = ['gust_aero_complete4.gif'];
            steps=30;
            
            %% (1-cos)gust
            obj.Uinf_x(1,:)=obj.Uinf(1);
            Uds=30;
            H=10;
            s=0:(2*H/(20)):2*H;
            gust_shape=Uds/2*(1-cos(pi*s/H));
            
            x_min=min(obj.grid(1,:));
            x_max=max(obj.grid_wake(1,:));
            
            obj.x_kin=x_min:(x_max-x_min)/(steps-1):x_max;
            obj.Uinf_x=zeros(3,steps);
            obj.Uinf_x(1,:)=obj.Uinf(1);
            obj.Uinf_x(2,:)=obj.Uinf(2);
            obj.Uinf_x(3,:)=obj.Uinf(3);
            
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            t_vec=0:obj.t_step:0.45;
            
            delta_x=obj.Uinf(1)*obj.t_step;
            gust_shape_interp=interp1(s,gust_shape,0:delta_x:2*H);
            
            start=30;
            n=1;
            k=1;
            m=1;
            for i=1:length(t_vec)
                if i>=start
                    obj.Uinf_x(1,:)=obj.Uinf(1);
                    obj.Uinf_x(2,:)=obj.Uinf(2);
                    obj.Uinf_x(3,:)=obj.Uinf(3);
                    if n<length(gust_shape_interp)
                        obj.Uinf_x(3,1:n)=obj.Uinf_x(3,k:k+n-1)+gust_shape_interp(n:-1:1);
                        n=n+1;
                    elseif k+n-1>length(obj.x_kin)
                        obj.Uinf_x(3,k:end)=obj.Uinf_x(3,k:end)+gust_shape_interp(1:length(obj.x_kin)-k+1);
                        k=k+1;
                    else
                        obj.Uinf_x(3,k:k+n-1)=obj.Uinf_x(3,k:k+n-1)+gust_shape_interp(1:n);
                        k=k+1;
                    end
                elseif k>=length(obj.x_kin)
                    obj.Uinf_x(1,:)=obj.Uinf(1);
                    obj.Uinf_x(2,:)=obj.Uinf(2);
                    obj.Uinf_x(3,:)=obj.Uinf(3);
                end
                obj=obj.determine_boundary_conditions_gust();
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                obj=obj.shed_wake_gust();
                obj=obj.perform_wake_rollup;
                % recompute influence coefficient matrices
                %obj=obj.compute_influence_coeff_matrix();
                t_vec(i)
                
                close all
                hfig=figure(33);
                set(hfig, 'Position', [0 0 1024*1.5 768*1.5])
                
                obj.plot_gamma;
                
                %plot(obj.x_kin,obj.Uinf_x(3,:))
                %view(3)
                %view(-30,-30);
                %xlim([-2 45])
                %ylim([-20 20]);
                
                ylim([0 60]);
                %zlim([-0.5 10]);
                %caxis([-30 30])
                drawnow
                frame = getframe(33);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if m == 1;
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
                end
                m=m+1;
            end
        end
        
        
        function obj=time_domain_flutter_analysis(obj,aircraft,aircraft_structure)
            
            n_beam=1;
            
            obj.t_step=obj.reference.c_ref/(norm(obj.Uinf)*16);
            
            n_wake=1*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
            
            obj=obj.initialize_vring_grid(obj.t_step);
            
            obj=obj.initialize_wake_grid(ceil(n_wake),obj.t_step);
            
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            
            %t_vec=0:obj.t_step:1;
            
            obj=obj.compute_influence_coeff_matrix();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.Abb,obj.b');
            obj.Gamma_wake=obj.Cbb*obj.Gamma;
            obj=obj.f_postprocess;
            
            % initialize newmark beta timemarching scheme
            
            beta=1/4;
            gamma=1/2;
            dt=obj.t_step;
            
            x_now=zeros(length(aircraft_structure.Ftest),1);
            xdot_now=zeros(length(aircraft_structure.Ftest),1);
            xdotdot_now=zeros(length(aircraft_structure.Ftest),1);
            x_nxt=zeros(length(aircraft_structure.Ftest),1);
            xdot_nxt=zeros(length(aircraft_structure.Ftest),1);
            xdotdot_nxt=zeros(length(aircraft_structure.Ftest),1);
            
            xdotdot_now=linsolve(aircraft_structure.Mff,aircraft_structure.Ftest*0-aircraft_structure.Kff*x_now);
            
            n=1;
            k=1;
            m=1;
            
            stop=0;
            dUinf=0.00*obj.Uinf;
            startUinf=obj.Uinf;
            
            while stop==0
                tic
                
                %% UPDATE XKIN
                if m>20
                    obj.Uinf= startUinf+1*[rand(1)-0.5 rand(1)-0.5 rand(1)-0.5];
                end
                
                % obj.t_step=obj.reference.c_ref/(norm(obj.Uinf)*6);
                %dt=obj.t_step;
                prv_colloc=obj.colloc;
                
                aircraft=aircraft.compute_beam_forces(obj.F_body,aircraft_structure);
                for bb=1:n_beam
                    aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
                end
                aircraft_structure=aircraft_structure.f_solve();
                
                x_nxt=linsolve(aircraft_structure.Mff+beta*dt^2*aircraft_structure.Kff,beta*dt^2*aircraft_structure.Ftest+aircraft_structure.Mff*x_now+dt*aircraft_structure.Mff*xdot_now+dt^2*aircraft_structure.Mff*(1/2-beta)*xdotdot_now);
                
                xdotdot_nxt=1/(beta*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta)*xdotdot_now);
                xdot_nxt=xdot_now+dt*((1-gamma)*xdotdot_now+gamma*xdotdot_nxt);
                
                aircraft_structure.nodal_deflections=x_nxt;
                
                if max(x_nxt(3:6:end))>3
                    stop=1;
                end
                
                aircraft_structure=aircraft_structure.f_postprocess();
                
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                
                obj=obj.update_grid();
                
                obj=obj.determine_boundary_conditions_deflection(prv_colloc,obj.t_step);
                %                obj=obj.determine_boundary_conditions_gust_deflection(prv_colloc,obj.t_step);
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj=obj.f_postprocess;
                obj=obj.shed_wake();
                
                obj.write_tecplot(['results/flutter/flutter_response_' num2str(m)],1);
                %                obj=obj.shed_wake_gust();
                % obj=obj.perform_wake_rollup;
                
                %                 close all
                %                 hfig=figure(33);
                %                 set(hfig, 'Position', [0 0 1024*1.5 768*1.5])
                %
                %                 obj.plot_gamma;
                %                 %plot(obj.x_kin,obj.Uinf_x(3,:))
                %                 view(3)
                %                 %view(-30,-30);
                %                 %                 xlim([-2 8])
                % %                 ylim([-2 2]);
                % %                 ylim([0 60]);
                % %                 zlim([-0.5 2]);
                %
                %                 xlim([-2 10])
                %                 ylim([-7 7]);
                %                 zlim([-0.5 10]);
                %                 caxis([-100 -5])
                %                 drawnow
                %                 frame = getframe(33);
                %                 im = frame2im(frame);
                %                 [imind,cm] = rgb2ind(im,256);
                %                 if m == 1;
                %                     imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
                %                 else
                %                     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
                %                 end
                m=m+1;
                %
                
                x_now=x_nxt;
                xdot_now=xdot_nxt;
                xdotdot_now=xdotdot_nxt;
                t_step=toc
            end
        end
        
        
        
        function obj=solve_unsteady_gust_mode(obj,amplitude,k,varargin)
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
            end
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            % compute pitching amplitude in [rad]
            gust_amplitude=amplitude;
            % initialize pitch rate for each timestep
            obj=obj.initialize_unsteady_computation();
            
            fprintf(['UVLM Gust Mode Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['gust_k=' num2str(k)]);
            end
            
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure);
            end
            
            x_min=min(obj.grid(1,:));
            x_max=max(obj.grid_wake(1,:));
            
            t_period=2*pi/omega;
            obj.x_kin=x_min:obj.t_step*obj.Uinf(1):x_max;
            t_kin=obj.x_kin/obj.Uinf(1);
            gust_shape=gust_amplitude*sin(omega*t_kin);
            
            obj.Uinf_x=zeros(3,length(obj.x_kin));
            obj.Uinf_x(1,:)=obj.Uinf(1);
            obj.Uinf_x(2,:)=obj.Uinf(2);
            obj.Uinf_x(3,:)=obj.Uinf(3)+gust_shape;
            % figure
            
            for i=1:length(obj.t_vec)
                gust_shape=gust_amplitude*sin(omega*t_kin-2*pi/t_period*obj.t_step*i);
                obj.Uinf_x(3,:)=obj.Uinf(3)+gust_shape;
                
                %                 plot(t_kin*obj.Uinf(1),gust_shape);
                %                 drawnow
                % compute boundary conditions in presence of gust mode
                obj=obj.determine_boundary_conditions_gust();
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                obj=obj.shed_wake_gust();
                % display current timestep
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                % postprocess timestep, compute aerodynamic coefficients
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/gust_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/obj.reference.S_ref;
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_roll(obj,amplitude,k,varargin)
            
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                shape=varargin{4};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end
            
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            % compute pitching amplitude in [rad]
            rolling_amplitude=amplitude*pi/180;
            % initialize pitch rate for each timestep
            phidot_vec=rolling_amplitude*omega*cos(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            
            alpha=obj.state.alpha*pi/180;
            fprintf(['UVLM Roll Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            %% COMPUTATION
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['roll_k=' num2str(k)]);
            end
            
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            % UVLM time marching computation
            for i=1:length(obj.t_vec)
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(phidot_vec(i)+phidot_vec(i+1))/2;
                else
                   bc=phidot_vec(i);
                end 
%                 bc=phidot_vec(i);
                obj.Uinf=norm(obj.Uinf)*[cos(alpha) 0 sin(alpha)];
                
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(bc,0,0);
                end
                % compute boundary conditions in presence of pitch rate
                obj=obj.determine_boundary_conditions_p(bc);
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                
                %alpha=alpha+thetadot_vec(i)*obj.t_step;
                % display current timestep
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/roll_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                % compute modal coupling coefficients if required
                
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                             %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                       %     gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end

            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        %% this function performs the UVLM computation for a pitching oscillation around p_ref
        
        function  obj=f_compute_timestep(obj,omega)

            %ts1=2*pi/(omega*obj.settings.spp);
%             if omega<10
%                 ts1=2*pi/(omega*obj.settings.spp*2);
%                 ts2=ts1;
%             else
            %ts2=obj.reference.c_ref/(norm(obj.Uinf)*16);
%             end
            %if ts1<ts2
              %  obj.t_step=ts1;
            %else
               
            %end
            
            %Method by simon
            %choose time constant to ensure similar area in wakepanels
            obj.t_step=diff(obj.grid(1,obj.panels(2:3,1)))/norm(obj.Uinf);
            % if time step is too coarse - choose finer and display warning
            minTStep=30;
            if 2*pi/omega/obj.t_step < minTStep
                disp('WARNING frequency too high for such a coarse grid on the wing in streamwise direction')
                disp('timestep is being changed, results not so good try finer mesh to overcome issues')
                obj.t_step=2*pi/omega/minTStep;
            end
            n_steps_min=obj.settings.wakelength_factor*obj.reference.b_ref*2/(obj.t_step*norm(obj.Uinf));
            n_steps_norm=length(0:obj.t_step:obj.settings.n_osc*pi/omega);
            n_steps=max(n_steps_norm,n_steps_min);
            obj.settings.n_osc=round(obj.t_step*n_steps*omega/pi);
        end
        
        function obj=solve_unsteady_pitch(obj,amplitude,k,varargin)
            
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                shape=varargin{4};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end 
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            % compute pitching amplitude in [rad]
            pitching_amplitude=amplitude*pi/180;
            % initialize pitch rate for each timestep
            thetadot_vec=pitching_amplitude*omega*cos(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            
            alpha=obj.state.alpha*pi/180;
            fprintf(['UVLM Pitch Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            %% COMPUTATION
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['pitch_k=' num2str(k)]);
            end
            
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            % UVLM time marching computation
            for i=1:length(obj.t_vec)
                obj.Uinf=norm(obj.Uinf)*[cos(alpha*obj.Ma_corr) 0 sin(alpha*obj.Ma_corr)];
                % compute boundary conditions in presence of pitch rate
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(thetadot_vec(i)+thetadot_vec(i+1))/2;
                else
                    bc=thetadot_vec(i);
                end
%                 bc=thetadot_vec(i);
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,bc,0);
                end
                obj=obj.determine_boundary_conditions_q(bc);
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                %obj=obj.compute_influence_coeff_matrix_wake();
                obj.alpha_s(i)=alpha;
                

                alpha=alpha+bc*obj.t_step;

                
                % display current timestep
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess(bc);
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/pitch_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                % compute modal coupling coefficients if required
                
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                         %   gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_yaw(obj,amplitude,k,varargin)
            
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                shape=varargin{4};
                                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end
            
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            % compute pitching amplitude in [rad]
            yawing_amplitude=amplitude*pi/180;
            % initialize pitch rate for each timestep
            qsidot_vec=yawing_amplitude*omega*cos(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            
            alpha=obj.state.alpha*pi/180;
            fprintf(['UVLM Yaw Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            %% COMPUTATION
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['yaw_k=' num2str(k)]);
            end
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            beta=0;
            % UVLM time marching computation
            for i=1:length(obj.t_vec)
                
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(qsidot_vec(i)+qsidot_vec(i+1))/2;
                else
                    bc=qsidot_vec(i);
                end 
                
%                  bc=qsidot_vec(i);
                obj.Uinf=a2bf(norm(obj.Uinf),alpha,beta,obj.Ma_corr);
                
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,bc);
                end
                % compute boundary conditions in presence of pitch rate
                obj=obj.determine_boundary_conditions_r(bc);
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
               
                 beta=beta+bc*obj.t_step;


                % display current timestep
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/yaw_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                % compute modal coupling coefficients if required
                
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                            %gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        function obj=solve_unsteady_sideheave(obj,amplitude,k,varargin)
            
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                 shape=varargin{4};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end
            
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            heaving_amplitude=amplitude;
            
            vmax=heaving_amplitude*omega;
            
            u_vec=vmax*cos(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            u=obj.Uinf;
            fprintf(['UVLM Sideheave Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['sideheave_k=' num2str(k)]);
            end
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            for i=1:length(obj.t_vec)
                
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(u_vec(i)+u_vec(i+1))/2;
                else
                    bc=u_vec(i);
                end
%                  bc=u_vec(i);
                obj.Uinf=u+[0 -bc 0];
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                obj=obj.determine_boundary_conditions();
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/sideheave_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                         %   gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_heave(obj,amplitude,k,varargin)
            
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                shape=varargin{4};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end
            
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            heaving_amplitude=amplitude;
            
            vmax=heaving_amplitude*omega;
            
            u_vec=vmax*cos(omega*obj.t_vec);
            
           % u_vec=vmax*sin(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            u=obj.Uinf;
            fprintf(['UVLM Heave Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['heave_k=' num2str(k)]);
            end
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            for i=1:length(obj.t_vec)
                %method 1
%                 bc=(u_vec(i));
                
                %method orig
                if (i<length(obj.t_vec))
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(u_vec(i)+u_vec(i+1))/2;
                else
                    bc=(u_vec(i));
                end

                %method 2
%                 if i==1
%                     bc=(u_vec(i)+u_vec(i+1))/2;
%                 elseif (i>1) && (i<length(obj.t_vec))
%                     bc=(u_vec(i-1)+u_vec(i+1)+u_vec(i))/3;
%                 else
%                     bc=(u_vec(i)+u_vec(i-1))/2;
%                 end
  
                obj.Uinf=u+[0 0 -bc*obj.Ma_corr];
                
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                obj.alpha_s(i)=atan(obj.Uinf(3)/obj.Uinf(1));
                obj=obj.determine_boundary_conditions();
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();

                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
                
                %obj=obj.compute_influence_coeff_matrix();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/heave_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                          %  gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        
        function obj=solve_unsteady_foraft(obj,amplitude,k,varargin)
            %% PREPROCESSING
            add_args=varargin;
            obj=obj.initialize_unsteady_computation_settings(nargin,add_args);
            shape=0;
            % set variables for computation of mode coupling coefficients
            if obj.settings.modal_data==1
                aircraft=varargin{2};
                aircraft_structure=varargin{3};
                shape=varargin{4};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            end
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            foraft_amplitude=amplitude;
            
            vmax=foraft_amplitude*omega;
            
            u_vec=vmax*cos(omega*obj.t_vec);
            
            obj=obj.initialize_unsteady_computation();
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['foraft_k=' num2str(k)]);
            end
            u=obj.Uinf;
            fprintf(['UVLM ForAft Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            
            for i=1:length(obj.t_vec)
                
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    bc=(u_vec(i)+u_vec(i+1))/2;
                else
                    bc=u_vec(i);
                end
%                  bc=u_vec(i);
                 
                obj.Uinf=u+[-bc 0 0];
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                obj=obj.determine_boundary_conditions();
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/foraft_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                        %    gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*obj.reference.b_ref;
            obj.CM=obj.CM*obj.reference.c_ref;
            obj.CN=obj.CN*obj.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function [dx_panel,dy_panel,dz_panel]=compute_modal_deflectionfield(obj,n_mode,aircraft,aircraft_structure,shape)
            dx_panel=zeros(obj.settings.n_mode,length(obj.panels));
            dy_panel=zeros(obj.settings.n_mode,length(obj.panels));
            dz_panel=zeros(obj.settings.n_mode,length(obj.panels));

            for mod_ctr=1:n_mode
                
                aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,mod_ctr)*0.1+shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                
                dx_grid=obj.grid(1,:)-aircraft.grid_deflected(1,:)/obj.Ma_corr;
                dy_grid=obj.grid(2,:)-aircraft.grid_deflected(2,:);
                dz_grid=obj.grid(3,:)-aircraft.grid_deflected(3,:);
                %for unsteady force which acts on the middle of the  panel, NOT the middle of the vortex ring
                for xi=1:length(obj.panels)
                    dx_panel(mod_ctr,xi)=(dx_grid(obj.panels(1,xi))+dx_grid(obj.panels(2,xi))+dx_grid(obj.panels(3,xi))+dx_grid(obj.panels(4,xi)))/4*10;
                    dy_panel(mod_ctr,xi)=(dy_grid(obj.panels(1,xi))+dy_grid(obj.panels(2,xi))+dy_grid(obj.panels(3,xi))+dy_grid(obj.panels(4,xi)))/4*10;
                    dz_panel(mod_ctr,xi)=(dz_grid(obj.panels(1,xi))+dz_grid(obj.panels(2,xi))+dz_grid(obj.panels(3,xi))+dz_grid(obj.panels(4,xi)))/4*10;
                end
            end
        end
        
        function [dx_panel,dy_panel,dz_panel]=compute_modal_deflectionfield_c4(obj,n_mode,aircraft,aircraft_structure,shape)
            dx_panel=zeros(obj.settings.n_mode,length(obj.panels));
            dy_panel=zeros(obj.settings.n_mode,length(obj.panels));
            dz_panel=zeros(obj.settings.n_mode,length(obj.panels));

            for mod_ctr=1:n_mode
                
                aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,mod_ctr)*0.1+shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                
                dx_grid=obj.grid(1,:)-aircraft.grid_deflected(1,:)/obj.Ma_corr;
                dy_grid=obj.grid(2,:)-aircraft.grid_deflected(2,:);
                dz_grid=obj.grid(3,:)-aircraft.grid_deflected(3,:);
                
                %for steady force which acts on the one quarter chord of
                %the panel
                for xi=1:length(obj.panels)
                    dx_panel(mod_ctr,xi)=(0.375*dx_grid(obj.panels(1,xi))+0.375*dx_grid(obj.panels(2,xi))+0.125*dx_grid(obj.panels(3,xi))+0.125*dx_grid(obj.panels(4,xi)))*10;
                    dy_panel(mod_ctr,xi)=(0.375*dy_grid(obj.panels(1,xi))+0.375*dy_grid(obj.panels(2,xi))+0.125*dy_grid(obj.panels(3,xi))+0.125*dy_grid(obj.panels(4,xi)))*10;
                    dz_panel(mod_ctr,xi)=(0.375*dz_grid(obj.panels(1,xi))+0.375*dz_grid(obj.panels(2,xi))+0.125*dz_grid(obj.panels(3,xi))+0.125*dz_grid(obj.panels(4,xi)))*10;
                end
            end
        end

        function obj=solve_unsteady_mode(obj,amplitude,k,nmod,aircraft,aircraft_structure,mode,varargin)
            %% PREPROCESSING
            shape=0;
            if ~isempty(varargin)
                add_args={nmod,aircraft,aircraft_structure};
                shape=varargin{1};
                aircraft_structure.nodal_deflections=shape;
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
            else
                add_args={nmod,aircraft,aircraft_structure};
            end
            obj=obj.initialize_unsteady_computation_settings(nargin-1,add_args);
            
            % compute actual frequency from reduced frequency TODO: check
            % formula blow omega/2*pi/n
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            obj=obj.initialize_unsteady_computation();
            fprintf(['UVLM Mode Oscillation computing for mode ' num2str(mode) ' k=' num2str(k) '\n']);
            fprintf('processed:     ');
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name '/mode',num2str(mode),'_k=' num2str(k)]);
            end
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end
            
            %calculate shape Deflection of grid
            aircraft_structure.nodal_deflections=shape;                
            aircraft_structure=aircraft_structure.f_postprocess();
            aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            shapeDeflection=aircraft.grid_deflected-aircraft.grid;
            
            %calculate positive Deflection grid
            aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,mode)*amplitude*1;                
            aircraft_structure=aircraft_structure.f_postprocess();
            aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            positiveModeDeflection=aircraft.grid_deflected-aircraft.grid;
            
            
            for i=1:length(obj.t_vec)
                prv_colloc=obj.colloc;
                %calculate current deflection
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    aircraft.grid_deflected=aircraft.grid+shapeDeflection+sin(omega*(obj.t_vec(i)+obj.t_vec(i+1))/2)*positiveModeDeflection;
                else
                    aircraft.grid_deflected=aircraft.grid+shapeDeflection+sin(omega*obj.t_vec(i))*positiveModeDeflection;
                end

                
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                obj=obj.update_grid();
                
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                
                obj=obj.determine_boundary_conditions_deflection(prv_colloc,obj.t_step);
                
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                % postprocess timestep, compute aerodynamic coefficients from solution
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/mode',num2str(mode),'_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                            %gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            obj.CL=obj.CL*aircraft.reference.b_ref;
            obj.CM=obj.CM*aircraft.reference.c_ref;
            obj.CN=obj.CN*aircraft.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_controlmode(obj,amplitude,k,nmod,aircraft,aircraft_structure,cs_name,varargin)
            shape=0;
            if ~isempty(varargin)
                add_args={nmod,aircraft,aircraft_structure};
                ref_def=varargin{1};
                
                shape=varargin{2};
                aircraft_structure.nodal_deflections=shape;
                if obj.settings.modal_data==1
                    aircraft_structure=aircraft_structure.f_postprocess();
                    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                    obj.grid(2,:)=aircraft.grid_deflected(2,:);
                    obj.grid(3,:)=aircraft.grid_deflected(3,:);
                    obj=obj.update_grid();
                end
                
            else
                add_args={nmod,aircraft};
            end
            obj=obj.initialize_unsteady_computation_settings(nargin-1,add_args);
            
            % compute actual frequency from reduced frequency TODO: check
            % formula blow omega/2*pi/n
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            obj=obj.initialize_unsteady_computation();
            fprintf(['UVLM ControlMode Oscillation computing for control surface ' cs_name ' k=' num2str(k) '\n']);
            fprintf('processed:     ');
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name '/control_mode_',cs_name,'_k=' num2str(k)]);
            end
            
      
            if obj.settings.modal_data==1
                [dx_panel,dy_panel,dz_panel]=obj.compute_modal_deflectionfield(obj.settings.n_mode,aircraft,aircraft_structure,shape);
                [dx_panel_c4,dy_panel_c4,dz_panel_c4]=obj.compute_modal_deflectionfield_c4(obj.settings.n_mode,aircraft,aircraft_structure,shape);
            end

            aircraft=aircraft.f_set_control_surface(cs_name,ref_def);
            aircraft=aircraft.compute_grid();
            grid=aircraft.grid;
            aircraft=aircraft.f_set_control_surface(cs_name,amplitude+ref_def);
            aircraft=aircraft.compute_grid();
            grid_p=aircraft.grid;
            aircraft=aircraft.f_set_control_surface(cs_name,-amplitude+ref_def);
            aircraft=aircraft.compute_grid();
            grid_n=aircraft.grid;
            
            for i=1:length(obj.t_vec)
                prv_colloc=obj.colloc;
%                 fact=sin(omega*obj.t_vec(i));
                
                if i<length(obj.t_vec)
                    % Not sure if this average is right, but it showed better result
                    % for heave and pitch motion compared to the Theodorsen
                    % method
                    aircraft.grid=aircraft.grid+(grid_p-grid)*sin(omega*(obj.t_vec(i)+obj.t_vec(i+1))/2);
                else
                    aircraft.grid=aircraft.grid+(grid_p-grid)*sin(omega*obj.t_vec(i));
                end
%                 if fact>0
%                     aircraft.grid=grid+(grid_p-grid)*fact;
%                 elseif fact==0;
%                     aircraft.grid=grid;
%                 elseif fact<0
%                     aircraft.grid=grid+(grid_n-grid)*abs(fact);
%                 end
                    
                if shape~=0
                    aircraft_structure.nodal_deflections=shape;
                    aircraft_structure=aircraft_structure.f_postprocess();
                    aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                    obj.grid(2,:)=aircraft.grid_deflected(2,:);
                    obj.grid(3,:)=aircraft.grid_deflected(3,:);
                    obj=obj.update_grid();
                else
                  %  aircraft=aircraft.compute_grid();
                    obj.grid(1,:)=aircraft.grid(1,:)/obj.Ma_corr;
                    obj.grid(2,:)=aircraft.grid(2,:);
                    obj.grid(3,:)=aircraft.grid(3,:);
                    obj=obj.update_grid();
                end

                
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                obj=obj.determine_boundary_conditions_deflection(prv_colloc,obj.t_step);
                
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                obj=obj.shed_wake();
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1  
                    obj.write_tecplot_wake(['results/' obj.case_name '/control_mode_',cs_name,'_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    for mod_ctr=1:obj.settings.n_mode
                        gaf=0;
                        for xi=1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*dz_panel_c4(mod_ctr,xi)+obj.F_body2(2,xi)*dy_panel_c4(mod_ctr,xi)+obj.F_body2(1,xi)*dx_panel_c4(mod_ctr,xi)+...       
                            +obj.F_body_unsteady(3,xi)*dz_panel(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*dy_panel(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*dx_panel(mod_ctr,xi);
                        %    gaf=gaf+obj.colloc_nvec(3,xi).*obj.cp(xi)*dz_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi).*obj.cp(xi)*dy_panel(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi).*obj.cp(xi)*dx_panel(mod_ctr,xi)*obj.area(xi);
                        end
                        obj.C_modes(mod_ctr,i)=gaf/(obj.reference.S_ref*obj.qinf);
                    end
                end
            end
            
            aircraft=aircraft.f_set_control_surface(cs_name,0);
            obj.CL=obj.CL*aircraft.reference.b_ref;
            obj.CM=obj.CM*aircraft.reference.c_ref;
            obj.CN=obj.CN*aircraft.reference.b_ref;
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_mode_mesh(obj,amplitude,k,mod,nmod,BCs,X,Y,Z,eigm)
            %% PREPROCESSING
            %obj.settings.modal_data=0;
            obj.settings.n_osc=15;
            obj.settings.spp=45;
            obj.settings.modal_data=1;
            
            obj.settings.n_mode=nmod;
            % compute actual frequency from reduced frequency TODO: check
            % formula blow omega/2*pi/n
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            
            obj=obj.initialize_unsteady_computation();
            fprintf(['UVLM Mode Oscillation computing for mode ' num2str(mod) ' k=' num2str(k) '\n']);
            fprintf('processed:     ');
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name '/mode',num2str(mod),'_k=' num2str(k)]);
            end
            
            [r,II]=sortrows([X;Y;Z]',2);
            
            for i=1:length(obj.panels(1,:))
                Eigm_x0(i)=0.375*obj.grid(1,obj.panels(1,i))+0.375*obj.grid(1,obj.panels(2,i))+0.125*obj.grid(1,obj.panels(3,i))+0.125*obj.grid(1,obj.panels(4,i));
                Eigm_y0(i)=0.375*obj.grid(2,obj.panels(1,i))+0.375*obj.grid(2,obj.panels(2,i))+0.125*obj.grid(2,obj.panels(3,i))+0.125*obj.grid(2,obj.panels(4,i));
                Eigm_z0(i)=0.375*obj.grid(3,obj.panels(1,i))+0.375*obj.grid(3,obj.panels(2,i))+0.125*obj.grid(3,obj.panels(3,i))+0.125*obj.grid(3,obj.panels(4,i));
                Eigm_x0_center(i)=0.25*obj.grid(1,obj.panels(1,i))+0.25*obj.grid(1,obj.panels(2,i))+0.25*obj.grid(1,obj.panels(3,i))+0.25*obj.grid(1,obj.panels(4,i));
                Eigm_y0_center(i)=0.25*obj.grid(2,obj.panels(1,i))+0.25*obj.grid(2,obj.panels(2,i))+0.25*obj.grid(2,obj.panels(3,i))+0.25*obj.grid(2,obj.panels(4,i));
                Eigm_z0_center(i)=0.25*obj.grid(3,obj.panels(1,i))+0.25*obj.grid(3,obj.panels(2,i))+0.25*obj.grid(3,obj.panels(3,i))+0.25*obj.grid(3,obj.panels(4,i));
            end
            for mod_ctr=1:5
                ck=1;
                jj=1;
                for j=1:length(X)
                    if BCs(ck)~=j
                        XX(j)=X(j)+eigm(jj,mod_ctr)*amplitude;
                        YY(j)=Y(j)+eigm(jj+1,mod_ctr)*amplitude;
                        ZZ(j)=Z(j)+eigm(jj+2,mod_ctr)*amplitude;
                        jj=jj+5;
                    else
                        XX(j)=X(j);
                        YY(j)=Y(j);
                        ZZ(j)=Z(j);
                        ck=ck+1;
                    end
                end
                coords_AGARD=[XX;YY;ZZ];
                grid_deflected=coords_AGARD(:,II);
                grid_deflected=[grid_deflected [grid_deflected(1,:);-grid_deflected(2,:);grid_deflected(3,:)]];
                
                obj.grid(1,:)=grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=grid_deflected(2,:);
                obj.grid(3,:)=grid_deflected(3,:);
                obj=obj.update_grid();
                for i=1:length(obj.panels(1,:))
                    Eigm_x(mod_ctr,i)=0.375*obj.grid(1,obj.panels(1,i))+0.375*obj.grid(1,obj.panels(2,i))+0.125*obj.grid(1,obj.panels(3,i))+0.125*obj.grid(1,obj.panels(4,i))-Eigm_x0(i);
                    Eigm_y(mod_ctr,i)=0.375*obj.grid(2,obj.panels(1,i))+0.375*obj.grid(2,obj.panels(2,i))+0.125*obj.grid(2,obj.panels(3,i))+0.125*obj.grid(2,obj.panels(4,i))-Eigm_y0(i);
                    Eigm_z(mod_ctr,i)=0.375*obj.grid(3,obj.panels(1,i))+0.375*obj.grid(3,obj.panels(2,i))+0.125*obj.grid(3,obj.panels(3,i))+0.125*obj.grid(3,obj.panels(4,i))-Eigm_z0(i);
                    
                    Eigm_Theta(mod_ctr,i)= atan(((0.5*obj.grid(3,obj.panels(1,i))+0.5*obj.grid(3,obj.panels(2,i)))+(0.5*obj.grid(3,obj.panels(3,i))+0.5*obj.grid(3,obj.panels(4,i))))/((0.5*obj.grid(1,obj.panels(4,i))+0.5*obj.grid(1,obj.panels(3,i)))-(0.5*obj.grid(1,obj.panels(1,i))+0.5*obj.grid(1,obj.panels(2,i)))));
                end
                for i=1:length(obj.panels(1,:))
                    Eigm_x_center(mod_ctr,i)=0.25*obj.grid(1,obj.panels(1,i))+0.25*obj.grid(1,obj.panels(2,i))+0.25*obj.grid(1,obj.panels(3,i))+0.25*obj.grid(1,obj.panels(4,i))-Eigm_x0_center(i);
                    Eigm_y_center(mod_ctr,i)=0.25*obj.grid(2,obj.panels(1,i))+0.25*obj.grid(2,obj.panels(2,i))+0.25*obj.grid(2,obj.panels(3,i))+0.25*obj.grid(2,obj.panels(4,i))-Eigm_y0_center(i);
                    Eigm_z_center(mod_ctr,i)=0.25*obj.grid(3,obj.panels(1,i))+0.25*obj.grid(3,obj.panels(2,i))+0.25*obj.grid(3,obj.panels(3,i))+0.25*obj.grid(3,obj.panels(4,i))-Eigm_z0_center(i);
                    Eigm_Theta_center(mod_ctr,i)= atan(((0.5*obj.grid(3,obj.panels(1,i))+0.5*obj.grid(3,obj.panels(2,i)))+(0.5*obj.grid(3,obj.panels(3,i))+0.5*obj.grid(3,obj.panels(4,i))))/((0.5*obj.grid(1,obj.panels(4,i))+0.5*obj.grid(1,obj.panels(3,i)))-(0.5*obj.grid(1,obj.panels(1,i))+0.5*obj.grid(1,obj.panels(2,i)))));
                end
            end
          %  for i=1:5
%     figure
%     hold on
%     plot3(Eigm_x0,Eigm_y0,Eigm_z0,'bx');
%     k=1;
%     jj=1;
%     for j=1:length(Eigm_x0)
%         if BCs(k)~=j
%             plot3(Eigm_x0(j)+amplitude*Eigm_x(3,j),Eigm_y0(j)+amplitude*Eigm_y(3,j),Eigm_z0(j)+amplitude*Eigm_z(3,j),'ro');   
%         else
%             k=k+1;
%         end
%     end  
%                 fullstructure.nodal_deflections=fullstructure.modeshapes(:,(i*2))*-0.5;
%             fullstructure=fullstructure.f_postprocess();
%             aircraft=aircraft.compute_deflected_grid(fullstructure.f_get_deflections);
%             aircraft.plot_grid_deflected;
%end
            for i=1:length(obj.t_vec)
                prv_colloc=obj.colloc;
                
                ck=1;
                jj=1;
                for j=1:length(X)
                    if BCs(ck)~=j
                        XX(j)=X(j)+eigm(jj,mod)*amplitude*sin(omega*obj.t_vec(i));
                        YY(j)=Y(j)+eigm(jj+1,mod)*amplitude*sin(omega*obj.t_vec(i));
                        ZZ(j)=Z(j)+eigm(jj+2,mod)*amplitude*sin(omega*obj.t_vec(i));
                        jj=jj+5;
                    else
                        XX(j)=X(j);
                        YY(j)=Y(j);
                        ZZ(j)=Z(j);
                        ck=ck+1;
                    end
                end
                coords_AGARD=[XX;YY;ZZ];
                grid_deflected=coords_AGARD(:,II);
                grid_deflected=[grid_deflected [grid_deflected(1,:);-grid_deflected(2,:);grid_deflected(3,:)]];
                
                obj.grid(1,:)=grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=grid_deflected(2,:);
                obj.grid(3,:)=grid_deflected(3,:);
                obj=obj.update_grid();
                
                %obj=obj.compute_influence_coeff_matrix();
                %obj.grid_wake(:,1:obj.row_length)=aircraft.grid_wake(:,1:obj.row_length);
                
                obj=obj.determine_boundary_conditions_deflection(prv_colloc,obj.t_step);
                
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                %obj=obj.shed_wake();
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                %                 if obj.settings.movie==1
                %                     close all
                %                     hfig=figure(33);
                %                     set(hfig, 'Position', [0 0 1024*1.5 768*1.5])
                %                 end
                % postprocess timestep, compute aerodynamic coefficients
                % from solution
                obj=obj.f_postprocess();
               % obj=obj.f_postprocess_exact();
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/mode',num2str(mod),'_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    gaf=0;
                    for mod_ctr=1:5
                        gaf=0;
                        for xi=length(obj.cp)/2:1:length(obj.cp)
                            %correction simon: simulations of simple
                            %heave/pitch motion show that CM only agrees
                            %when the unsteady force is also applied to the
                            %one quarter point - so here it should be also
                            %the one quarter point which is used to
                            %determine the virtual work
                            % update: rafael checked and it seems that mid of panel is correct fvap for unsteady part
                            gaf=gaf+obj.F_body2(3,xi)*Eigm_z(mod_ctr,xi)+obj.F_body2(2,xi)*Eigm_y(mod_ctr,xi)+obj.F_body2(1,xi)*Eigm_x(mod_ctr,xi)+...
                            +obj.F_body_unsteady(3,xi)*Eigm_z_center(mod_ctr,xi)+obj.F_body_unsteady(2,xi)*Eigm_y_center(mod_ctr,xi)+obj.F_body_unsteady(1,xi)*Eigm_x_center(mod_ctr,xi);%...
                       %-Eigm_Theta(mod_ctr,xi)*Eigm_z_center(mod_ctr,xi)*(Eigm_x(mod_ctr,xi)-Eigm_x_center(mod_ctr,xi));
                           % gaf=gaf+obj.colloc_nvec(3,xi)*obj.cp(xi)*Eigm_z(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi)*obj.cp(xi)*Eigm_y(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi)*obj.cp(xi)*Eigm_x(mod_ctr,xi)*obj.area(xi);%*cos(Eigm_Theta(mod_ctr,xi));
                        end
                        obj.C_modes(mod_ctr,i)=gaf/obj.qinf;
                    end
                end
            end
            obj.CX=obj.CX/abs(amplitude);
            obj.CY=obj.CY/abs(amplitude);
            obj.CZ=obj.CZ/abs(amplitude);
            obj.CL=obj.CL/abs(amplitude);
            obj.CM=obj.CM/abs(amplitude);
            obj.CN=obj.CN/abs(amplitude);
            obj.C_modes=obj.C_modes./abs(amplitude^2);
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=solve_unsteady_gust_mode_mesh(obj,amplitude,k,nmod,BCs,X,Y,Z,eigm)
            obj.settings.n_osc=20;
            obj.settings.spp=45;
            obj.settings.modal_data=1;
            
            obj.settings.n_mode=nmod;
            % compute actual frequency from reduced frequency
            omega=k*2*norm(obj.Uinf)/obj.reference.c_ref;
            % compute required timestep
            obj=obj.f_compute_timestep(omega);
            % initialize time marching vector
            obj.t_vec=0:obj.t_step:obj.settings.n_osc*pi/omega;
            % compute pitching amplitude in [rad]
            gust_amplitude=amplitude;
            % initialize pitch rate for each timestep
            obj=obj.initialize_unsteady_computation();
            
            fprintf(['UVLM Gust Mode Oscillation computing for: k=' num2str(k) '\n']);
            fprintf('processed:     ');
            
            if obj.settings.movie==1
                mkdir(['results/' obj.case_name],['gust_k=' num2str(k)]);
            end
            
              
            [r,II]=sortrows([X;Y;Z]',2);
            
            for i=1:length(obj.panels(1,:))
                Eigm_x0(i)=0.375*obj.grid(1,obj.panels(1,i))+0.375*obj.grid(1,obj.panels(2,i))+0.125*obj.grid(1,obj.panels(3,i))+0.125*obj.grid(1,obj.panels(4,i));
                Eigm_y0(i)=0.375*obj.grid(2,obj.panels(1,i))+0.375*obj.grid(2,obj.panels(2,i))+0.125*obj.grid(2,obj.panels(3,i))+0.125*obj.grid(2,obj.panels(4,i));
                Eigm_z0(i)=0.375*obj.grid(3,obj.panels(1,i))+0.375*obj.grid(3,obj.panels(2,i))+0.125*obj.grid(3,obj.panels(3,i))+0.125*obj.grid(3,obj.panels(4,i));
            end
            for mod_ctr=1:5
                ck=1;
                jj=1;
                for j=1:length(X)
                    if BCs(ck)~=j
                        XX(j)=X(j)+eigm(jj,mod_ctr);
                        YY(j)=Y(j)+eigm(jj+1,mod_ctr);
                        ZZ(j)=Z(j)+eigm(jj+2,mod_ctr);
                        jj=jj+5;
                    else
                        XX(j)=X(j);
                        YY(j)=Y(j);
                        ZZ(j)=Z(j);
                        ck=ck+1;
                    end
                end
                coords_AGARD=[XX;YY;ZZ];
                grid_deflected=coords_AGARD(:,II);
                grid_deflected=[grid_deflected [grid_deflected(1,:);-grid_deflected(2,:);grid_deflected(3,:)]];
                
                for i=1:length(obj.panels(1,:))
                    Eigm_x(mod_ctr,i)=0.375*grid_deflected(1,obj.panels(1,i))+0.375*grid_deflected(1,obj.panels(2,i))+0.125*grid_deflected(1,obj.panels(3,i))+0.125*grid_deflected(1,obj.panels(4,i))-Eigm_x0(i);
                    Eigm_y(mod_ctr,i)=0.375*grid_deflected(2,obj.panels(1,i))+0.375*grid_deflected(2,obj.panels(2,i))+0.125*grid_deflected(2,obj.panels(3,i))+0.125*grid_deflected(2,obj.panels(4,i))-Eigm_y0(i);
                    Eigm_z(mod_ctr,i)=0.375*grid_deflected(3,obj.panels(1,i))+0.375*grid_deflected(3,obj.panels(2,i))+0.125*grid_deflected(3,obj.panels(3,i))+0.125*grid_deflected(3,obj.panels(4,i))-Eigm_z0(i);
                end
            end
            
            x_min=min(obj.grid(1,:));
            x_max=max(obj.grid_wake(1,:));
            
            t_period=2*pi/omega;
            obj.x_kin=x_min:obj.t_step*obj.Uinf(1):x_max;
            t_kin=obj.x_kin/obj.Uinf(1);
            gust_shape=gust_amplitude*sin(omega*t_kin);
            
            obj.Uinf_x=zeros(3,length(obj.x_kin));
            obj.Uinf_x(1,:)=obj.Uinf(1);
            obj.Uinf_x(2,:)=obj.Uinf(2);
            obj.Uinf_x(3,:)=obj.Uinf(3)+gust_shape;
            % figure
            for i=1:length(obj.t_vec)
                gust_shape=gust_amplitude*sin(omega*t_kin-2*pi/t_period*obj.t_step*i);
                
                obj.Uinf_x(3,:)=obj.Uinf(3)+gust_shape;
                % compute boundary conditions in presence of gust mode
                obj=obj.determine_boundary_conditions_gust();
                % compute new vorticity of wake
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                % add wake effect
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                % save current timestep solution
                obj.Gamma_prv=obj.Gamma;
                % compute new solution
                obj.Gamma=linsolve(obj.Abb,bx');
                % shed wake
                obj=obj.shed_wake_gust();
                % display current timestep
                fprintf('\b\b\b\b %03d',round(obj.t_vec(i)/obj.t_vec(end)*100));
                % postprocess timestep, compute aerodynamic coefficients
                obj=obj.f_postprocess();
                
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/gust_k=' num2str(k) '/time_series' num2str(i)],1);
                end
                
                % compute modal coupling coefficients if required
                if obj.settings.modal_data==1
                    gaf=0;
                    for mod_ctr=1:5
                        gaf=0;
                        for xi=length(obj.cp)/2:1:length(obj.cp)
                            gaf=gaf+obj.F_body2(3,xi)*Eigm_z(mod_ctr,xi)+obj.F_body2(2,xi)*Eigm_y(mod_ctr,xi)+obj.F_body2(1,xi)*Eigm_x(mod_ctr,xi);%*cos(Eigm_Theta(mod_ctr,xi));
                           % gaf=gaf+obj.colloc_nvec(3,xi)*obj.cp(xi)*Eigm_z(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(2,xi)*obj.cp(xi)*Eigm_y(mod_ctr,xi)*obj.area(xi)+obj.colloc_nvec(1,xi)*obj.cp(xi)*Eigm_x(mod_ctr,xi)*obj.area(xi);%*cos(Eigm_Theta(mod_ctr,xi));
                        end
                        obj.C_modes(mod_ctr,i)=gaf/obj.qinf;
                    end
                end
            end
            obj.CX=obj.CX/abs(gust_amplitude);
            obj.CY=obj.CY/abs(gust_amplitude);
            obj.CZ=obj.CZ/abs(gust_amplitude);
            obj.CL=obj.CL/abs(gust_amplitude)*obj.reference.b_ref;
            obj.CM=obj.CM/abs(gust_amplitude)*obj.reference.c_ref;
            obj.CN=obj.CN/abs(gust_amplitude)*obj.reference.b_ref;
            obj.C_modes=obj.C_modes./abs(gust_amplitude);
            obj=obj.complex_coeffs_from_timeseries(omega);
        end
        
        function obj=f_postprocess(obj,varargin)
            
            U=norm(obj.Uinf')-(obj.Abw*obj.Gamma_wake)'*0;
            if nargin==2
                qq=varargin{1};
            else
                qq=0;
            end
           
            if obj.Ma_corr<1
                beta_inf=obj.Ma_corr;
            else
                beta_inf=1;
            end
           
            obj.F_body2=zeros(3,length(obj.Gamma));
            obj.F_body_unsteady=zeros(3,length(obj.Gamma));
            obj.cp_st=zeros(1,length(obj.Gamma));
            obj.cp_un=zeros(1,length(obj.Gamma));
           
            obj.v_ind=+[(obj.Abw_x*obj.Gamma_wake)';
                (obj.Abw_y*obj.Gamma_wake)';
                (obj.Abw_z*obj.Gamma_wake)']+[(obj.Abb_x*obj.Gamma)';
                (obj.Abb_y*obj.Gamma)';
               (obj.Abb_z*obj.Gamma)'];
% %             
%               obj.v_ind=+[(obj.Abw_x*obj.Gamma_wake)';
%                 (obj.Abw_y*obj.Gamma_wake)';
%                 (obj.Abw_z*obj.Gamma_wake)']+[(obj.Awb_x*obj.Gamma)';
%                 (obj.Awb_y*obj.Gamma)';
%                (obj.Awb_z*obj.Gamma)'];

            dq=qq*obj.Ma_corr;
            for i=1:1:length(obj.colloc)
                rrr=obj.r(:,i);
                rot=[0 dq 0];
                dV=[rrr(2)*rot(3)-rrr(3)*rot(2),rrr(3)*rot(1)-rrr(1)*rot(3),rrr(1)*rot(2)-rrr(2)*rot(1)];   %r1xr2=cross(r1,r2);
                obj.v_ind(:,i)=obj.v_ind(:,i)+dV';
            end
            Uinf=obj.Uinf;
            
            sci=(obj.grid(:,obj.panels(4,:))-obj.grid(:,obj.panels(1,:)));
            sco=(obj.grid(:,obj.panels(3,:))-obj.grid(:,obj.panels(2,:)));
            for i=1:length(obj.panels)
                sc=0.5*norm(sci(:,i))+0.5*norm(sco(:,i));
                
                
                %  ss=0.5*norm(obj.grid(:,obj.panels(2,i))-obj.grid(:,obj.panels(1,i)))+0.5*norm(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(3,i)));
                pi=0.75*obj.grid(:,obj.panels(1,i))+0.25*obj.grid(:,obj.panels(4,i));
                po=0.75*obj.grid(:,obj.panels(2,i))+0.25*obj.grid(:,obj.panels(3,i));
                sp=po-pi;
                %             if nargin==2 %related to linearization
                %                 sc=0.5*norm(obj.grid_init(:,obj.panels(4,i))-obj.grid_init(:,obj.panels(1,i)))+0.5*norm(obj.grid_init(:,obj.panels(3,i))-obj.grid_init(:,obj.panels(2,i)));
                %                 ss=0.5*norm(obj.grid_init(:,obj.panels(2,i))-obj.grid_init(:,obj.panels(1,i)))+0.5*norm(obj.grid_init(:,obj.panels(4,i))-obj.grid_init(:,obj.panels(3,i)));
                %             else
                %                 sc=0.5*norm(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i)))+0.5*norm(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(2,i)));
                %                 ss=0.5*norm(obj.grid(:,obj.panels(2,i))-obj.grid(:,obj.panels(1,i)))+0.5*norm(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(3,i)));
                %             end
                % new average method by rafael
                if i==1
                    isLE=1;
                else
                    isLE=obj.is_te(i-1);
                end
                vel=Uinf'+obj.v_ind(:,i);
                if  isLE      %leading edge panels
                    obj.cp_st(i)=((obj.Gamma(i))*U(i)/sc)/(0.5*norm(Uinf')^2);
                    
                    spxvel=[sp(2)*vel(3)-sp(3)*vel(2),sp(3)*vel(1)-sp(1)*vel(3),sp(1)*vel(2)-sp(2)*vel(1)];   %r1xr2=cross(r1,r2);
                    obj.F_body2(:,i)=obj.rho*(obj.Gamma(i))*spxvel*1/beta_inf;
                    if ~isempty(obj.Gamma_prv)
                        obj.cp_un(i)=((obj.Gamma(i)-obj.Gamma_prv(i))/2)/(0.5*norm(Uinf')^2*obj.t_step);
                        obj.F_body_unsteady(:,i)=obj.rho*((obj.Gamma(i)-obj.Gamma_prv(i))/2)*obj.area(i).*obj.colloc_nvec(:,i)/(obj.t_step); 
                    end
                else                    %all other bound panels
                    obj.cp_st(i)=((obj.Gamma(i)-obj.Gamma(i-1))*U(i)/sc)/(0.5*norm(Uinf')^2);
                    spxvel=[sp(2)*vel(3)-sp(3)*vel(2),sp(3)*vel(1)-sp(1)*vel(3),sp(1)*vel(2)-sp(2)*vel(1)];   %r1xr2=cross(r1,r2);
                    obj.F_body2(:,i)=obj.rho*(obj.Gamma(i)-obj.Gamma(i-1))*spxvel*1/beta_inf;

                    if ~isempty(obj.Gamma_prv)
                        obj.cp_un(i)=((obj.Gamma(i)+obj.Gamma(i-1))/2-(obj.Gamma_prv(i)+obj.Gamma_prv(i-1))/2)/(0.5*norm(Uinf')^2*obj.t_step);
                        obj.F_body_unsteady(:,i)=obj.rho*((obj.Gamma(i)+obj.Gamma(i-1))/2-(obj.Gamma_prv(i)+obj.Gamma_prv(i-1))/2)*obj.area(i).*obj.colloc_nvec(:,i)/(obj.t_step); 
                    end
                end
                
                    
                
            end
            
            if ~isempty(obj.cp_un)
                obj.cp=1/beta_inf*(obj.cp_st)+obj.cp_un;
             %   obj.F_body2=obj.F_body2+obj.F_body_unsteady;
            else
                obj.cp=1/beta_inf*obj.cp_st;
            end
            
            alpha=atan(Uinf(3)/Uinf(1));
            beta=atan(Uinf(2)/Uinf(1));
            
            obj.F_body=(obj.F_body2+obj.F_body_unsteady);
            obj.M_body=cross((obj.r_mid-obj.r),obj.F_body_unsteady); %this seems wrong! is M_body used somewhere??
            
%            obj.F_body(1,:)=obj.qinf*obj.area.*obj.colloc_nvec(1,:).*obj.cp;
%            obj.F_body(2,:)=obj.qinf*obj.area.*obj.colloc_nvec(2,:).*obj.cp;
%            obj.F_body(3,:)=obj.qinf*obj.area.*obj.colloc_nvec(3,:).*obj.cp;
            
            %                 %alpha2=atan((obj.Uinf(3)+uvwind(3))/(obj.Uinf(1)+uvwind(1)));
            %                 %beta=atan((obj.Uinf(2)+uvwind(2))/(obj.Uinf(1)+uvwind(1)));
            %
            %                 alpha=atan((obj.Uinf(3))/(obj.Uinf(1)));
            %                 beta=atan((obj.Uinf(2))/(obj.Uinf(1)));
            %                 %                alphai(i)=alpha-alpha2;
            %
            M_BA=[  cos(alpha)*cos(beta) , -cos(alpha)*sin(beta), -sin(alpha);
                sin(beta)             , cos(beta)            ,    0;
                sin(alpha)*cos(beta) , -sin(alpha)*sin(beta), cos(alpha);];
            
                        for i=1:length(obj.panels)
                            obj.F_aero(:,i)=M_BA'*obj.F_body(:,i);
                        end
            CX=sum(obj.F_body(1,:))/(obj.qinf*obj.reference.S_ref);
            CY=sum(obj.F_body(2,:))/(obj.qinf*obj.reference.S_ref);
            CZ=sum(obj.F_body(3,:))/(obj.qinf*obj.reference.S_ref);
            Cl=sum(obj.F_aero(3,:))/(obj.qinf*obj.reference.S_ref);
            Cdi=sum(obj.F_aero(1,:))/(obj.qinf*obj.reference.S_ref);
            Cdi2=sum(obj.F_body(1,:))/(obj.qinf*obj.reference.S_ref); %Cdi body fixed -> unsteady part has no influence since force always normal to panel this is only used for garrick comparison
           % Cdi=sum(sqrt(obj.F_drag(1,:).^2+obj.F_drag(2,:).^2+obj.F_drag(3,:).^2))/(obj.qinf*obj.reference.S_ref);
            %             Cl=sum(obj.F_aero(3,:))/(obj.qinf*obj.reference.S_ref);
            %             Cy=sum(obj.F_aero(2,:))/(obj.qinf*obj.reference.S_ref);
         
          %  obj.r(1,:)=obj.r(1,:)*obj.Ma_corr;
            CL=sum(+obj.F_body(3,:).*obj.r(2,:)-obj.F_body(2,:).*obj.r(3,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
            CM=(sum(-obj.F_body2(3,:).*obj.r(1,:)+obj.F_body2(1,:).*obj.r(3,:))+sum(-obj.F_body_unsteady(3,:).*obj.r_mid(1,:)+obj.F_body_unsteady(1,:).*obj.r_mid(3,:)))/(obj.qinf*obj.reference.S_ref*obj.reference.c_ref);
            CN=sum(obj.F_body(2,:).*obj.r(1,:)-obj.F_body(1,:).*obj.r(2,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
         %   obj.r(1,:)=obj.r(1,:)/obj.Ma_corr;
            obj.CX=[obj.CX CX];
            obj.CY=[obj.CY CY];
            obj.CZ=[obj.CZ CZ];
            obj.Cl=[obj.Cl Cl];
            
            obj.Cdi=[obj.Cdi Cdi];
            
            obj.Cdi2=[obj.Cdi2 Cdi2];
            obj.CL=[obj.CL CL];
            obj.CM=[obj.CM CM];
            obj.CN=[obj.CN CN];
        end
        
        function obj=f_solve_quick(obj)
            obj.CX=[];
            obj.CY=[];
            obj.CZ=[];
            
            obj.CL=[];
            obj.CM=[];
            obj.CN=[];
            n_wake=2*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
            for i=1:n_wake
                obj=obj.determine_boundary_conditions();
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                obj=obj.shed_wake();
            end
            obj=obj.f_postprocess;
        end
        
        function obj=f_solve(obj)
            obj.CX=[];
            obj.CY=[];
            obj.CZ=[];
            
            obj.CL=[];
            obj.CM=[];
            obj.CN=[];
            n_wake=obj.settings.wakelength_factor*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));          
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.Abb,obj.b');
            obj.Gamma_wake=obj.Cbb*obj.Gamma;
            obj=obj.shed_wake();
            while abs((sum(abs(obj.Gamma_prv))-sum(abs(obj.Gamma)))/sum(abs(obj.Gamma)))>10^(-10)
                obj=obj.determine_boundary_conditions();
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                obj=obj.shed_wake();
            end
            obj=obj.f_postprocess;
        end
        
        function obj=f_solve_for_Cl(obj,Cl_target)
            %% TODO: validate, enable entering of alpha start and tolerance
            alpha_0=0;
            tolerance=0.0001;
            
            
            obj.Uinf=a2bf(norm(obj.Uinf),alpha_0,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve();
            Cl_1=obj.CZ;
            alpha_0=5;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha_0,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve();
            Cl_prev=obj.CZ;
            dCLalpha=Cl_prev-Cl_1;
            dalpha=(Cl_target-Cl_prev)/dCLalpha;
            alpha=alpha_0+dalpha;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve();
            dCLalpha=(obj.CZ-Cl_prev)/dalpha;
            i=1;
            while abs(obj.CZ-Cl_target)>tolerance
                Cl_prev=obj.CZ;
                dalpha=(Cl_target-Cl_prev)/dCLalpha;
                alpha=alpha+dalpha;
                obj.Uinf=a2bf(norm(obj.Uinf),alpha,obj.state.beta,obj.Ma_corr);
                obj=obj.f_solve();
                i=i+1;
                dCLalpha=(obj.CZ-Cl_prev)/dalpha;
                if i>30
                    sprintf('f_solve_for_Cl did not converge!')
                    break;
                end
            end
            obj=obj.f_solve();
            obj.state=obj.state.set_alpha(alpha);
        end
        
        
        function obj=f_postprocess_exact(obj,varargin)
            
            U_x=(obj.Abb_x*obj.Gamma)';
            U_y=(obj.Abb_y*obj.Gamma)';
            U=norm(obj.Uinf')-(obj.Abw*obj.Gamma_wake)'*0;
            
            if obj.Ma_corr<1
                beta_inf=obj.Ma_corr;
            else
                beta_inf=1;
            end
            
            obj.cp_st=zeros(1,length(obj.Gamma));
            obj.cp_un=zeros(1,length(obj.Gamma));
            
            for i=1:length(obj.panels)
                
                sc=0.5*norm(obj.grid(1,obj.panels(4,i))-obj.grid(1,obj.panels(1,i)))+0.5*norm(obj.grid(1,obj.panels(3,i))-obj.grid(1,obj.panels(2,i)));
                ss=0.5*norm(obj.grid(2,obj.panels(2,i))-obj.grid(2,obj.panels(1,i)))+0.5*norm(obj.grid(2,obj.panels(4,i))-obj.grid(2,obj.panels(3,i)));
                
                if i>1
                    if obj.is_te(i-1)==1
                        bwd=i;
                        while (~obj.is_te(fwd))&&(~(fwd==length(obj.panels)))
                            fwd=fwd+1;
                        end
                        nbr=i+(fwd-bwd+1);
                        if nbr>length(obj.panels)
                            nbr=i;
                        end
                        obj.cp_st(i)=((obj.Gamma(i))*U(i)/sc)/(0.5*norm(obj.Uinf')^2)+((obj.Gamma(i)-obj.Gamma(nbr))*(U(i)+U_x(i)+U_y(i))*0.7071/ss)/(0.5*norm(obj.Uinf')^2);
                        if ~isempty(obj.Gamma_prv)
                            obj.cp_un(i)=(obj.Gamma(i)-obj.Gamma_prv(i))/(0.5*norm(obj.Uinf')^2*obj.t_step);
                        end
                    else
                        bwd=i;
                        fwd=i;
                        while (~(bwd==1))
                            if bwd>1
                                if(~obj.is_te(bwd-1))
                                    bwd=bwd-1;
                                else
                                    break;
                                end
                            else
                                break;
                            end
                        end
                        
                        while (~obj.is_te(fwd))&&(~(fwd==length(obj.panels)))
                            fwd=fwd+1;
                        end
                        nbr=i+(fwd-bwd+1);
                        if nbr>length(obj.panels)
                            nbr=i;
                        end
                        % distinguish left or right
                        obj.cp_st(i)=((obj.Gamma(i)-obj.Gamma(i-1))*U(i)/sc)/(0.5*norm(obj.Uinf')^2)+((obj.Gamma(i)-obj.Gamma(nbr))*(U(i)+U_x(i)+U_y(i))*0.7071/ss)/(0.5*norm(obj.Uinf')^2);
                        if ~isempty(obj.Gamma_prv)
                            %obj.cp_un(i)=((obj.Gamma(i)-obj.Gamma(i-1))-(obj.Gamma_prv(i)-obj.Gamma_prv(i-1)))/(0.5*norm(obj.Uinf')^2*obj.t_step);
                            obj.cp_un(i)=(obj.Gamma(i)-obj.Gamma_prv(i))/(0.5*norm(obj.Uinf')^2*obj.t_step);
                        end
                    end
                else
                    fwd=i;
                    while (~obj.is_te(fwd))&&(~(fwd==length(obj.panels)))
                        fwd=fwd+1;
                    end
                    nbr=i+(fwd);
                    if nbr>length(obj.panels)
                        nbr=i;
                    end
                    obj.cp_st(i)=((obj.Gamma(i))*U(i)/sc)/(0.5*norm(obj.Uinf')^2)+((obj.Gamma(i)-obj.Gamma(nbr))*(U(i)+U_x(i)+U_y(i))*0.7071/ss)/(0.5*norm(obj.Uinf')^2);
                    if ~isempty(obj.Gamma_prv)
                        obj.cp_un(i)=(obj.Gamma(i)-obj.Gamma_prv(i))/(0.5*norm(obj.Uinf')^2*obj.t_step);
                    end
                end
            end
            
            if ~isempty(obj.cp_un)
                obj.cp=1/beta_inf*(obj.cp_st)+obj.cp_un;
            else
                obj.cp=1/beta_inf*obj.cp_st;
            end
            
            alpha=atan(obj.Uinf(3)/obj.Uinf(1));
            beta=atan(obj.Uinf(2)/obj.Uinf(1));
            
            obj.F_body(1,:)=obj.qinf*obj.area.*obj.colloc_nvec(1,:).*obj.cp;
            obj.F_body(2,:)=obj.qinf*obj.area.*obj.colloc_nvec(2,:).*obj.cp;
            obj.F_body(3,:)=obj.qinf*obj.area.*obj.colloc_nvec(3,:).*obj.cp;
            
            %                 %alpha2=atan((obj.Uinf(3)+uvwind(3))/(obj.Uinf(1)+uvwind(1)));
            %                 %beta=atan((obj.Uinf(2)+uvwind(2))/(obj.Uinf(1)+uvwind(1)));
            %                 alpha=atan((obj.Uinf(3))/(obj.Uinf(1)));
            %                 beta=atan((obj.Uinf(2))/(obj.Uinf(1)));
            %                 %                alphai(i)=alpha-alpha2;
            
            M_BA=[  cos(alpha)*cos(beta) , -cos(alpha)*sin(beta), -sin(alpha);
                sin(beta)             , cos(beta)            ,    0;
                sin(alpha)*cos(beta) , -sin(alpha)*sin(beta), cos(alpha);];
            
            %             for i=1:length(obj.panels)
            %                 obj.F_aero(:,i)=M_BA'*obj.F_body(:,i);
            %             end
            
            CX=sum(obj.F_body(1,:))/(obj.qinf*obj.reference.S_ref);
            CY=sum(obj.F_body(2,:))/(obj.qinf*obj.reference.S_ref);
            CZ=sum(obj.F_body(3,:))/(obj.qinf*obj.reference.S_ref);
            
            %             Cl=sum(obj.F_aero(3,:))/(obj.qinf*obj.reference.S_ref);
            %             Cy=sum(obj.F_aero(2,:))/(obj.qinf*obj.reference.S_ref);
            
            Cdi=0;
            
            CL=sum(+obj.F_body(3,:).*obj.r(2,:)-obj.F_body(2,:).*obj.r(3,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
            CM=sum(-obj.F_body(3,:).*obj.r(1,:)+obj.F_body(1,:).*obj.r(3,:))/(obj.qinf*obj.reference.S_ref*obj.reference.c_ref);
            CN=sum(obj.F_body(2,:).*obj.r(1,:)-obj.F_body(1,:).*obj.r(2,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
            
            obj.CX=[obj.CX CX];
            obj.CY=[obj.CY CY];
            obj.CZ=[obj.CZ CZ];
            
            obj.CL=[obj.CL CL];
            obj.CM=[obj.CM CM];
            obj.CN=[obj.CN CN];
        end
        %          function obj=plot_nvecs(obj)
        % %             for i=1:length(obj.panels)
        % %                 handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.flag(i));
        % %                 hold on
        % %             end
        %
        %             quiver3(obj.colloc(1,:),obj.colloc(2,:),obj.colloc(3,:),obj.colloc_nvec(1,:),obj.colloc_nvec(2,:),obj.colloc_nvec(3,:),'r')
        %         end
        
        function grid=get_grid(obj)
            grid(1,:)=obj.grid(1,:)*obj.Ma_corr;
            grid(2,:)=obj.grid(2,:);
            grid(3,:)=obj.grid(3,:);
        end
        function obj=f_solve_std(obj)
            obj=obj.solve_MEX_vor3();
            %obj=obj.solve();
            obj=obj.f_postprocess();
        end
       
        function obj=f_solve_exp(obj)
            obj=obj.solve();
            obj=obj.f_postprocess();
        end
        
        
        function obj=plot_grid(obj)
            for i=1:length(obj.panels)
                if obj.is_te(i)==1
                    handle= fill3(obj.grid_vring(1,obj.panels(:,i)), obj.grid_vring(2,obj.panels(:,i)),obj.grid_vring(3,obj.panels(:,i)),'b');
                else
                    handle= fill3(obj.grid_vring(1,obj.panels(:,i)), obj.grid_vring(2,obj.panels(:,i)),obj.grid_vring(3,obj.panels(:,i)),'c');
                end
                hold on
            end
            
                        for i=1:length(obj.panels)
                            if obj.is_te(i)==1
                                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'y');
                            else
                                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'m');
                            end
                            hold on
                        end
            
            axis equal
            for i=1:length(obj.panels_wake)
                handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),'r');
                hold on
            end
            
            %             for i=1:length(obj.panels)
            %                 handle= fill3(obj.grid_vring(1,obj.panels(1:4,i)), obj.grid_vring(2,obj.panels(1:4,i)),obj.grid_vring(3,obj.panels(1:4,i)),'g');
            %                 hold on
            %             end
        end
        
        function obj=plot_gamma(obj)
            %   for i=1:length(obj.panels)
            %        %handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),0);
            %        handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.Gamma(i));
            %        hold on
            %   end
            axis equal
            for i=1:length(obj.panels_wake)
                %handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),uind);
                handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i))*obj.Ma_corr, obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),obj.Gamma_wake(i));
                hold on
            end
            
            for i=1:length(obj.panels)
                %handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),0);
                handle= fill3(obj.grid(1,obj.panels(:,i))*obj.Ma_corr, obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.Gamma(i));
                hold on
            end
            axis equal
            %             for i=1:length(obj.panels_wake)
            %                 %handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),uind);
            %                 handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),'w');
            %                 hold on
            %             end
        end
        
        function obj=plot_cp_wake(obj)
            for i=1:length(obj.panels)
                %handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),0);
                handle= fill3(obj.grid(1,obj.panels(:,i))*obj.Ma_corr, obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cp(i));
                hold on
                
            end
            axis equal
            for i=1:length(obj.panels_wake)
                %handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),uind);
                handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i))*obj.Ma_corr, obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),'w');
                hold on
            end
        end
        
        function obj=plot_uind(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),0);
                %  handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.Gamma(i));
                hold on
            end
            axis equal
            for i=1:length(obj.panels_wake)
                uind=(obj.uind(obj.panels_wake(1,i))+obj.uind(obj.panels_wake(2,i))+obj.uind(obj.panels_wake(3,i))+obj.uind(obj.panels_wake(4,i)))/4;
                handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),uind);
                %handle= fill3(obj.grid_wake(1,obj.panels_wake(1:4,i)), obj.grid_wake(2,obj.panels_wake(1:4,i)),obj.grid_wake(3,obj.panels_wake(1:4,i)),obj.Gamma_wake(i));
                hold on
            end
        end
        
        function obj=plot_L(obj)
            uwake=(obj.Abw*obj.Gamma_wake)';
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),uwake(i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_cz(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cz(i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_cy(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cy(i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_Fy(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.F_body(3,i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_cx(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cx(i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_vind(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),0);
                hold on
                quiver3(obj.fvap(1,i),obj.fvap(2,i),obj.fvap(3,i),obj.v_ind(1,i),obj.v_ind(2,i),obj.v_ind(3,i),5);
            end
            axis equal
        end
        
        function obj=plot_cp(obj)
            
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i))*obj.Ma_corr, obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cp(i)*obj.area(i));
                hold on
            end
            axis equal
            
            fileID = fopen('UVLM_vis.tp','w');
            
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z"\n');
            fprintf(fileID,'ZONE T="BIG ZONE", I=3, J=3, DATAPACKING=POINT');
            
            % determine nodal Cp
            nodal_cp=zeros(1,length(obj.grid));
            nodal_factor=zeros(1,length(obj.grid));
            
            for i=1:length(obj.panels)
                for j=1:4
                    nodal_cp(obj.panels(j,i))=nodal_cp(obj.panels(j,i))+obj.cp(i)*obj.area(i);
                    nodal_factor(obj.panels(j,i))=nodal_factor(obj.panels(j,i))+1;
                end
            end
            nodal_cp=nodal_cp./nodal_factor;
            
            fileID = fopen('VLM_vis.tp','w');
            
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            for i =1:length(obj.grid)
                fprintf(fileID,'%f %f %f %f\n',obj.grid(1,i),obj.grid(2,i),obj.grid(3,i),nodal_cp(i));
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end
            
            fclose(fileID);
        end
        
        function obj=plot_bc(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i))*obj.Ma_corr, obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.b(i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_wind(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cp(i));
                hold on
                quiver3(obj.colloc(1,i),obj.colloc(2,i),obj.colloc(3,i),obj.wind(i)*obj.colloc_nvec(1,i)/50,obj.wind(i)*obj.colloc_nvec(2,i)/50,obj.wind(i)*obj.colloc_nvec(3,i)/50)
            end
            axis equal
        end
        
        function obj=plot_flag(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.flag(i));
                hold on
            end
            axis equal
        end
        
        function obj=write_tecplot_wake(obj,filename,varargin)
            if nargin==3
                if varargin{1}==1
                    mode='w';
                    append=0;
                else
                    mode='a';
                    append=1;
                end
                
            end
            % determine nodal Cp
            nodal_cp=zeros(1,length(obj.grid));
            nodal_factor=zeros(1,length(obj.grid));
            
            for i=1:length(obj.panels)
                for j=1:4
                    nodal_cp(obj.panels(j,i))=nodal_cp(obj.panels(j,i))+obj.cp(i);
                    nodal_factor(obj.panels(j,i))=nodal_factor(obj.panels(j,i))+1;
                end
            end
            
            nodal_cp=nodal_cp./nodal_factor;
            
            if exist([filename '.tp'], 'file')==2
                delete([filename '.tp']);
            end
            fileID = fopen([filename '.tp'],mode);
            
            if append~=1
                fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
                fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            for i =1:length(obj.grid)
                fprintf(fileID,'%f %f %f %f\n',obj.grid(1,i)*obj.Ma_corr,obj.grid(2,i),obj.grid(3,i),nodal_cp(i));
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid_wake),length(obj.panels_wake));
            
            
            for i =1:length(obj.grid_wake)
                fprintf(fileID,'%f %f %f %f\n',obj.grid_wake(1,i)*obj.Ma_corr,obj.grid_wake(2,i),obj.grid_wake(3,i),0);
            end
            
            for i =1:length(obj.panels_wake)
                fprintf(fileID,'%i %i %i %i \n',obj.panels_wake(1,i),obj.panels_wake(2,i),obj.panels_wake(3,i),obj.panels_wake(4,i));
            end
            fclose(fileID);
        end
        
        function obj=write_tecplot_free_flying_wake(obj,filename,pos,varargin)
            mode='W';
            append=0;
            
            pos(1)=0;
            % determine nodal Cp
            nodal_cp=zeros(1,length(obj.grid));
            nodal_factor=zeros(1,length(obj.grid));
            
            for i=1:length(obj.panels)
                for j=1:4
                    nodal_cp(obj.panels(j,i))=nodal_cp(obj.panels(j,i))+obj.cp(i);
                    nodal_factor(obj.panels(j,i))=nodal_factor(obj.panels(j,i))+1;
                end
            end
            
            nodal_cp=nodal_cp./nodal_factor;
            
            if exist([filename '.tp'], 'file')==2
                delete([filename '.tp']);
            end
            fileID = fopen([filename '.tp'],mode);
            
            if append~=1
                fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
                fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            Euler=pos(4:6);
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
            
            for i =1:length(obj.grid)
                ri=obj.grid(:,i)-obj.reference.p_ref';

                pt=M_BI'*ri+pos(1:3);
                fprintf(fileID,'%f %f %f %f\n',pt(1)*obj.Ma_corr,pt(2),pt(3),nodal_cp(i));
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end

%             nodal_gamma_wake=zeros(1,length(obj.grid_wake));
%             nodal_factor_wake=zeros(1,length(obj.grid_wake));
%             
%             for i=1:length(obj.panels_wake)
%                 for j=1:4
%                     nodal_gamma_wake(obj.panels_wake(j,i))=nodal_gamma_wake(obj.panels_wake(j,i))+obj.Gamma_wake(i);
%                     nodal_factor_wake(obj.panels_wake(j,i))=nodal_factor_wake(obj.panels_wake(j,i))+1;
%                 end
%             end
%             
%             nodal_gamma_wake=nodal_gamma_wake./nodal_factor_wake;
%             
%             
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid_wake),length(obj.panels_wake));
            
            for i =1:length(obj.grid_wake)
                ri=obj.grid_wake(:,i)-obj.reference.p_ref';
                pt=M_BI'*ri+pos(1:3);
                fprintf(fileID,'%f %f %f %f\n',pt(1)*obj.Ma_corr,pt(2),pt(3),0);
            end
            
            for i =1:length(obj.panels_wake)
                fprintf(fileID,'%i %i %i %i \n',obj.panels_wake(1,i),obj.panels_wake(2,i),obj.panels_wake(3,i),obj.panels_wake(4,i));
            end
            fclose(fileID);
        end
        
        function obj=write_tecplot(obj,filename,varargin)
            if nargin==3
                if varargin{1}==1
                    mode='w';
                    append=0;
                else
                    mode='a';
                    append=1;
                end
                
            end
            % determine nodal Cp
            nodal_cp=zeros(1,length(obj.grid));
            nodal_factor=zeros(1,length(obj.grid));
            
            for i=1:length(obj.panels)
                for j=1:4
                    nodal_cp(obj.panels(j,i))=nodal_cp(obj.panels(j,i))+obj.cp(i);
                    nodal_factor(obj.panels(j,i))=nodal_factor(obj.panels(j,i))+1;
                end
            end
            
            nodal_cp=nodal_cp./nodal_factor;
            
            %             if exist([filename '.tp'], 'file')==2
            %                 delete([filename '.tp']);
            %             end
            
            fileID = fopen([filename '.tp'],mode);
            
            if append~=1
                fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
                fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            for i =1:length(obj.grid)
                fprintf(fileID,'%f %f %f %f\n',obj.grid(1,i)*obj.Ma_corr,obj.grid(2,i),obj.grid(3,i),nodal_cp(i));
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end
            
            fclose(fileID);
        end
        function obj=generate_state_space_matrices(obj,Vref,omega)
            %this function generates the state space matrices for a
            %continuous time state space model
            %input vector is Vn and Vndot, (normal velocities on
            %collocation point)
            %output is Gammab, Gammabdot,vsindx vsindy vsindz (induced
            %velocity on segments)
            
            obj=obj.f_compute_timestep(omega);
            
            %first the grid needs to be initialized
            obj=obj.initialize_unsteady_computation();
            
            
            %number of spanwise panels
            nS=sum(obj.is_te);
            %number of chordwise panels
            nC=length(obj.Abb)/nS;
            %number of bound panels
            nB=nC*nS;
            %xW describes the spacing of the wake grid
            xW=obj.grid_wake(1,obj.panels_wake(1,1+nS))-obj.grid_wake(1,obj.panels_wake(1,1));
            
            
            
            %the grid needs to be modified (small first wake panel)
            obj.grid_wake(1,:)=obj.grid_wake(1,:)-xW*0.25;
            obj.grid_wake(1,nS+2+1:end)=obj.grid_wake(1,nS+2+1:end)-xW*0.75;
            obj.grid_vring(1,nC+2:nC+2:end)=obj.grid_vring(1,nC+2:nC+2:end)-0.25*xW;
            obj=obj.compute_influence_coeff_matrix();

            % number of ALL wake panels
            nW=size(obj.Abw,2);
            % number of REST wake panels
            nWR=nW-nS;
            
            %K1 is the influence coefficients of bound on bound panels
            K1=-obj.Abb;
            %K2 is the influence coefficients of small first wake row on bound
            K2=-obj.Abw(:,1:nS);
            %K3 is the influence coefficients of wake on bound
            K3=-obj.Abw(:,nS+1:end);
            % same for M1X M1Y M1Z M2X M2Y M2Z ...
            M1X=-obj.Abb_x;
            M2X=-obj.Abw_x(:,1:nS);
            M3X=-obj.Abw_x(:,nS+1:end);
            M1Y=-obj.Abb_y;
            M2Y=-obj.Abw_y(:,1:nS);
            M3Y=-obj.Abw_y(:,nS+1:end);
            M1Z=-obj.Abb_z;
            M2Z=-obj.Abw_z(:,1:nS);
            M3Z=-obj.Abw_z(:,nS+1:end);
            %K4 maps all body gamma to vector of gammas of last row
            K4=diag(obj.is_te)'*1;
            K4(all(K4==0,2),:)=[];
            %K5 is an eye matrix of the number of spanwise panels 
            K5=-eye(nS); 
            %K6 is transport within rest wake
            K6=([zeros(nS,nWR) ;eye(nWR-nS) zeros(nWR-nS,nS)]-eye(nWR))*Vref/xW;
            %K7 is transport within first wake row
            K7=[eye(nS); zeros(nWR-nS,nS)]*Vref/xW;
            %%
            K8=+K6+K7*(K5-K4*K1^-1*K2)^-1*K4*K1^-1*K3;
            K9=K7*(K5-K4*K1^-1*K2)^-1*K4*K1^-1;

            L1=[zeros(1,nB); -diag(~obj.is_te(1:end-1)') zeros(nB-1,1)]+ eye(nB,nB);
            L2=1/2*([zeros(1,nB); diag(~obj.is_te(1:end-1)') zeros(nB-1,1)]+ eye(nB,nB));

            L3=(K2*K5^-1*K4-K1)^-1*K3;
            L4=(K2*K5^-1*K4-K1)^-1;
            L5=L3*K8;
            L6=L3*K9;

            L7=(M1X-M2X*K5^-1*K4)*L3+M3X;
            L8=(M1X-M2X*K5^-1*K4)*L4;
            L9=(M1Y-M2Y*K5^-1*K4)*L3+M3Y;
            L10=(M1Y-M2Y*K5^-1*K4)*L4;
            L11=(M1Z-M2Z*K5^-1*K4)*L3+M3Z;
            L12=(M1Z-M2Z*K5^-1*K4)*L4;
            
            Ass=K8;
            Bss=[K9 zeros(nWR,nB)];
            Css=[L1*L3; L2*L5; L7; L9; L11];
            Dss=[L1*L4 zeros(nB,nB); L2*L6 L2*L4; L8 zeros(nB,nB); L10 zeros(nB,nB); L12 zeros(nB,nB)];
            obj.sys=ss(Ass,Bss,Css,Dss);
            %restore original grid
            obj=obj.initialize_unsteady_computation();
        end 
        function obj=generate_reduced_parameter_varying_ssm(obj,nOrder,omega)
            % this function creates the ssm and the matrices to make it an parameter varying ssm which is valid for all speeds
            %therefore a and b need to be multiplied with the actual
            %velocity
            
            %c and d need to be multiplied by the sum of the product of the
            %actual velocity and cf +cf2
            %c_res=V*(cf+cf2)*c
            %d_res=V*(df+df2)*d
            %for nOrder==0 no reduction will be carried out
            %create state space for unit velocity
            obj=obj.generate_state_space_matrices(1,omega);
            %reduce state space model
            if nOrder~=0
                obj.sys=balred(obj.sys,nOrder);
            end
            
            obj.lpvSys.a=obj.sys.a;
            obj.lpvSys.b=obj.sys.b;
            obj.lpvSys.c=obj.sys.c;
            obj.lpvSys.d=obj.sys.d;
            nB=size(obj.lpvSys.c,1)/5;
            nStates=size(obj.lpvSys.a,1);
            obj.lpvSys.cF=[zeros(nB,nStates); ones(nB,nStates); zeros(3*nB,nStates)];
            obj.lpvSys.cF2=(obj.lpvSys.cF==0);
            obj.lpvSys.dF=[zeros(nB) zeros(nB); ones(nB) zeros(nB); zeros(3*nB,2*nB)];
            obj.lpvSys.dF2=(obj.lpvSys.dF==0);
        end
    end
    
end

