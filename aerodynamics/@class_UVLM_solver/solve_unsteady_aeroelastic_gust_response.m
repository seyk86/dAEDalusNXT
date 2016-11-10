%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj=solve_unsteady_aeroelastic_gust_response(obj,aircraft,aircraft_structure,Uds,H,t_start,t_end,init_def)
            n_beam=length(aircraft_structure.beam);
            
            obj.t_step=obj.reference.c_ref/(norm(obj.Uinf)*32);
            
            n_wake=obj.settings.wakelength_factor*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
            
            obj=obj.initialize_vring_grid(obj.t_step);
            
            obj=obj.initialize_wake_grid(ceil(n_wake),obj.t_step);
            
            steps=60;
            
            %t_start=0.5;
            
            %filename = ['unsteady_time_domain_solution1.gif'];
            %% (1-cos)gust
            s=0:(2*H/(20)):2*H;
            gust_shape=Uds/2*(1-cos(pi*s/H));
            x_min=min(obj.grid(1,:));
            x_max=max(obj.grid_wake(1,:));
            
            
            delta_x=obj.Uinf(1)*obj.t_step;
            obj.x_kin=x_min:delta_x:x_max;
            obj.Uinf_x=zeros(3,length(obj.x_kin));
            obj.Uinf_x(1,:)=obj.Uinf(1);
            obj.Uinf_x(2,:)=obj.Uinf(2);
            obj.Uinf_x(3,:)=obj.Uinf(3);
            gust_shape_interp=interp1(s,gust_shape,0:delta_x:2*H);
            d_gust_shape_interp=[ diff(gust_shape_interp) 0];
            
            obj.Gamma_wake=zeros(length(obj.panels_wake),1);
            t_vec=0:obj.t_step:t_end;
            obj.t_vec=t_vec;
            obj=obj.compute_influence_coeff_matrix();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.Abb,obj.b');
            obj.Gamma_wake=obj.Cbb*obj.Gamma;
            obj=obj.f_postprocess;
            
            % initialize newmark beta timemarching scheme
            
            beta=1/4;
            gamma=1/2;
            dt=obj.t_step;
            %% check Ma
            x_now=init_def;%zeros(length(aircraft_structure.Ftest),1);
            xdot_now=zeros(length(aircraft_structure.Ftest),1);
            xdotdot_now=zeros(length(aircraft_structure.Ftest),1);
            x_nxt=zeros(length(aircraft_structure.Ftest),1);
            xdot_nxt=zeros(length(aircraft_structure.Ftest),1);
            xdotdot_nxt=zeros(length(aircraft_structure.Ftest),1);
            
            xdotdot_now=linsolve(aircraft_structure.Mff,aircraft_structure.Ftest-aircraft_structure.Kff*x_now);
            
            n=1;
            k=1;
            l=1;
            m=1;
            
            mkdir(['results/' obj.case_name],'gust');
            
            obj.u_test=zeros(6,length(t_vec));
            acc_z=0;
            C=0.0035*aircraft_structure.Kff+0.0035*aircraft_structure.Mff;
            max_loads.loads=aircraft_structure.beam(1).node_loadings_loc*0;
%             for bb=2:n_beam
%                 max_loads(bb).loads=aircraft_structure.beam(bb).node_loadings_loc*0;
%             end
            gust=0;
            for i=1:length(t_vec)
                if i==1
                    obj=obj.initialize_wake_circulation_pqr(0,0,0);
                end
                
                acc_z=0;
                if t_vec(i)>=t_start
                    gust=1;
                    obj.Uinf_x(1,:)=obj.Uinf(1);
                    obj.Uinf_x(2,:)=obj.Uinf(2);
                    obj.Uinf_x(3,:)=obj.Uinf(3);
                    if length(obj.Uinf_x)<length(gust_shape_interp)
                        if n<=length(obj.Uinf_x)
                            obj.Uinf_x(3,1:n)=obj.Uinf(3)+gust_shape_interp(end-n+1:end);
                          %  acc_z=d_gust_shape_interp(end-n+1:end);
                            n=n+1;
                        elseif k+n<=length(gust_shape_interp)
                            obj.Uinf_x(3,1:end)=obj.Uinf(3)+gust_shape_interp(end-k-length(obj.x_kin):end-k-1);
                           % acc_z=d_gust_shape_interp(end-k-length(obj.x_kin):end-k-1);
                            k=k+1;
                        else
                            obj.Uinf_x(3,l+1:end)=obj.Uinf(3)+gust_shape_interp(1:length(obj.x_kin)-l);
                          %  acc_z=d_gust_shape_interp(1:length(obj.x_kin)-l);
                            l=l+1;
                        end
                    else
                        if n<length(gust_shape_interp)
                            obj.Uinf_x(3,1:n)=obj.Uinf_x(3,k:k+n-1)+gust_shape_interp(n:-1:1);
                            %acc_z=d_gust_shape_interp(n:-1:1);
                            n=n+1;
                        elseif k+n-1>length(obj.x_kin)
                            obj.Uinf_x(3,k:end)=obj.Uinf_x(3,k:end)+gust_shape_interp(1:length(obj.x_kin)-k+1);
                            %acc_z=d_gust_shape_interp(1:length(obj.x_kin)-k+1);
                            k=k+1;
                        else
                            obj.Uinf_x(3,k:k+n-1)=obj.Uinf_x(3,k:k+n-1)+gust_shape_interp(1:n);
                           % acc_z=d_gust_shape_interp(1:n);
                            k=k+1;
                        end
                    end
                elseif k>=length(obj.x_kin)
                    gust=0;
                     acc_z=0;
                    obj.Uinf_x(1,:)=obj.Uinf(1);
                    obj.Uinf_x(2,:)=obj.Uinf(2);
                    obj.Uinf_x(3,:)=obj.Uinf(3);
                end
                
%                 plot(obj.x_kin,obj.Uinf_x(3,:),'x')
%                 drawnow
                
                
                if isempty(acc_z)
                    acc_z=0;
                end
                
                prv_colloc=obj.colloc;
                
                aircraft=aircraft.compute_beam_forces(obj.F_body,aircraft_structure);
                
                if n_beam>1
                    for bb=1:n_beam-1
                        aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
                    end
                else
                    for bb=1:n_beam
                        aircraft_structure.beam(bb)=aircraft_structure.beam(bb).f_set_aeroloads(aircraft.wings(bb));
                    end
                end
                
                if i>1
                    aircraft_structure=aircraft_structure.f_set_acceleration([0,0,acc_z-9.81,0,0,0],[0 0 0]);
                    aircraft_structure=aircraft_structure.f_assemble(1,0);
                    
                    x_nxt=linsolve(aircraft_structure.Mff+beta*dt^2*aircraft_structure.Kff+gamma*dt*C,beta*dt^2*aircraft_structure.Ftest+aircraft_structure.Mff*x_now+dt*aircraft_structure.Mff*xdot_now+dt^2*aircraft_structure.Mff*(1/2-beta)*xdotdot_now...
                        +C*beta*dt^2*(gamma/(beta*dt)*x_now+(gamma/beta-1)*xdot_now+1/2*dt*(gamma/beta-2)*xdotdot_now));
                    
                    %aircraft_structure=aircraft_structure.f_assemble(1,0);
                    gamma=1/2;
                    if i<10
                        xdotdot_nxt=1/(beta*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta)*xdotdot_now)*0;
                    else
                        xdotdot_nxt=1/(beta*dt^2)*(x_nxt-x_now-dt*xdot_now-dt^2*(1/2-beta)*xdotdot_now);
                    end
                    xdot_nxt=xdot_now+dt*((1-gamma)*xdotdot_now+gamma*xdotdot_nxt);
                else
                    x_nxt=x_now;
                    xdot_nxt=xdot_now;
                    xdotdot_nxt=xdotdot_now;
                end

                aircraft_structure.nodal_deflections=x_nxt;
                
                % for debugging only
                if n_beam>1
                    obj.u_test(:,i)=x_nxt(19:24);
                else
                    obj.u_test(:,i)=x_nxt(1:6);
                end
                %
                
                aircraft_structure=aircraft_structure.f_postprocess();
                
%                 for bb=1:n_beam
%                     max_loads(bb).loads=max(max_loads(bb).loads,aircraft_structure.beam(bb).node_loadings_loc);
%                 end
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                
                obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
                obj.grid(2,:)=aircraft.grid_deflected(2,:);
                obj.grid(3,:)=aircraft.grid_deflected(3,:);
                
                obj=obj.update_grid();
                
                obj=obj.determine_boundary_conditions_gust_deflection(prv_colloc,obj.t_step);
                
                obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
                bx=obj.b-(obj.Abw*obj.Gamma_wake)';
                obj.Gamma_prv=obj.Gamma;
                obj.Gamma=linsolve(obj.Abb,bx');
                %obj.Gamma_wake=obj.Cbb*obj.Gamma+obj.Cbw*obj.Gamma_wake;
                obj=obj.f_postprocess;
                if gust==1
                    obj=obj.shed_wake_gust();
                else
                    obj=obj.shed_wake();
                end
                % obj=obj.perform_wake_rollup;
                if obj.settings.movie==1
                    obj.write_tecplot_wake(['results/' obj.case_name '/gust/gust_response_wake' num2str(m)],1);
                end
                m=m+1;
                x_now=x_nxt;
                xdot_now=xdot_nxt;
                xdotdot_now=xdotdot_nxt;
                m*obj.t_step
                drawnow
            end
            obj.max_gust_loads=max_loads;
            
        end
