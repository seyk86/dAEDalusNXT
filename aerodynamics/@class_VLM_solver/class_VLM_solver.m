%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_VLM_solver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> aerodynamic state used for solution (type class_aerodynamic_state)
        state;
        %> grid point coordinates
        grid;
        %> panel edge indices of grid
        panels;
        %> trailing edge indices
        te_idx;
        %> collocation point coordinates
        colloc;
        %> normal vector on collocation point
        colloc_nvec;
        %> aerodynamic influence coefficient matrix
        C;
        
        Cind;
        wind;
        b;
        %> vorticity vector
        Gamma;
        F_body;
        M_body;
        F_aero;
        Uinf;
        cl;
        cdi;
        cp;
        
        cz;
        cy;
        cx;
        cl_n;
        rho;
        qinf;
        A;
        r;% lever from 1/4 midpoint to p_ref
        r_dwn; % lever from 3/4 midpoint to p_ref
        dvortex;
        
        Ma_corr;
        
        %%%% RESULTS
        %% aerodynamic coefficients
        %> Coefficients in body axis
        CX;
        CY;
        CZ;
        Cl;
        Cy;
        Cdi;
        CL;
        CM;
        CN;
        %> u Derivatives in body axis
        CXu;
        CYu;
        CZu;
        Clu;
        Cyu;
        Cdiu;
        CLu;
        CMu;
        CNu;
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
        
        %% Reference Values
        
        reference;
        
        flag;
        
        crossp;
        
    end
    
    methods
        function obj=class_VLM_solver(grid,te_idx,panels,state,reference)
            obj.state=state;
            obj.rho=state.rho_air;
            obj.Ma_corr=sqrt(1-state.Ma^2);
            obj.state.p_ref=reference.p_ref;
            obj.grid=grid;
            obj.panels=panels;

            obj.te_idx=te_idx;
            obj.reference=reference;
            obj.r=zeros(3,length(obj.panels));
            for i=1:length(obj.panels)
                    obj.r(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.75+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.25-obj.state.p_ref';
                obj.r_dwn(:,i)=0.5*(obj.grid(:,obj.panels(2,i))+obj.grid(:,obj.panels(1,i)))*0.25+0.5*(obj.grid(:,obj.panels(3,i))+obj.grid(:,obj.panels(4,i)))*0.75-obj.state.p_ref';
            end
            if obj.Ma_corr<1
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,obj.Ma_corr);
                obj.grid(1,:)=grid(1,:)/obj.Ma_corr;
                obj.grid(2,:)=grid(2,:);
                obj.grid(3,:)=grid(3,:);
            else
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,1);
            end
            obj.qinf=1/2*state.rho_air*norm(obj.Uinf)^2;
            
            obj=obj.compute_colloc_points();
        end
        
        function obj=f_set_state(obj,state)
            pgprev=obj.Ma_corr;
            obj.Ma_corr=sqrt(1-state.Ma^2);
            if obj.Ma_corr<1
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,obj.Ma_corr);
                obj.grid(1,:)=obj.grid(1,:)*pgprev/obj.Ma_corr;
%                 obj.grid(2,:)=obj.grid(2,:);
%                 obj.grid(3,:)=grid(3,:);
            else
                obj.Uinf=a2bf(state.V_A,state.alpha,state.beta,1);
            end
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
        
        function new = f_copy(this)
            % Instantiate new object of the same class.
            class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,wingaero.state,aircraft.reference);
            new = eval([class(this) '(' this.grid ',' this.te_idx ',' this.panels, ',' this.state ',' this.referencen ')']);
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
        
        function obj=compute_colloc_points(obj)
            for i=1:1:length(obj.panels)
                obj.colloc(:,i)=0.125*obj.grid(:,obj.panels(1,i))+0.125*obj.grid(:,obj.panels(2,i))+0.375*obj.grid(:,obj.panels(3,i))+0.375*obj.grid(:,obj.panels(4,i));
                n_vec=cross(obj.grid(:,obj.panels(3,i))-obj.grid(:,obj.panels(1,i)),obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(2,i)));
                obj.colloc_nvec(:,i)=n_vec/norm(n_vec);
            end
        end
        
        function vortex=compute_vortex_coords(obj,panel_idx)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] V*1E10];
            if obj.flag(panel_idx)==0
                pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
                po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
                vortex=[pi pi po po]+obj.dvortex;
            else
                p_i=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
                p_o=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
                vortex=[obj.grid(:,obj.panels(4,panel_idx)) obj.grid(:,obj.panels(4,panel_idx)) p_i p_o obj.grid(:,obj.panels(3,panel_idx)) obj.grid(:,obj.panels(3,panel_idx))];
                
                p_i=obj.grid(:,obj.panels(4,panel_idx));
                p_o=obj.grid(:,obj.panels(3,panel_idx));
                test=1;
                while 1
                    i_min_i=1;
                    j_min_i=1;
                    min_i=1e7;
                    for i=1:length(obj.panels)
                        for j=1:4
                            vec=obj.grid(:,obj.panels(j,i))-p_i;
                            
                            dp=sum(conj(vec).*obj.Uinf')/(norm(vec)*norm(obj.Uinf));
                            %dp=dot(vec,obj.Uinf)/(norm(vec)*norm(obj.Uinf));
                            if acosd(dp)<80 && acosd(dp)>-80
                                if ((norm(vec)+(1-dp)*norm(vec))<min_i) && (norm(vec)>0)
                                    min_i=norm(vec)+(1-dp)*norm(vec);
                                    i_min_i=i;
                                    j_min_i=j;
                                end
                            end
                        end
                    end
                    if test==2
                        vortex=[obj.grid(:,obj.panels(j_min_i,i_min_i)) vortex];
                        break;
                    end
                    if obj.flag(i_min_i)==0
                        test=2;
                    else
                        vortex=[obj.grid(:,obj.panels(j_min_i,i_min_i)) vortex];
                        p_i=obj.grid(:,obj.panels(j_min_i,i_min_i));
                    end
                end
                test=1;
                while 1
                    i_min_o=1;
                    j_min_o=1;
                    min_o=1e7;
                    for i=1:length(obj.panels)
                        for j=1:4
                            vec=obj.grid(:,obj.panels(j,i))-p_o;
                            dp=sum(conj(vec).*obj.Uinf')/(norm(vec)*norm(obj.Uinf));
                            %dp=dot(vec,obj.Uinf)/(norm(vec)*norm(obj.Uinf));
                            if acosd(dp)<80 && acosd(dp)>-80
                                if ((norm(vec)+(1-dp)*norm(vec))<min_o) && (norm(vec)>0)
                                    min_o=norm(vec)+(1-dp)*norm(vec);
                                    i_min_o=i;
                                    j_min_o=j;
                                end
                            end
                        end
                    end
                    if test==2
                        vortex=[vortex obj.grid(:,obj.panels(j_min_o,i_min_o)) ];
                        break;
                    end
                    if obj.flag(i_min_o)==0
                        test=2;
                    else
                        vortex=[vortex obj.grid(:,obj.panels(j_min_o,i_min_o)) ];
                        p_o=obj.grid(:,obj.panels(j_min_o,i_min_o));
                    end
                end
                
                vortex=[vortex(:,1)+V*1E10 vortex vortex(:,end)+V*1E10];
            end
            plot3(vortex(1,:),vortex(2,:),vortex(3,:),'LineWidth',1)
            hold on
        end
        
        %         function vfield=compute_velocity_field(obj)
        %             n_x=100;
        %             n_y=100;
        %             n_z=100;
        %             min_x=1.5*min(obj.grid(1,:));
        %             max_x=1.5*max(obj.grid(1,:));
        %             min_y=1.5*min(obj.grid(2,:));
        %             max_y=1.5*max(obj.grid(2,:));
        %             min_z=1.5*min(obj.grid(3,:));
        %             max_z=1.5*max(obj.grid(3,:));
        %             dx=(x_max-x_min)/n_x;
        %             dy=(y_max-y_min)/n_y;
        %             dz=(z_max-z_min)/n_z;
        %             field=meshgrid(x_min:dx:x_max,y_min:dy:y_max,z_min:dz:z_max);
        %
        %             for xx=1:n_x
        %                 for yy=1:n_y
        %                     for zz:1:n_z
        %                         w=[0 0 0];
        %                         for colloc_idx=1:length(obj.panels)
        %                             vortex=[obj.grid(:,obj.te_idx(obj.panels(4,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(4,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(3,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(3,colloc_idx)))]+obj.dvortex;
        %                             for i=1:length(vortex)-1
        %                                 r1=field(xx,yy,zz)-vortex(:,i);
        %                                 r2=coords(:,j)-vortex(:,i+1);
        %                                 r0=vortex(:,i+1)-vortex(:,i);
        %                                 r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
        %                                 r1dr2=(r1/norm(r1)-r2/norm(r2));
        %                                 wAB=-obj.Gamma(colloc_idx)*r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
        %                                 coredist=norm((norm(cross(coords(:,j)-vortex(:,i),r0)))/(norm(r0)));
        %
        %                                 if (~isnan(wAB(1))) && (coredist>1e-8)
        %                                     w=w+wAB;
        %                                 end
        %                             end
        %                         end
        %                         wstream(:,j)=w;
        %                     end
        %                 end
        %             end
        %
        %
        %         end
        
        function streamlines=compute_streamlines(obj,coords,timestep,nsteps)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] V*1E10];
            wstream=zeros(3,length(coords));
            streamlines=zeros(nsteps,3,length(coords));
            streamlines(1,:,:)=coords;
            for k=2:nsteps
                for j=1:length(coords)
                    w=[0 0 0];
                    for colloc_idx=1:length(obj.panels)
                        vortex=[obj.grid(:,obj.te_idx(obj.panels(4,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(4,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(3,colloc_idx))) obj.grid(:,obj.te_idx(obj.panels(3,colloc_idx)))]+obj.dvortex;
                        for i=1:length(vortex)-1
                            r1=coords(:,j)-vortex(:,i);
                            r2=coords(:,j)-vortex(:,i+1);
                            r0=vortex(:,i+1)-vortex(:,i);
                            r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                            r1dr2=(r1/norm(r1)-r2/norm(r2));
                            wAB=-obj.Gamma(colloc_idx)*r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                            coredist=norm((norm(cross(coords(:,j)-vortex(:,i),r0)))/(norm(r0)));
                            
                            if (~isnan(wAB(1))) && (coredist>1e-8)
                                w=w+wAB;
                            end
                        end
                    end
                    wstream(:,j)=w;
                end
                
                for i=1:length(coords)
                    coords_nxt(:,i)=coords(:,i)+(V+wstream(:,i))*timestep;
                end
                coords=coords_nxt;
                streamlines(k,:,:)=coords_nxt;
            end
            
        end
        
        function Cdi=compute_trefftz_drag_mex(obj)
            %                         dist=-200;
            %                         u=obj.Uinf(1);
            %                         v=obj.Uinf(2);
            %                         w=obj.Uinf(3);
            %                         k=1;
            %                         D=0;
            %
            %                         for i=1:length(obj.panels)
            %                                 te_x=obj.grid(1,obj.te_idx(obj.panels(4,i)));
            %                                 te_y=obj.grid(2,obj.te_idx(obj.panels(4,i)));
            %                                 te_z=obj.grid(3,obj.te_idx(obj.panels(4,i)));
            %                                 lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
            %
            %                                 treffz_coords_in(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
            %
            %                                 te_x=obj.grid(1,obj.te_idx(obj.panels(3,i)));
            %                                 te_y=obj.grid(2,obj.te_idx(obj.panels(3,i)));
            %                                 te_z=obj.grid(3,obj.te_idx(obj.panels(3,i)));
            %                                 lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
            %                                 treffz_coords_out(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
            %
            %
            %                                 k=k+1;
            %                         end
            
            [D,x,y,u,v]=compute_trefftz_drag(obj.grid,obj.panels,obj.te_idx,obj.Gamma,obj.Uinf,230);
            %
            %                                     figure
            %                          contour(x,y,sqrt(u.^2+v.^2),100)
            %                          hold on
            %                         quiver(x,y,u,v);
            %                         hold on
            %                         plot(treffz_coords_in(2,:),treffz_coords_in(3,:),'ro');
            %                         plot(treffz_coords_out(2,:),treffz_coords_out(3,:),'cx');
            Cdi=D/(norm(obj.Uinf)^2*obj.reference.S_ref);
        end
        
                
        function Cdi=plot_trefftz_drag_mex(obj)
                                    dist=-200;
                                    u=obj.Uinf(1);
                                    v=obj.Uinf(2);
                                    w=obj.Uinf(3);
                                    k=1;
                                    D=0;
            
                                    for i=1:length(obj.panels)
                                            te_x=obj.grid(1,obj.te_idx(obj.panels(4,i)));
                                            te_y=obj.grid(2,obj.te_idx(obj.panels(4,i)));
                                            te_z=obj.grid(3,obj.te_idx(obj.panels(4,i)));
                                            lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
            
                                            treffz_coords_in(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
            
                                            te_x=obj.grid(1,obj.te_idx(obj.panels(3,i)));
                                            te_y=obj.grid(2,obj.te_idx(obj.panels(3,i)));
                                            te_z=obj.grid(3,obj.te_idx(obj.panels(3,i)));
                                            lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
                                            treffz_coords_out(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
            
            
                                            k=k+1;
                                    end
            
            [D,x,y,u,v]=compute_trefftz_drag(obj.grid,obj.panels,obj.te_idx,obj.Gamma,obj.Uinf,500);
            
                                                figure
                                     contourf(x,y,sqrt(u.^2+v.^2),200)
                                     hold on
                                    %quiver(x,y,u,v);
                                    hold on
                                    plot(treffz_coords_in(2,:),treffz_coords_in(3,:),'ro');
                                    plot(treffz_coords_out(2,:),treffz_coords_out(3,:),'cx');
            Cdi=D/(norm(obj.Uinf)^2*obj.reference.S_ref);
        end
        
        function   Cdi=compute_treffz_drag(obj)
            
            dist=-0.12;
            u=obj.Uinf(1);
            v=obj.Uinf(2);
            w=obj.Uinf(3);
            k=1;
            D=0;
            
            
            
            for i=1:length(obj.panels)
                te_x=obj.grid(1,obj.te_idx(obj.panels(4,i)));
                te_y=obj.grid(2,obj.te_idx(obj.panels(4,i)));
                te_z=obj.grid(3,obj.te_idx(obj.panels(4,i)));
                lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
                
                treffz_coords_in(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
                
                te_x=obj.grid(1,obj.te_idx(obj.panels(3,i)));
                te_y=obj.grid(2,obj.te_idx(obj.panels(3,i)));
                te_z=obj.grid(3,obj.te_idx(obj.panels(3,i)));
                lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
                treffz_coords_out(:,k)=[te_x;te_y;te_z]+lambda*[u;v;w];
                
                
                k=k+1;
            end
            
            bmax=3*(max(treffz_coords_in(2,:))-min(treffz_coords_in(2,:)));
            start_y=min(treffz_coords_in(2,:))-0.333*bmax;
            start_z=min(treffz_coords_in(3,:))-0.5*bmax;
            
            
            n_pan=400;
            span_multiplier=1;
            
            dl=bmax/n_pan;
            V=0;
            for i=1:n_pan
                for j=1:n_pan
                    
                    V=0;
                    %ingoing
                    ry=(start_y+dl/2)+i*dl-treffz_coords_in(2,:);
                    rz=(start_z+dl/2)+j*dl-treffz_coords_in(3,:);
                    %                     hold on
                    %                     for kk=1:length(ry)
                    %
                    %                     line([treffz_coords_in(2,kk) treffz_coords_in(2,kk)+ry(kk)],[treffz_coords_in(3,kk) treffz_coords_in(3,kk)+rz(kk)]);
                    %                     end
                    r=[rz; -ry]';
                    lr=sqrt(sum(r.^2,2));
                    er=r(:,1)./lr;
                    er(:,2)=r(:,2)./lr;
                    
                    v1=obj.Gamma(:,1)./(2*pi*lr).*er(:,1);
                    v1(:,2)=obj.Gamma(:,1)./(2*pi*lr).*er(:,2);
                    
                    %                     hold on
                    %                     for kk=1:length(v1)
                    %
                    %                     quiver(treffz_coords_in(2,kk),treffz_coords_in(3,kk),v1(kk,1),v1(kk,2));
                    %                     end
                    
                    
                    
                    %outgoing
                    ry=(start_y+dl/2)+i*dl-treffz_coords_out(2,:);
                    rz=(start_z+dl/2)+j*dl-treffz_coords_out(3,:);
                    
                    r=[rz; -ry]';
                    lr=sqrt(sum(r.^2,2));
                    er=r(:,1)./lr;
                    er(:,2)=r(:,2)./lr;
                    
                    v2=-obj.Gamma(:,1)./(2*pi*lr).*er(:,1);
                    v2(:,2)=-obj.Gamma(:,1)./(2*pi*lr).*er(:,2);
                    
                    V=sum([v1;v2]);
                    
                    
                    u(i,j)=V(1);
                    v(i,j)=V(2);
                    x(i,j)=(start_y+dl/2)+i*dl;
                    y(i,j)=(start_z+dl/2)+j*dl;
                    D(i,j)=0.5*obj.rho*(sum(V.^2))*dl^2;
                    
                    
                    % VVEC(i+1,j+1,:)=V;
                    % VTOT(i+1,j+1)=sqrt(sum(V.^2));
                    
                end
            end
            figure
            quiver(x,y,u,v);
            hold on
            plot(treffz_coords_in(2,:),treffz_coords_in(3,:),'ro');
            plot(treffz_coords_out(2,:),treffz_coords_out(3,:),'cx');
            
            %           results.treffts.x=x;
            %           results.treffts.y=y;
            %           results.treffts.D=D;
            %           results.treffts.VTOT=VTOT;
            %           results.treffts.VVEC=VVEC;
            Cdi=sum(sum(D))/(0.5*norm(obj.Uinf)^2*obj.S_ref*obj.rho)
            figure
                        hold on
                        %plot3(treffz_coords_in(1,:),treffz_coords_in(2,:),treffz_coords_in(3,:),'x');
                        plot3(dist*ones(n_pan,n_pan),x,y,'rx')
                        plot3(treffz_coords_in(1,:),treffz_coords_in(2,:),treffz_coords_in(3,:),'x');
                        plot3(treffz_coords_out(1,:),treffz_coords_out(2,:),treffz_coords_out(3,:),'ro');
            quiver3(treffz_coords(1,:),treffz_coords(2,:),treffz_coords(3,:),wtreffz(1,:),wtreffz(2,:),wtreffz(3,:));
        end
        
        
        function wij=compute_downwash(obj,panel_idx,colloc_idx)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] V*1E10];
            
            pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
            po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
            w=0;
            vortex=[pi pi po po]+obj.dvortex;
            for i=1:length(vortex)-1
                r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                r0=vortex(:,i+1)-vortex(:,i);
                %r1xr2=cross(r1,r2);
                r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                r1dr2=(r1/norm(r1)-r2/norm(r2));
                wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                if ~isnan(wAB(1))
                    w=w+wAB;
                end
            end
            %wij=dot(w,obj.colloc_nvec(:,panel_idx));
            wij=sum(conj(w).*obj.colloc_nvec(:,panel_idx)');
        end
        
        function [wij,wijind]=compute_downwash_ad(obj,panel_idx,colloc_len)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] [0;0;0] [0;0;0] V*1E10];
            obj.dvortex=[V*1E10  [0;0;0] [0;0;0] V*1E10];
            pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
            po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
            vortex=[pi pi po po]+obj.dvortex;
            %vortex=[obj.grid(:,obj.panels(4,panel_idx)) obj.grid(:,obj.panels(4,panel_idx)) pi po obj.grid(:,obj.panels(3,panel_idx)) obj.grid(:,obj.panels(3,panel_idx))]+obj.dvortex;
            wijind=zeros(1,colloc_len);
            wij=zeros(1,colloc_len);
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                wind=[0 0 0];
                for i=1:length(vortex)-1
                    r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                    r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                    r0=vortex(:,i+1)-vortex(:,i);
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB;
                        if i~=2
                            wind=wind+wAB;
                        end
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wijind(1,colloc_idx)=sum(conj(wind).*obj.colloc_nvec(:,panel_idx)');
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,panel_idx)');
            end
        end
        
        
        function [wij,wijind]=compute_downwash_vor2(obj,panel_idx,colloc_len)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] [0;0;0] [0;0;0] V*1E10];
            pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
            po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
            vortex=[obj.grid(:,obj.panels(4,panel_idx)) obj.grid(:,obj.panels(4,panel_idx)) pi po obj.grid(:,obj.panels(3,panel_idx)) obj.grid(:,obj.panels(3,panel_idx))]+obj.dvortex;
            wijind=zeros(1,colloc_len);
            wij=zeros(1,colloc_len);
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                wind=[0 0 0];
                for i=1:length(vortex)-1
                    r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                    r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                    r0=vortex(:,i+1)-vortex(:,i);
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB;
                        if i~=3
                            wind=wind+wAB;
                        end
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wijind(1,colloc_idx)=sum(conj(wind).*obj.colloc_nvec(:,panel_idx)');
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,panel_idx)');
            end
        end
        
        function [wij,wijind]=compute_downwash_vor3(obj,panel_idx,colloc_len)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0]  [0;0;0] [0;0;0] [0;0;0] V*1E10];
            pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
            po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
            vortex=[obj.grid(:,obj.te_idx(obj.panels(4,panel_idx))) obj.grid(:,obj.te_idx(obj.panels(4,panel_idx))) pi po obj.grid(:,obj.te_idx(obj.panels(3,panel_idx))) obj.grid(:,obj.te_idx(obj.panels(3,panel_idx)))]+obj.dvortex;
            %             hold on
            %             plot3(vortex(1,:),vortex(2,:),vortex(3,:),'r');
            
            wijind=zeros(1,colloc_len);
            wij=zeros(1,colloc_len);
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                wind=[0 0 0];
                for i=1:length(vortex)-1
                    r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                    r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                    r0=vortex(:,i+1)-vortex(:,i);
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    
                    coredist=norm((norm(cross(r1,r0)))/(norm(r0)));
                    
                    if norm((0.5*pi+0.5*po-obj.colloc(:,colloc_idx)))>8
                        coredist_lim=0.5;
                    else
                        coredist_lim=1e-8;
                    end
                    
                    
                    if (~isnan(wAB(1))) && (coredist>coredist_lim)
                        
                        w=w+wAB;
                        if (i==1)||(i==6)
                            wind=wind+wAB;
                        end
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wijind(1,colloc_idx)=sum(conj(wind).*obj.colloc_nvec(:,panel_idx)');
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,panel_idx)');
            end
        end
        
        
        
        function wij=compute_downwash_new(obj,panel_idx,colloc_len)
            vortex=obj.compute_vortex_coords(panel_idx);
            wij=zeros(1,colloc_len);
            
            for colloc_idx=1:1:colloc_len
                w=[0 0 0];
                for i=1:length(vortex)-1
                    r1=obj.colloc(:,colloc_idx)-vortex(:,i);
                    r2=obj.colloc(:,colloc_idx)-vortex(:,i+1);
                    r0=vortex(:,i+1)-vortex(:,i);
                    %r1xr2=cross(r1,r2);
                    r1xr2=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
                    r1dr2=(r1/norm(r1)-r2/norm(r2));
                    wAB=-r1xr2/norm(r1xr2)^2*(sum(conj(r0).*r1dr2));
                    %wAB=-r1xr2/norm(r1xr2)^2*(dot(r0,(r1/norm(r1)-r2/norm(r2)),1)/(4*3.141592));
                    if ~isnan(wAB(1))
                        w=w+wAB;
                    end
                end
                %wij=dot(w,obj.colloc_nvec(:,panel_idx));
                wij(1,colloc_idx)=sum(conj(w).*obj.colloc_nvec(:,panel_idx)');
            end
            %             if mod(colloc_len,10)==0
            %                 hold on
            %                 plot3(vortex(1,:),vortex(2,:),vortex(3,:),'LineWidth',3);
            %                 drawnow
            %             end
        end
        
        function obj=detect_luv_lee(obj)
            obj.flag=size(1,length(obj.panels));
            for i=1:length(obj.panels)
                costheta=dot(obj.Uinf,obj.colloc_nvec(:,i));
                if costheta>0
                    obj.flag(1,i)=1;
                else
                    obj.flag(1,i)=0; 
                end
            end
        end
        
        function wij=compute_downwash_mex(obj,panel_idx)
            V=obj.Uinf';
            obj.dvortex=[V*1E10 [0;0;0] [0;0;0] V*1E10];
            obj.grid(:,obj.panels(:,panel_idx))
            pi=0.75*obj.grid(:,obj.panels(1,panel_idx))+0.25*obj.grid(:,obj.panels(4,panel_idx));
            po=0.75*obj.grid(:,obj.panels(2,panel_idx))+0.25*obj.grid(:,obj.panels(3,panel_idx));
            vortex=[pi pi po po]+obj.dvortex;
            wij=compute_downwash(vortex(1:3,1:4),obj.colloc(:,:),obj.colloc_nvec(:,:));
        end
        
        function obj=compute_influence_coeff_matrix(obj)
            %obj.plot_flag;
            n=length(obj.colloc);
            C=zeros(n,n);
            Cind=zeros(n,n);
            for i=1:n
                
                [C(1:n,i),Cind(1:n,i) ]=obj.compute_downwash_ad(i,n);
            end
            obj.C=C;
            obj.Cind=Cind;
            %                         for i=1:n
            %                             for j=1:n
            %                                  C2(i,j)=obj.compute_downwash(i,j);
            %                             end
            %                         end
            %                         obj.C=C2;
            
        end
        
        function obj=compute_influence_coeff_matrix_vor2(obj)
            %obj.plot_flag;
            n=length(obj.colloc);
            C=zeros(n,n);
            Cind=zeros(n,n);
            for i=1:n
                [C(1:n,i),Cind(1:n,i) ]=obj.compute_downwash_vor2(i,n);
            end
            obj.C=C;
            obj.Cind=Cind;
        end
        
        function obj=compute_influence_coeff_matrix_vor3(obj)
            %obj.plot_flag;
            n=length(obj.colloc);
            C=zeros(n,n);
            Cind=zeros(n,n);
            for i=1:n
                [C(1:n,i),Cind(1:n,i) ]=obj.compute_downwash_vor3(i,n);
            end
            obj.C=C;
            obj.Cind=Cind;
        end
        
        function obj=determine_boundary_conditions(obj)
            for i=1:1:length(obj.colloc)
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf)*(4*pi);
            end
        end
        
        function obj=determine_boundary_conditions_deflection_pqr(obj,p,q,r,deflection_induced_speed)
            for i=1:length(obj.colloc)
                dV=cross(obj.r(:,i),[p,q,r]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV+deflection_induced_speed(:,i)')*(4*pi);
            end
        end
        
        function obj=determine_boundary_conditions_q(obj)
            dpqr=deg2rad(0.001);
            %dq=0.05*obj.Ma_corr; %% TODO: check if correct?
            dq=dpqr*obj.Ma_corr; %% TODO: check if correct?
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[0 dq 0]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV)*(4*pi);
            end
        end
        
        function obj=determine_boundary_conditions_p(obj)
            dpqr=deg2rad(0.001);
            %dp=0.05;
            dp=dpqr*obj.Ma_corr; %% TODO: check if correct?
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[dp 0 0]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV)*(4*pi);
            end
        end
        
        function obj=determine_boundary_conditions_r(obj)
            dpqr=deg2rad(0.001);
            %dr=0.05;
            dr=dpqr*obj.Ma_corr; %% TODO: check if correct?
            for i=1:1:length(obj.colloc)
                dV=cross(obj.r_dwn(:,i),[0 0 dr]);
                obj.b(i)=-dot(obj.colloc_nvec(:,i),obj.Uinf+dV)*(4*pi);
            end
        end
        
        function obj=solve(obj)
            obj=obj.compute_influence_coeff_matrix();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        
        function obj=solve_vor2(obj)
            obj=obj.compute_influence_coeff_matrix_vor2();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        
        function obj=solve_vor3(obj)
            obj=obj.compute_influence_coeff_matrix_vor3();
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        
        
        function obj=solve_MEX_vor3(obj)
            [obj.C,obj.Cind]=compute_influence_coefficients_vor3(obj.grid,obj.panels,obj.te_idx',obj.colloc(:,:),obj.colloc_nvec(:,:),obj.Uinf);
            obj.Cind=[zeros(size(obj.C)),zeros(size(obj.C)),zeros(size(obj.C))];
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        
        function obj=solve_MEX_vor4(obj)
            [obj.C,cc1,cc2,cc3]=compute_influence_coefficients_vor4(obj.grid,obj.panels,obj.te_idx',obj.colloc(:,:),obj.colloc_nvec(:,:),obj.Uinf);
            obj.Cind=[cc1,cc2,cc3];
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        function obj=solve_MEX_vor5(obj) % in vor5, the Cind Matrix contains the influence coefficients on the 1/4 point of each panel
            [obj.C,cc1,cc2,cc3]=compute_influence_coefficients_vor5(obj.grid,obj.panels,obj.te_idx',obj.colloc(:,:),obj.colloc_nvec(:,:),obj.Uinf);
            obj.Cind=[cc1,cc2,cc3];
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
        end
        function obj=f_postprocess(obj,CD_f)
        % as second input either 'quick' for no
            alpha=atan(obj.Uinf(3)/obj.Uinf(1));
            beta=atan(obj.Uinf(2)/obj.Uinf(1));
            if obj.Ma_corr<1
                beta_inf=obj.Ma_corr;
            else
                beta_inf=1;
            end
            M_BA=[  cos(alpha)*cos(beta) , -cos(alpha)*sin(beta), -sin(alpha);
                sin(beta)             , cos(beta)            ,    0;
                sin(alpha)*cos(beta) , -sin(alpha)*sin(beta), cos(alpha);];
            
			%determine downwash at the one quarter point
            nPanels = size(obj.panels,2);
            u=(obj.Cind(:,1:nPanels)*obj.Gamma)/(4*pi);
            v=obj.Cind(:,nPanels+1:2*nPanels)*obj.Gamma/(4*pi);
            w=obj.Cind(:,2*nPanels+1:3*nPanels)*obj.Gamma/(4*pi);





            for i=1:1:length(obj.Gamma)

                p_i=0.75*obj.grid(:,obj.panels(1,i))+0.25*obj.grid(:,obj.panels(4,i));
                p_o=0.75*obj.grid(:,obj.panels(2,i))+0.25*obj.grid(:,obj.panels(3,i));
                sp=p_o-p_i;

                %obj.crossp(:,i)=(obj.rho*obj.Gamma(i)*cross(sp',obj.Uinf'))';
                
%                 obj.Uinf(3)=obj.Uinf(3)+uvwind(3); 
%                 obj.Uinf(2)=obj.Uinf(2)+uvwind(2);
%                 obj.Uinf(1)=obj.Uinf(1)+uvwind(1);
%                 
                obj.F_body(:,i)=1/beta_inf*(obj.rho*obj.Gamma(i)*cross(sp',(obj.Uinf+[u(i),v(i),w(i)])'))';
                obj.M_body(:,i)=zeros(3,1);
               % obj.F_body(2,i)=-obj.F_body(2,i);
                
                %alpha2=atan((obj.Uinf(3)+uvwind(3))/(obj.Uinf(1)+uvwind(1)));
                %beta=atan((obj.Uinf(2)+uvwind(2))/(obj.Uinf(1)+uvwind(1)));
                
                alpha=atan((obj.Uinf(3))/(obj.Uinf(1)));
                beta=atan((obj.Uinf(2))/(obj.Uinf(1)));
                %                alphai(i)=alpha-alpha2;
                
                M_BA=[  cos(alpha)*cos(beta) , -cos(alpha)*sin(beta), -sin(alpha);
                    sin(beta)             , cos(beta)            ,    0;
                    sin(alpha)*cos(beta) , -sin(alpha)*sin(beta), cos(alpha);];
                obj.F_aero(:,i)=M_BA'*obj.F_body(:,i);
                % obj.F_aero(:,i)=(obj.rho*obj.Gamma(i)*cross(sp',obj.Uinf'))';
                % obj.cdi(:,i)=-(obj.rho*obj.Gamma(i)*obj.wind(i)*cross(sp',obj.Uinf'/norm(obj.Uinf)))';
                % obj.cdi(:,i)=M_BA'*(obj.rho*obj.Gamma(i)*cross(sp',(obj.Uinf'+uvwind)/norm(obj.Uinf+uvwind')))';
                
                obj.cp(i)=1/beta_inf*2*obj.Gamma(i)/(norm(obj.Uinf)*norm(obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i))));
                
                %obj.A(i)=norm(cross(obj.grid(:,obj.panels(2,i))-obj.grid(:,obj.panels(1,i)),obj.grid(:,obj.panels(4,i))-obj.grid(:,obj.panels(1,i))));
                obj.cl(:,i)=obj.F_body(:,i)/(obj.qinf*obj.reference.S_ref);
            end
            
            %obj.Clq=sum(obj.r(1,:)*obj.A(i))/(obj.reference.S_ref);
            
            obj.cx=obj.F_body(1,:)/(obj.qinf*obj.reference.S_ref);
            obj.cy=obj.F_body(2,:)/(obj.qinf*obj.reference.S_ref);
            obj.cz=obj.F_body(3,:)/(obj.qinf*obj.reference.S_ref);
            
            obj.CX=sum(obj.F_body(1,:))/(obj.qinf*obj.reference.S_ref);
            obj.CY=sum(obj.F_body(2,:))/(obj.qinf*obj.reference.S_ref);
            obj.CZ=sum(obj.F_body(3,:))/(obj.qinf*obj.reference.S_ref);
            
            obj.Cl=sum(obj.F_aero(3,:))/(obj.qinf*obj.reference.S_ref);
            obj.Cy=sum(obj.F_aero(2,:))/(obj.qinf*obj.reference.S_ref);
            obj.Cdi=sum(obj.F_aero(1,:))/(obj.qinf*obj.reference.S_ref);

                        
            C_CD_f=M_BA*[CD_f; 0; 0];
            
            obj.CX=obj.CX+C_CD_f(1);
            obj.CY=obj.CY+C_CD_f(2);
            obj.CZ=obj.CZ+C_CD_f(3);
            
            obj.CL=sum(+obj.F_body(3,:).*obj.r(2,:)-obj.F_body(2,:).*obj.r(3,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
            obj.CM=sum(-obj.F_body(3,:).*obj.r(1,:)+obj.F_body(1,:).*obj.r(3,:))/(obj.qinf*obj.reference.S_ref*obj.reference.c_ref);
            obj.CN=sum(obj.F_body(2,:).*obj.r(1,:)-obj.F_body(1,:).*obj.r(2,:))/(obj.qinf*obj.reference.S_ref*obj.reference.b_ref);
            
        end
        
        
        
        
        function grid=get_grid(obj)
            grid(1,:)=obj.grid(1,:)*obj.Ma_corr;
            grid(2,:)=obj.grid(2,:);
            grid(3,:)=obj.grid(3,:);
        end
        function obj=f_solve_std(obj)
            obj=obj.solve_MEX_vor5();
            %obj=obj.solve();
            obj=obj.f_postprocess(0);
        end
        
        function obj=f_solve_dynamic(obj,p,q,r,deflection_induced_speed)
            obj=obj.solve_MEX_vor3_dynamic(p,q,r,deflection_induced_speed);
            %obj=obj.solve();
            obj=obj.f_postprocess('quick');
        end
        
        function obj=f_solve_vor3(obj)
            %vor 3 -> no computation of Cind ->Lift normal to panel and
            %Uinf
            obj=obj.solve_MEX_vor3(); 
            obj=obj.f_postprocess(0);
        end
                
        function obj=f_solve_vor4(obj)
            %vor 4 -> Cind for colloc points
            obj=obj.solve_MEX_vor4(); 
            obj=obj.f_postprocess(0);
        end
        
        function obj=f_solve_vor5(obj)
            %vor 4 -> Cind for 1/4 points 
            obj=obj.solve_MEX_vor5(); 
            obj=obj.f_postprocess(0);
        end
        
        function obj=f_solve_full(obj,varargin)
            
            if nargin==2
                CD_f=varargin{1};
            else
                CD_f=0;
            end
            
            dalpha=deg2rad(0.001);
            dbeta=deg2rad(0.001);
            dpqr=deg2rad(0.001);
            du=0.1;
            origUinf=obj.Uinf;
            alpha=obj.state.alpha;
            beta=obj.state.beta;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,beta,obj.Ma_corr);
            obj.Uinf=obj.Uinf+[du 0 0];
            obj=obj.determine_boundary_conditions();
            %Wenn hier nicht Uinf wieder r�ckg�ngig korrigiert wird
            %beinhalten alle folgenden Derivative den Einfluss der Uinf
            %Perturbation 
            obj.Uinf=origUinf;

            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXdu=obj.CX;
            CYdu=obj.CY;
            CZdu=obj.CZ;
            Cldu=obj.Cl;
            Cdidu=obj.Cdi;
            Cydu=obj.Cy;
            CLdu=obj.CL;
            CMdu=obj.CM;
            CNdu=obj.CN;
            
            alpha=obj.state.alpha+dalpha;
            beta=obj.state.beta;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,beta,obj.Ma_corr);
            obj=obj.determine_boundary_conditions();
            obj.Uinf=origUinf;
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXda=obj.CX;
            CYda=obj.CY;
            CZda=obj.CZ;
            Clda=obj.Cl;
            Cdida=obj.Cdi;
            Cyda=obj.Cy;
            CLda=obj.CL;
            CMda=obj.CM;
            CNda=obj.CN;
            
            alpha=obj.state.alpha;
            beta=obj.state.beta+dbeta;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,beta,obj.Ma_corr);
            obj=obj.determine_boundary_conditions();
            obj.Uinf=origUinf;
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXdb=obj.CX;
            CYdb=obj.CY;
            CZdb=obj.CZ; 
            Cldb=obj.Cl;
            Cdidb=obj.Cdi;
            Cydb=obj.Cy;
            CLdb=obj.CL;
            CMdb=obj.CM;
            CNdb=obj.CN;
            
            alpha=obj.state.alpha;
            beta=obj.state.beta;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,beta,obj.Ma_corr);
            obj=obj.solve_MEX_vor3();
            obj=obj.f_postprocess(CD_f);
            
            
            obj=obj.determine_boundary_conditions_p();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXdp=obj.CX;
            CYdp=obj.CY;
            CZdp=obj.CZ;
            Cldp=obj.Cl;
            Cdidp=obj.Cdi;
            Cydp=obj.Cy;
            CLdp=obj.CL;
            CMdp=obj.CM;
            CNdp=obj.CN;
            
            obj=obj.determine_boundary_conditions_q();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXdq=obj.CX;
            CYdq=obj.CY;
            CZdq=obj.CZ;
            Cldq=obj.Cl;
            Cdidq=obj.Cdi;
            Cydq=obj.Cy;
            CLdq=obj.CL;
            CMdq=obj.CM;
            CNdq=obj.CN;
            
            obj=obj.determine_boundary_conditions_r();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            CXdr=obj.CX;
            CYdr=obj.CY;
            CZdr=obj.CZ;
            Cldr=obj.Cl;
            Cdidr=obj.Cdi;
            Cydr=obj.Cy;
            CLdr=obj.CL;
            CMdr=obj.CM;
            CNdr=obj.CN;
            
            obj=obj.determine_boundary_conditions();
            obj.Gamma=linsolve(obj.C,obj.b');
            %obj.wind=obj.Cind*obj.Gamma/(4*pi);
            obj=obj.f_postprocess(CD_f);
            
            %> u derivatives
            obj.CXu=(CXdu-obj.CX)/(du/norm(obj.Uinf));
            obj.CYu=(CYdu-obj.CY)/(du/norm(obj.Uinf));
            obj.CZu=(CZdu-obj.CZ)/(du/norm(obj.Uinf));
            obj.Clu=(Cldu-obj.Cl)/(du/norm(obj.Uinf));
            obj.Cyu=(Cydu-obj.Cy)/(du/norm(obj.Uinf));
            obj.Cdiu=(Cdidu-obj.Cdi)/(du/norm(obj.Uinf));
            obj.CLu=(CLdu-obj.CL)/(du/norm(obj.Uinf));
            obj.CMu=(CMdu-obj.CM)/(du/norm(obj.Uinf));
            obj.CNu=(CNdu-obj.CN)/(du/norm(obj.Uinf));
            %> Alpha Derivatives in body axis
            obj.CXa=(CXda-obj.CX)/(dalpha*pi/180);
            obj.CYa=(CYda-obj.CY)/(dalpha*pi/180);
            obj.CZa=(CZda-obj.CZ)/(dalpha*pi/180);
            obj.Cla=(Clda-obj.Cl)/(dalpha*pi/180);
            obj.Cya=(Cyda-obj.Cy)/(dalpha*pi/180);
            obj.Cdia=(Cdida-obj.Cdi)/(dalpha*pi/180);
            obj.CLa=(CLda-obj.CL)/(dalpha*pi/180);
            obj.CMa=(CMda-obj.CM)/(dalpha*pi/180);
            obj.CNa=(CNda-obj.CN)/(dalpha*pi/180);
            %> Beta Derivatives in body axis
            obj.CXb=(CXdb-obj.CX)/(dbeta*pi/180);
            obj.CYb=(CYdb-obj.CY)/(dbeta*pi/180);
            obj.CZb=(CZdb-obj.CZ)/(dbeta*pi/180);
            obj.Clb=(Cldb-obj.Cl)/(dbeta*pi/180);
            obj.Cyb=(Cydb-obj.Cy)/(dbeta*pi/180);
            obj.Cdib=(Cdidb-obj.Cdi)/(dbeta*pi/180);
            obj.CLb=(CLdb-obj.CL)/(dbeta*pi/180);
            obj.CMb=(CMdb-obj.CM)/(dbeta*pi/180);
            obj.CNb=(CNdb-obj.CN)/(dbeta*pi/180);
            %> p Derivatives in body axis
            obj.CXp=(CXdp-obj.CX)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CYp=(CYdp-obj.CY)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CZp=(CZdp-obj.CZ)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Clp=(Cldp-obj.Cl)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Cyp=(Cydp-obj.Cy)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Cdip=0;
            obj.CLp=(CLdp-obj.CL)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CMp=(CMdp-obj.CM)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CNp=(CNdp-obj.CN)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            %> q Derivatives in body axis
            obj.CXq=(CXdq-obj.CX)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.CYq=(CYdq-obj.CY)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.CZq=(CZdq-obj.CZ)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.Clq=(Cldq-obj.Cl)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.Cyq=(Cydq-obj.Cy)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.Cdiq=0;
            obj.CLq=(CLdq-obj.CL)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.CMq=(CMdq-obj.CM)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            obj.CNq=(CNdq-obj.CN)/(dpqr*obj.reference.c_ref/(2*norm(obj.Uinf)));
            %> r Derivatives in body axis
            obj.CXr=(CXdr-obj.CX)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CYr=(CYdr-obj.CY)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CZr=(CZdr-obj.CZ)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Clr=(Cldr-obj.Cl)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Cyr=(Cydr-obj.Cy)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.Cdir=0;
            obj.CLr=(CLdr-obj.CL)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CMr=(CMdr-obj.CM)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
            obj.CNr=(CNdr-obj.CN)/(dpqr*obj.reference.b_ref/(2*norm(obj.Uinf)));
        end
        
        function obj=f_solve_exp(obj)
            obj=obj.solve();
            obj=obj.f_postprocess();
        end
        
        function obj=f_solve_for_Cl(obj,Cl_target)
            %% TODO: validate, enable entering of alpha start and tolerance
            alpha_0=0;
            tolerance=0.0001;
            
            
            obj.Uinf=a2bf(norm(obj.Uinf),alpha_0,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve_std();
            Cl_1=obj.CZ;
            alpha_0=5;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha_0,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve_std();
            Cl_prev=obj.CZ;
            dCLalpha=Cl_prev-Cl_1;
            dalpha=(Cl_target-Cl_prev)/dCLalpha;
            alpha=alpha_0+dalpha;
            obj.Uinf=a2bf(norm(obj.Uinf),alpha,obj.state.beta,obj.Ma_corr);
            obj=obj.f_solve_std();
            dCLalpha=(obj.CZ-Cl_prev)/dalpha;
            i=1;
            while abs(obj.CZ-Cl_target)>tolerance
                Cl_prev=obj.CZ;
                dalpha=(Cl_target-Cl_prev)/dCLalpha;
                alpha=alpha+dalpha;
                obj.Uinf=a2bf(norm(obj.Uinf),alpha,obj.state.beta,obj.Ma_corr);
                obj=obj.f_solve_std();
                i=i+1;
                dCLalpha=(obj.CZ-Cl_prev)/dalpha;
                if i>30
                    sprintf('f_solve_for_Cl did not converge!')
                    break;
                end
            end
            obj=obj.f_solve_std();
            obj.state=obj.state.set_alpha(alpha);
        end
        
        function obj=plot_grid(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(squeeze(obj.panels(i,:)),1), obj.grid(squeeze(obj.panels(i,:)),2),obj.grid(squeeze(obj.panels(i,:)),3),'b');
                alpha(handle,0.4);
            end
        end
        
        function obj=plot_L(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'r');
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
        
        function obj=plot_cx(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cx(i));
                hold on
            end
            axis equal
        end
               function obj=plot_Fy(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.F_body(2,i));
                hold on
            end
            axis equal
        end
        
        function obj=plot_cp(obj)
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i))*obj.Ma_corr, obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),obj.cp(i));
                hold on
            end
            axis equal
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
%             fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid_wake),length(obj.panels_wake));
%             
%             for i =1:length(obj.grid_wake)
%                 ri=obj.grid_wake(:,i)-obj.reference.p_ref';
%                 pt=M_BI'*ri+pos(1:3);
%                 fprintf(fileID,'%f %f %f %f\n',pt(1)*obj.Ma_corr,pt(2),pt(3),nodal_gamma_wake(i));
%             end
%             
%             for i =1:length(obj.panels_wake)
%                 fprintf(fileID,'%i %i %i %i \n',obj.panels_wake(1,i),obj.panels_wake(2,i),obj.panels_wake(3,i),obj.panels_wake(4,i));
%             end
            fclose(fileID);
        end
        function obj=update_nvec(obj, aircraft)
            
        %for i=all wings
        %temporarily only wing 1
            i=1;
                p=1;
                nPanels=size(aircraft.wings(i).panels,2);
                improvedCollocNvec(3,nPanels)=0;
                for j=1:size(aircraft.wings(i).wing_segments,2)
                    nSegSpan=aircraft.wings(i).wing_segments(j).n_span;
                    nSegChord=aircraft.wings(i).wing_segments(j).n_chord;
                    %sweep due to taper - twist is defined around c4 line
                    c4SweepDueTaper=atan((aircraft.wings(i).wing_segments(j).c_r-aircraft.wings(i).wing_segments(j).c_t)/(4*aircraft.wings(i).wing_segments(j).b));

                    for iSegSpan=1:nSegSpan
                        %spanwise relative position
                        relativeSpanwisePosition=((iSegSpan-1)+0.5)/nSegSpan;
                        %twist
                        localTwist=-((1-relativeSpanwisePosition)*aircraft.wings(i).wing_segments(j).Theta_r+relativeSpanwisePosition*aircraft.wings(i).wing_segments(j).Theta_t)*pi/180;
                        for iSegChord=1:nSegChord
                        %chordwise relative position position
                        chordwise_relative_position=((iSegChord-1)+0.75)/nSegChord;
                        %profile camber
                        skeleton_angle=aircraft.wings(i).wing_segments(j).compute_skeleton_angle(chordwise_relative_position,relativeSpanwisePosition);
                        improvedCollocNvec(:,p)=angle2dcm(localTwist,-c4SweepDueTaper,skeleton_angle,'YZY')*[0 0 -1]';
                        p=p+1;
                        end
                    end
                    %dihedral
                    improvedCollocNvec=angle2dcm(0,0,-aircraft.wings(i).wing_segments(j).dihed*pi/180)*improvedCollocNvec;
                end
                if aircraft.wings(i).symmetric
                    improvedCollocNvec(1,end/2+1:end)=improvedCollocNvec(1,1:end/2);
                    improvedCollocNvec(2,end/2+1:end)=-improvedCollocNvec(2,1:end/2);
                    improvedCollocNvec(3,end/2+1:end)=improvedCollocNvec(3,1:end/2);
                end
                obj.colloc_nvec(:,aircraft.wings(i).panel_start_idx:aircraft.wings(i).panel_start_idx+nPanels-1)=improvedCollocNvec;
            %end
                
        end
        
        function obj=update_nvec2(obj, aircraft)
            if obj.Ma_corr<1
                grid_flat(1,:)=aircraft.grid_flat(1,:)/obj.Ma_corr;
                grid_flat(2,:)=aircraft.grid_flat(2,:);
                grid_flat(3,:)=aircraft.grid_flat(3,:);
            else
                grid_flat(1,:)=aircraft.grid_flat(1,:);
                grid_flat(2,:)=aircraft.grid_flat(2,:);
                grid_flat(3,:)=aircraft.grid_flat(3,:);                
            end
            for iPanel=1:length(obj.panels)
                %compute normal vectors from grid flat as in compute_colloc_points     
                obj.colloc_nvec(:,iPanel)=cross(grid_flat(:,aircraft.panels(3,iPanel))-grid_flat(:,aircraft.panels(1,iPanel)),grid_flat(:,aircraft.panels(4,iPanel))-grid_flat(:,aircraft.panels(2,iPanel)));
                obj.colloc_nvec(:,iPanel)=obj.colloc_nvec(:,iPanel)./norm(obj.colloc_nvec(:,iPanel));
            end
            %rotate flat grid colloc_nvec around global y axis with local skeleton angle
            currentPanel=1;
            improvedCollocNvec(3,length(obj.panels))=0;
            for iWings=1:size(aircraft.wings,2)
                for iSegments=1:size(aircraft.wings(iWings).wing_segments,2)
                    nPanelsSpanwise=aircraft.wings(iWings).wing_segments(iSegments).n_span;
                    nPanelsChordwise=aircraft.wings(iWings).wing_segments(iSegments).n_chord;
                    if aircraft.wings(iWings).wing_segments(iSegments).has_le_cs
                        nPanelsChordwise=nPanelsChordwise+aircraft.wings(iWings).wing_segments(iSegments).n_le_panels;
                    end
                    if aircraft.wings(iWings).wing_segments(iSegments).has_te_cs
                        nPanelsChordwise=nPanelsChordwise+aircraft.wings(iWings).wing_segments(iSegments).n_te_panels;
                    end
                    %sweep due to taper - twist is defined around c4 line
                    for iPanelsSpanwise=1:nPanelsSpanwise
                        %spanwise relative position
                        relativeSpanwisePosition=((iPanelsSpanwise-1)+0.5)/nPanelsSpanwise;
                        for iPanelChordwise=1:nPanelsChordwise
                            %chordwise relative position position
                            chordwise_relative_position=((iPanelChordwise-1)+0.75)/nPanelsChordwise;
                            chordwise_relative_position_0=((iPanelChordwise-1)+0)/nPanelsChordwise;
                            chordwise_relative_position_1=((iPanelChordwise-1)+1)/nPanelsChordwise;
                            %profile camber
                            skeletonAngle=aircraft.wings(iWings).wing_segments(iSegments).compute_skeleton_angle(chordwise_relative_position,chordwise_relative_position_0,chordwise_relative_position_1,relativeSpanwisePosition);
                            improvedCollocNvec(:,currentPanel)=angle2dcm(0,0,skeletonAngle,'YZY')*obj.colloc_nvec(:,currentPanel);
                            currentPanel=currentPanel+1;
                        end
                    end
                end
                if aircraft.wings(iWings).symmetric
                    nWingPanels=size(aircraft.wings(iWings).panels,2)/2;
                    improvedCollocNvec(1,currentPanel:currentPanel+nWingPanels-1)=improvedCollocNvec(1,currentPanel-nWingPanels:currentPanel-1);
                    improvedCollocNvec(2,currentPanel:currentPanel+nWingPanels-1)=-improvedCollocNvec(2,currentPanel-nWingPanels:currentPanel-1);
                    improvedCollocNvec(3,currentPanel:currentPanel+nWingPanels-1)=improvedCollocNvec(3,currentPanel-nWingPanels:currentPanel-1);
                    currentPanel=currentPanel+nWingPanels;
                end
            end
            %overwrite colloc_nvec
            obj.colloc_nvec=improvedCollocNvec;
        end
    end
    
end

