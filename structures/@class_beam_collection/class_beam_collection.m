%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_beam_collection.m
%> @brief File containing the class for a finite element beam collection
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
%> @brief class for a finite element beam collection
%> class_beam_collection: couple and solve systems of beams
%> This class collects several classes of class_beam and derivates for
%> example class_wing and has the capability of coupling beams together
%> with certain coupling/boundary conditions
%> class can solve the coupled system equations and updates the
%> subclasses with respective information
%> this class has the same public methods as a single beam
% ======================================================================

classdef class_beam_collection
    
    properties
        %> name of beam collection
        identifier;
        %> number of beams in collection
        n;
        %> array of beams
        beam;
        %> number of coupling conditions
        nc=0;
        %> coupling conditions
        coupling_condition;
        
        sort_vec;
        deflections;
        dof_node_beam;
        Ktest;
        Mtest; % Mass Matrix test
        
        Mff_lumped;
        Ftest;
        Ptest;
        
        Mff;
        
        Kff;
        
        n_dof_sys;
        
        Pglob;
        Kglob;
        %> for consistent mass matrix
        Mglob;
        
        %> lumped mass matrix
        Mglob_lumped; %
        Fglob;
        m_total=0;
        
        def;
        
        nodal_deflections;
        modal_deflections;
        
        J;
        
        err=100;
        
        T_glob;
        db_glob;
         
        node_coords;
        node_coords_full;
        node_orient;
        nodal_deflections_ordered;
        
        %% dynamic solution results
        modeshapes;
        modefrequencies;
        
        modeshapes_lumped;
        modefrequencies_lumped;
        
        %% settings
        
        settings
        %> index of loadcase
        loadcase_idx
    end
    
    methods
        %> @brief constructor
        function obj = class_beam_collection(beam1,varargin)
            obj.beam=beam1;
            obj.n=1;
            
            if ~isempty(varargin)
                obj.identifier=cell2mat(varargin(1));
            end
        end
        
        function obj=f_set_state(obj,state)
            for i=1:length(obj.beam)
                obj.beam(i).load_factor=state.load_factor;
                obj.beam(i).loadcase_index=state.loadcase_index;
            end
        end
        
        function obj=f_set_solver_settings(obj,solver_settings)
            obj.settings=solver_settings;
            for i=1:length(obj.beam)
                obj.beam(i).settings=solver_settings;
            end
        end
                
        obj=f_assemble(obj,loadstep,eval_nonlin);
        
        obj=f_assemble_free(obj,loadstep,eval_nonlin);
        %> @brief deep copy of this class (since sub beams are handles)
        function new = f_copy(this)
            % Instantiate new object of the same class.
            new=this;
            for i=1:length(this.beam)
                new.beam(i)=this.beam(i).f_copy;
            end
        end
        
        %> @brief add beam
        function obj =f_add_beam(obj,beam)
            obj.n=obj.n+1;
            obj.beam(obj.n)=beam;
        end
        
        function name=f_get_name(obj)
            if ~isempty(obj.identifier)
                name=obj.identifier;
            else
                name='identifier not set';
            end
        end
        
        
        function mass=f_get_mass(obj)
            mass=0;
            for i=1:length(obj.beam)
                for j=1:size(obj.beam(i).node_coords,1)-1
                    obj.beam(i).f_calc_mass(obj.weights)
                    mass=mass+obj.beam(i).m_total;
                end
            end  
        end
        
        
        function obj =f_add_coupling_condition(obj,beam_idx,node_idx,dof)
            obj.nc=obj.nc+1;
            
            if(obj.nc==1)
                obj.coupling_condition=class_coupling_condition(node_idx,beam_idx,dof);
            else
                obj.coupling_condition(obj.nc)=class_coupling_condition(node_idx,beam_idx,dof);
            end
            
        end
        
        function mass_structure=f_get_masses(obj,weights)
            n_fus=1;
            OEW=0;
            fuel=0;
            payload=0;
            MTOW=0;
            MZFW=0;
            disp([obj.identifier '         has a MTOW of            ' num2str(obj.m_total) 'kg'])
            for i=1:length(obj.beam)
              disp([obj.beam(i).identifier '  has a mass of         ' num2str(obj.beam(i).m_total) 'kg'])
              if isa(obj.beam(i),'class_wing')
                  disp(['            Structural Mass                ' num2str(obj.beam(i).m_total) 'kg'])   
                  OEW=OEW+obj.beam(i).m_total-obj.beam(i).m_fuel_total;
                  disp(['                   Wingbox Mass            ' num2str(obj.beam(i).m_wingbox_total) 'kg'])
                  disp(['                   Other Mass              ' num2str(-obj.beam(i).m_wingbox_total+obj.beam(i).m_total-obj.beam(i).m_fuel_total) 'kg'])
                  disp(['            NonStructural Mass             ' num2str(obj.beam(i).m_fuel_total) 'kg'])
                  disp(['                   Fuel Mass               ' num2str(obj.beam(i).m_fuel_total) 'kg'])
                  fuel=fuel+obj.beam(i).m_fuel_total;
                  disp(['                   Fuel Volume             ' num2str(obj.beam(i).fuel_volume_total) 'm3'])
              elseif isa(obj.beam(i),'class_fuselage')
                  disp(['            Structural Mass                ' num2str(obj.beam(i).m_total-weights.FuselageNonStructuralEstimate(n_fus)) 'kg'])
                  disp(['                  Loadcarrying Mass        ' num2str(obj.beam(i).m_shell_total) 'kg'])  
                  disp(['                  Other Mass               ' num2str(obj.beam(i).m_total-weights.FuselageNonStructuralEstimate(n_fus)-obj.beam(i).m_shell_total) 'kg']) 
                  disp(['            Nonstructural Mass             ' num2str(weights.FuselageNonStructuralEstimate(n_fus)) 'kg'])      
                  disp(['                  Passengers               ' num2str(obj.beam(i).n_PAX*weights.PassengerWeight) 'kg'])     
                  disp(['                  Baggage                  ' num2str(obj.beam(i).n_PAX*weights.BaggageWeight) 'kg'])
                  payload=payload+obj.beam(i).n_PAX*weights.PassengerWeight+obj.beam(i).n_PAX*weights.BaggageWeight;
                  disp(['                  Aicraft Seats            ' num2str(obj.beam(i).n_PAX*weights.SeatWeight) 'kg']) 
                  rest=weights.FuselageNonStructuralEstimate(n_fus)-obj.beam(i).n_PAX*weights.PassengerWeight-obj.beam(i).n_PAX*weights.BaggageWeight-obj.beam(i).n_PAX*weights.SeatWeight;
                  OEW=OEW+obj.beam(i).m_total-weights.FuselageNonStructuralEstimate(n_fus)+obj.beam(i).n_PAX*weights.SeatWeight+rest;
                  disp(['              Other Interior, crew equipment and accommodation accessories' num2str(rest) 'kg']) 
                  n_fus=n_fus+1;
              else
                  OEW=OEW+obj.beam(i).m_total;
              end
            end
            disp([obj.identifier '         has a OEW of             ' num2str(OEW) 'kg'])
            disp([obj.identifier '         has a Payload of         ' num2str(payload) 'kg'])
            disp([obj.identifier '         has a MZFW of            ' num2str(OEW+payload) 'kg'])
            disp([obj.identifier '         has a Fuel of            ' num2str(fuel) 'kg'])
            disp([obj.identifier '         has a MTOW of            ' num2str(OEW+payload+fuel) 'kg'])
            mass_structure=0;
        end
        
        %%
        % Calculates and returns the total mass of the structure (only the
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
        
        %%
        % Calculates the center of gravity and returns its coordinates in
        % the format [CGx;CGy;CGz].
        % Uses the mass matrix Mglob_lumped
        function CG = f_compute_CG(obj)
            nNodes = length(obj.node_coords);
            Ti = eye(3);
            D = [];
            Tr=zeros(3,3,nNodes);
            di = zeros(6,6,nNodes);
            CG = zeros(3,1);
            pos = zeros(size(obj.node_coords));
            
            for i=1:nNodes
                if isempty(obj.nodal_deflections)
                    pos(i,:) = obj.node_coords(i,:);
                else
                    pos(i,:) = obj.node_coords(i,:) + obj.nodal_deflections((i-1)*6+1:(i-1)*6+3)';
                end
            end
            for i=1:nNodes
                Tr(:,:,i)=[0,pos(i,3),-pos(i,2);-pos(i,3),0,pos(i,1);pos(i,2),-pos(i,1),0];
            end

            for i=1:nNodes
                di(:,:,i)=[Ti',Ti'*Tr(:,:,i);zeros(3,3),Ti'];
            end

            for i=1:nNodes
                D = [D, di(:,:,i)'];
            end
            D = D';

            MO = D'*obj.Mff_lumped*D;

            M_t = MO(1:3,1:3);
            M_tr = MO(1:3,4:6);

            delta = sqrt(M_t(1,1)^2+M_t(2,2)^2+M_t(3,3)^2);
            epsilon = sqrt(M_t(1,2)^2+M_t(1,3)^2+M_t(2,3)^2);
            if (epsilon/delta)>0.001
                disp('WARNING: Excessive coupling exists. CG calculation might be wrong')
            end

            M_x = M_t(1,1);
            M_y = M_t(2,2);

            CG(1,1) = M_tr(2,3)/M_y;
            CG(2,1) = -M_tr(1,3)/M_x;
            CG(3,1) = M_tr(1,2)/M_x;
        end
        
        %%
        % Calculates and returns the moment of inertia matrix.
        % Uses the mass matrix obj.Mff_lumped and needs the reference point as input.
        % All calculations are done in the body fixed system.
        function Inertia = f_compute_moment_of_inertia(obj, referencePoint)
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
                if isempty(obj.nodal_deflections)
                    delta_inertial(k:k+2)=0;
                    delta_inertial(k+3:k+5)=0;
                else
                    delta_inertial(k:k+2)=obj.nodal_deflections((j-1)*6+1:(j-1)*6+3);
                    delta_inertial(k+3:k+5)=obj.nodal_deflections((j-1)*6+4:j*6);
                end
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

        %% old version
%         function J=f_compute_moment_of_inertia(obj,pref)
%             J=zeros(1,3);
%             mTotal=0;
%             for i=1:length(obj.beam)
% %                 for j=1:size(obj.beam(i).node_coords,1)-1 01718273305
% %                     mi=obj.beam(i).beamelement(j).m;
% %                     le=obj.beam(i).beamelement(j).le;
% %                     if  isa(obj.beam(i),'class_wing')
% %                     mi=mi+obj.beam(i).beamelement(j).el_m_fuel/le;
% %                     end
% %                     ri=pref-(0.5*obj.beam(i).node_coords(j+1,:)+0.5*obj.beam(i).node_coords(j,:)+obj.beam(i).nodal_deflections(1+(i-1)*6:3+(i-1)*6)');
% %                     % TODO: check!
% % %                     J(1)=J(1)+sqrt(ri(3)^2+ri(2)^2)^2*mi*le;
% % %                     J(2)=J(2)+sqrt(ri(1)^2+ri(3)^2)^2*mi*le;
% % %                     J(3)=J(3)+sqrt(ri(1)^2+ri(2)^2)^2*mi*le;
% %                     
% %                     J(1)=J(1)+(ri(3)^2+ri(2)^2)*mi*le;
% %                     J(2)=J(2)+(ri(1)^2+ri(3)^2)*mi*le;
% %                     J(3)=J(3)+(ri(1)^2+ri(2)^2)*mi*le;
% %                 end
%                 
%                 for j=1:size(obj.beam(i).node_coords,1)-1
%                     mi=obj.beam(i).beamelement(j).m;
%                     le=obj.beam(i).beamelement(j).le;
%                     if  isa(obj.beam(i),'class_wing')
%                     mi=mi+obj.beam(i).beamelement(j).el_m_fuel/le;
%                     end
%                     ri=pref-(0.5*obj.beam(i).node_coords(j+1,:)+0.5*obj.beam(i).node_coords(j,:)+(0.5*obj.beam(i).nodal_deflections(1+(j-1)*6:3+(j-1)*6)+0.5*obj.beam(i).nodal_deflections(1+(j)*6:3+(j)*6))');
%                     % TODO: check!
% %                     J(1)=J(1)+sqrt(ri(3)^2+ri(2)^2)^2*mi*le;
% %                     J(2)=J(2)+sqrt(ri(1)^2+ri(3)^2)^2*mi*le;
% %                     J(3)=J(3)+sqrt(ri(1)^2+ri(2)^2)^2*mi*le;
%                     mTotal=mTotal+mi*le;
%                     J(1)=J(1)+(ri(3)^2+ri(2)^2)*mi*le+obj.beam(i).beamelement(j).elM_lumped_global(4,4)+obj.beam(i).beamelement(j).elM_lumped_global(10,10);
%                     J(2)=J(2)+(ri(1)^2+ri(3)^2)*mi*le+obj.beam(i).beamelement(j).elM_lumped_global(5,5)+obj.beam(i).beamelement(j).elM_lumped_global(11,11);
%                     J(3)=J(3)+(ri(1)^2+ri(2)^2)*mi*le+obj.beam(i).beamelement(j).elM_lumped_global(6,6)+obj.beam(i).beamelement(j).elM_lumped_global(12,12);
%                 end
%                 
%             end  
%             for i=1:3
%                 if J(i)<=1
%                     J(i)=1;
%                 end
%             end
%             J=diag(J);
%         end
        
%%
        function obj=f_set_acceleration(obj,acc,pref)
            for i=1:length(obj.beam)
                % faster assignment
                [obj.beam(i).beamelement.ax]=deal(acc(1));
                [obj.beam(i).beamelement.ay]=deal(acc(2));
                [obj.beam(i).beamelement.az]=deal(acc(3));
%                 for j=1:size(obj.beam(i).beamelement,2)
%                     %ri=(0.5*obj.beam(i).node_coords(j+1,:)+0.5*obj.beam(i).node_coords(j,:))-pref;
%                     obj.beam(i).beamelement(j).ax=acc(1);%-sin(ri(3))*acc(5)+sin(ri(2))*acc(6);
%                     obj.beam(i).beamelement(j).ay=acc(2);%+sin(ri(3))*acc(4)-sin(ri(1))*acc(6);
%                     obj.beam(i).beamelement(j).az=acc(3);%-sin(ri(2))*acc(4)+sin(ri(1))*acc(5);
%                 end
                obj.beam(i).update_Q=1;
            end
        end
        
        function deflections=f_get_deflections(obj)
            for i=1:length(obj.beam)
                if  isa(obj.beam(i),'class_wing')
                    mem=obj.beam(i).f_get_deflections_c4();
                    deflections(i).def=mem.def;
                elseif isa (obj.beam(i),'class_fuselage')
                    mem=obj.beam(i).f_get_deflections();
                    deflections(i).def=mem.def;
                elseif  isa (obj.beam(i),'class_pylon')
                    mem=obj.beam(i).f_get_deflections();
                    deflections(i).def=mem.def;
                end
            end
        end
        
        function structmesh=get_struct_mesh(obj)
            for i=1:length(obj.beam)
                mem=obj.beam(i).get_struct_mesh;
                structmesh(i).y=mem.y;
            end
        end
        
        function obj=plot_geometry(obj,varargin)
            for i=1:length(obj.beam)
                if~isempty(varargin)
                    obj.beam(i).plot_geometry(varargin);
                else
                    obj.beam(i).plot_geometry;
                end
            end
        end
        
        function obj=plot_deformations(obj,varargin)
            for i=1:length(obj.beam)
                if~isempty(varargin)
                    obj.beam(i).plot_deformations(varargin);
                else
                    obj.beam(i).plot_deformations;
                end
            end
        end
        
        function obj=plot_internalforces(obj)
            node_coords=obj.beam(1).node_coords';
            
            aQx=obj.beam(1).node_loadings_loc(1:6:end)';
            aQy=obj.beam(1).node_loadings_loc(2:6:end)';
            aQz=obj.beam(1).node_loadings_loc(3:6:end)';
            aMx=obj.beam(1).node_loadings_loc(4:6:end)';
            aMy=obj.beam(1).node_loadings_loc(5:6:end)';
            aMz=obj.beam(1).node_loadings_loc(6:6:end)';
            for i=2:length(obj.beam)
                node_coords=[node_coords obj.beam(i).node_coords'];
                aQx=[aQx obj.beam(i).node_loadings_loc(1:6:end)'];
                aQy=[aQy obj.beam(i).node_loadings_loc(2:6:end)'];
                aQz=[aQz obj.beam(i).node_loadings_loc(3:6:end)'];
                aMx=[aMx obj.beam(i).node_loadings_loc(4:6:end)'];
                aMy=[aMy obj.beam(i).node_loadings_loc(5:6:end)'];
                aMz=[aMz obj.beam(i).node_loadings_loc(6:6:end)'];
            end
            
            %scale bending moments
            
            node_coords=node_coords';
            max_x=max(node_coords(:,1));
            min_x=min(node_coords(:,1));
            
            max_y=max(node_coords(:,2));
            min_y=min(node_coords(:,2));
            
            max_z=max(node_coords(:,3));
            min_z=min(node_coords(:,3));
            
            delta_x=max_x-min_x;
            delta_z=max_z-min_z;
            delta_y=max_y-min_y;
            
            med_y=(max_y+min_y)/2;
            med_z=(max_z+min_z)/2;
            med_x=(max_x+min_x)/2;
            scale=0.35*sqrt((delta_x)^2+(delta_z)^2+(delta_y)^2);
            
            max_Mb=max([abs(aMx) abs(aMz) abs(aMy)]);
            %scale torsional moments
            max_Mt=max([abs(aMy)]);
            max_Q =max([abs(aQx) abs(aQy) abs(aQz)]);
            
            maxMx=max(abs(aMx))
            maxMt=max(abs(aMy))
            maxMz=max(abs(aMz))
            maxQx=max(abs(aQx))
            maxQy=max(abs(aQy))
            maxQz=max(abs(aQz))
            
            
            scale_bendingmoment=(2*max_Mb)/scale;
            scale_torsionmoment=(2*max_Mb)/scale;
            scale_force=(2*max_Q)/scale;
            for i=1:length(obj.beam)
                obj.beam(i).plot_internalforces(delta_x,delta_y,delta_z,med_x,med_y,med_z,scale_bendingmoment,scale_torsionmoment,scale_force,maxMx,maxMt,maxMz,maxQx,maxQy,maxQz);
            end
        end
        
        function  plot_critical_case_idx(obj)
            for i=1:length(obj.beam)
                if  isa(obj.beam(i),'class_wing')
                    obj.beam(i).plot_critical_case_idx;
                end
            end
        end
        
        function obj=plot_structure(obj)
            t_max=0;
            t_min=1000;
            for j=1:1:length(obj.beam)
                if  isa(obj.beam(j),'class_wing')
                    
                    for i=1:1:length(obj.beam(j).beamelement)
                        t_sk_up(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_up});
                        t_sk_lo(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_lo});
                        t_sp_fr(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_fr});
                        t_sp_re(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_re});
                    end
                    
                elseif isa(obj.beam(j),'class_fuselage')
                    
                    for i=1:1:length(obj.beam(j).beamelement)
                        t_sk_eq(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_eq});
                    end
                    if sum(t_sk_eq)>0
                        tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re t_sk_eq];
                    else
                        tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re];
                    end
                    t_max=max([tk t_max]);
                    t_min=min([tk t_min]);
                    
                end
            end
            t_max
            t_min
            for i=1:length(obj.beam)
                if  isa(obj.beam(i),'class_wing')
                obj.beam(i).plot_structure(t_max,t_min);
                elseif isa(obj.beam(i),'class_fuselage')
                obj.beam(i).plot_structure(t_max,t_min); 
                end
            end
        end
        
        
       function obj=write_structure_tecplot(obj,filename)
           t_max=0;
           t_min=1000;
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
           
           for j=1:1:length(obj.beam)
                if  isa(obj.beam(j),'class_wing')
                    for i=1:1:length(obj.beam(j).beamelement)
                        t_sk_up(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_up});
                        t_sk_lo(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_lo});
                        t_sp_fr(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_fr});
                        t_sp_re(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_re});
                    end
                elseif isa(obj.beam(j),'class_fuselage')
                    for i=1:1:length(obj.beam(j).beamelement)
                        t_sk_eq(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_eq});
                    end
                    if sum(t_sk_eq)>0
                        tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re t_sk_eq];
                    else
                        tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re];
                    end
                    t_max=max([tk t_max]);
                    t_min=min([tk t_min]);
              
                end
           end
           
           for i=1:1:length(obj.beam)
               if  isa(obj.beam(i),'class_wing')
                    obj.beam(i).write_structure_tecplot(fileID,i);
               elseif isa(obj.beam(i),'class_fuselage')
                   obj.beam(i).write_structure_tecplot(fileID,i);
               end
           end
       end
        
       function obj=write_fuel_tecplot(obj,filename)
           t_max=0;
           t_min=1000;
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
           
%            for j=1:1:length(obj.beam)
%                 if  isa(obj.beam(j),'class_wing')
%                     for i=1:1:length(obj.beam(j).beamelement)
%                         t_sk_up(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_up});
%                         t_sk_lo(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_lo});
%                         t_sp_fr(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_fr});
%                         t_sp_re(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sp_re});
%                     end
%                 elseif isa(obj.beam(j),'class_fuselage')
%                     for i=1:1:length(obj.beam(j).beamelement)
%                         t_sk_eq(i)=cell2mat({obj.beam(j).beamelement(i).crosssection.t_sk_eq});
%                     end
%                     if sum(t_sk_eq)>0
%                         tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re t_sk_eq];
%                     else
%                         tk=[t_sk_lo t_sk_up t_sp_fr t_sp_re];
%                     end
%                     t_max=max([tk t_max]);
%                     t_min=min([tk t_min]);
%               
%                 end
%            end
           
           for i=1:1:length(obj.beam)
               if  isa(obj.beam(i),'class_wing')
                    obj.beam(i).write_fuel_tecplot(fileID,i);
               elseif isa(obj.beam(i),'class_fuselage')
               end
           end
        end 
       
       
        function obj=f_init_material_properties(obj,structure)
            for i=1:length(obj.beam)
                obj.beam(i)=obj.beam(i).f_init_material_properties(structure(i));
            end
        end
        
        function obj=f_init_aeroloads(obj,wingno,results)
            for i=1:length(obj.beam)
                obj.beam(i)=obj.beam(i).f_init_aeroloads(i,results);
            end
        end
        
        function obj=write_tecplot_modes(obj,path,nmodes,exaturation)
            mkdir(path,'Modes2')
            for i=1:nmodes
                obj.nodal_deflections=obj.modeshapes(:,i)*exaturation;
                obj=obj.f_postprocess();
                obj.write_tecplot([path '/Modes2/modep' num2str(i)])
                obj.nodal_deflections=-obj.modeshapes(:,i)*exaturation;
                obj=obj.f_postprocess();
                obj.write_tecplot([path '/Modes2/modem' num2str(i)])
            end
        end
        
        
        function obj=write_tecplot(obj,filename)
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
           
                      Euler=[0 0 0]*pi/180;
           M_BI=[   cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
                    cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
                    cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
           
           
%            for j=1:1:length(obj.beam)
%                 fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
%                 for i=1:(length(obj.beam(j).beamelement)+1)
%                     xb=obj.beam(j).node_coords(i,1);
%                     yb=obj.beam(j).node_coords(i,2);
%                     zb=obj.beam(j).node_coords(i,3);
%                     nc=M_BI'*[xb;yb;zb];
%                     
%                     fprintf(fileID,'%f %f %f %f\n',nc(1),nc(2),nc(3),1);%obj.beam(j).node_loadings(1+6*(i-1))
%                 end
%                 for i=1:(length(obj.beam(j).beamelement))
%                     fprintf(fileID,'%i %i \n',i,i+1);
%                 end
%            end 
           

           for j=1:1:length(obj.beam)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
                for i=1:(length(obj.beam(j).beamelement)+1)
                    xb=obj.beam(j).node_coords(i,1)+obj.beam(j).nodal_deflections(1+6*(i-1));
                    yb=obj.beam(j).node_coords(i,2)+obj.beam(j).nodal_deflections(2+6*(i-1));
                    zb=obj.beam(j).node_coords(i,3)+obj.beam(j).nodal_deflections(3+6*(i-1));
                    nc=M_BI'*[xb;yb;zb];
                    
                    fprintf(fileID,'%f %f %f %f\n',nc(1),nc(2),nc(3),1);%obj.beam(j).node_loadings(1+6*(i-1))
                end
                for i=1:(length(obj.beam(j).beamelement))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
           end 

        end
        
        function obj=plot_deflected_grid(obj)
           figure
           hold on
           for j=1:1:length(obj.beam)
                for i=1:(length(obj.beam(j).beamelement)+1)
                    xb=obj.beam(j).node_coords(i,1)+obj.beam(j).nodal_deflections(1+6*(i-1));
                    yb=obj.beam(j).node_coords(i,2)+obj.beam(j).nodal_deflections(2+6*(i-1));
                    zb=obj.beam(j).node_coords(i,3)+obj.beam(j).nodal_deflections(3+6*(i-1));
                    plot3(xb,yb,zb,'bx')
                end
           end  
        end
        
        function obj=write_tecplot_mass(obj,filename)
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Mass"\n');
           
                      Euler=[0 0 0]*pi/180;
           M_BI=[   cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
                    cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
                    cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];

           for j=1:1:length(obj.beam)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
                for i=1:(length(obj.beam(j).beamelement)+1)
                    xb=obj.beam(j).node_coords(i,1);%+obj.beam(j).nodal_deflections(1+6*(i-1));
                    yb=obj.beam(j).node_coords(i,2);%+obj.beam(j).nodal_deflections(2+6*(i-1));
                    zb=obj.beam(j).node_coords(i,3);%+obj.beam(j).nodal_deflections(3+6*(i-1));
                    nc=M_BI'*[xb;yb;zb];
                    fprintf(fileID,'%f %f %f %f\n',nc(1),nc(2),nc(3),obj.beam(j).M_lumped(1+6*(i-1),1+6*(i-1)));
                end
                for i=1:(length(obj.beam(j).beamelement))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
           end  
        end
        
        function obj=write_tecplot_mode(obj,filename,mode_id,exaturation_factor)
            fileID = fopen([filename '.tp'],'w');
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
            
            obj.nodal_deflections=obj.modeshapes(:,mode_id)*exaturation_factor;
            obj=obj.f_postprocess();

            for j=1:1:length(obj.beam)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
                for i=1:(length(obj.beam(j).beamelement)+1)
                    j
                    i
                    fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1)+obj.beam(j).nodal_deflections(1+6*(i-1)),obj.beam(j).node_coords(i,2)+obj.beam(j).nodal_deflections(2+6*(i-1)),obj.beam(j).node_coords(i,3)+obj.beam(j).nodal_deflections(3+6*(i-1)),1);%obj.beam(j).node_loadings(1+6*(i-1))
                end
                for i=1:(length(obj.beam(j).beamelement))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
            end
        end
        
        function obj=write_tecplot_acc(obj,filename)
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
           
           for j=1:1:length(obj.beam)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
                for i=1:(length(obj.beam(j).beamelement)+1)
                    fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),1);%obj.beam(j).node_loadings(1+6*(i-1))
                end
                for i=1:(length(obj.beam(j).beamelement))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
           end  
           
           for j=1:1:length(obj.beam)
               for i=1:(length(obj.beam(j).beamelement))
                   fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4*length(obj.beam(j).beamelement),length(obj.beam(j).beamelement));
                   for i =1:length(obj.beam(j).beamelement)
                       acc_z=obj.beam(j).beamelement(i).az; 
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+acc_z/10,1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+acc_z/10,1);
                   end
                   k=1;
                   for i =1:length(obj.beam(j).beamelement)
                       fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
                       k=k+4;
                   end
               end
           end
           
           for j=1:1:length(obj.beam)
               for i=1:(length(obj.beam(j).beamelement))
                   fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4*length(obj.beam(j).beamelement),length(obj.beam(j).beamelement));
                   for i =1:length(obj.beam(j).beamelement)
                       acc_y=obj.beam(j).beamelement(i).ay;
                       
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),2);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),2);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2)+acc_y/10,obj.beam(j).node_coords(i+1,3),2);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2)+acc_y/10,obj.beam(j).node_coords(i,3),2);
                   end
                   k=1;
                   for i=1:length(obj.beam(j).beamelement)
                       fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
                       k=k+4;
                   end
               end
           end
           
           for j=1:1:length(obj.beam)
               for i=1:(length(obj.beam(j).beamelement))
                   fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4*length(obj.beam(j).beamelement),length(obj.beam(j).beamelement));
                   for i=1:length(obj.beam(j).beamelement)
                       acc_x=obj.beam(j).beamelement(i).ax;
                       
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),3);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),3);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1)+acc_x/10,obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),3);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1)+acc_x/10,obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),3);
                   end
                   k=1;
                   for i =1:length(obj.beam(j).beamelement)
                       fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
                       k=k+4;
                   end
               end
           end
           
        end
        
        
        function obj=write_tecplot_extloads(obj,filename)
           ld=5000;
           fileID = fopen([filename '.tp'],'w');
           fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
           fprintf(fileID,'VARIABLES = "X", "Y", "Z","Equivalent Thickness"\n');
           for j=1:1:length(obj.beam)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',length(obj.beam(j).beamelement)+1,length(obj.beam(j).beamelement));
                for i=1:(length(obj.beam(j).beamelement)+1)
                    fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),1);%obj.beam(j).node_loadings(1+6*(i-1))
                end
                for i=1:(length(obj.beam(j).beamelement))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                if ~isempty(obj.beam(j).engine)
                    for ii=1:length(obj.beam(j).engine)
                        dist=zeros(obj.beam(j).nel,1);
                        for i=1:obj.beam(j).nel
                            dist(i)=sqrt(sum((obj.beam(j).node_coords(i,2)-obj.beam(j).engine(ii).cg_pos(2)').^2));
                        end
                        [Y,I] = min(dist);
                        fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',2,1);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).engine(ii).cg_pos(1),obj.beam(j).engine(ii).cg_pos(2) ,obj.beam(j).engine(ii).cg_pos(3),1);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(I,1),obj.beam(j).node_coords(I,2),obj.beam(j).node_coords(I,3),1);
                        fprintf(fileID,'1 2 \n');
                        fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4,1);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).engine(ii).cg_pos(1),obj.beam(j).engine(ii).cg_pos(2),obj.beam(j).engine(ii).cg_pos(3),6);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).engine(ii).cg_pos(1),obj.beam(j).engine(ii).cg_pos(2),obj.beam(j).engine(ii).cg_pos(3)-obj.beam(j).engine(ii).m*obj.beam(j).load_factor*9.81/ld/2,6);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).engine(ii).cg_pos(1)+2,obj.beam(j).engine(ii).cg_pos(2),obj.beam(j).engine(ii).cg_pos(3)-obj.beam(j).engine(ii).m*obj.beam(j).load_factor*9.81/ld/2,6);
                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).engine(ii).cg_pos(1)+2,obj.beam(j).engine(ii).cg_pos(2),obj.beam(j).engine(ii).cg_pos(3),6);
                        fprintf(fileID,'1 2 3 4\n');
                    end
                end 
           end  

           
           for j=1:1:length(obj.beam)     
               if  isa(obj.beam(j),'class_wing')
                   nq=0;
                   fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',3*4*length(obj.beam(j).beamelement),3*length(obj.beam(j).beamelement));
                   for i=1:length(obj.beam(j).beamelement)
                       %aerodynamic load
                       if (i==1)
                           q=obj.beam(j).beamelement(i).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i).qz];
                           qp1=obj.beam(j).beamelement(i+1).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i+1).qz];
                           q_z1=q(3);
                           q_z2=(q(3)+qp1(3))/2;
                       elseif  (i==length(obj.beam(j).beamelement))
                           qm1=obj.beam(j).beamelement(i-1).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i-1).qz];
                           q=obj.beam(j).beamelement(i).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i).qz];
                           q_z2=q(3);
                           q_z1=(q(3)+qm1(3))/2;
                       else
                           qm1=obj.beam(j).beamelement(i-1).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i-1).qz];
                           q=obj.beam(j).beamelement(i).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i).qz];
                           qp1=obj.beam(j).beamelement(i+1).T(1:3,1:3)*[0;0;obj.beam(j).beamelement(i+1).qz];
                           q_z1=(q(3)+qm1(3))/2;
                           q_z2=(q(3)+qp1(3))/2;
                       end

                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+q_z2/ld,1);
                       fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+q_z1/ld,1);
                       nq=nq+1;
                       % inertial mass load
                       if ~isequal(sign(q_z1),sign(obj.beam(j).beamelement(i).az))
                           inertial_load=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).m;
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,2);
                           nq=nq+1;
                           % inertial fuel load
                           if obj.beam(j).beamelement(i).is_fueled==1
                               fuel_load=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).el_m_fuel/obj.beam(j).beamelement(i).le;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld+fuel_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld+fuel_load/ld,3);
                               nq=nq+1;
                           else
                               fuel_load=0;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld+fuel_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld+fuel_load/ld,3);
                               nq=nq+1;
                           end
                       else
                           inertial_load=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).m;
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,2);
                           nq=nq+1;
                           % inertial fuel load
                           if obj.beam(j).beamelement(i).is_fueled==1
                               fuel_load=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).el_m_fuel/obj.beam(j).beamelement(i).le;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld+fuel_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld+fuel_load/ld,3);
                               nq=nq+1;
                           else
                               fuel_load=0;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_load/ld+fuel_load/ld,3);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_load/ld+fuel_load/ld,3);
                               nq=nq+1;
                           end
                       end
                   end
                   

                   k=1;
                   for i =1:nq
                       fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
                       k=k+4;
                   end
                   
               elseif isa(obj.beam(j),'class_fuselage')
                   
                    nq=0;
                   fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',3*4*length(obj.beam(j).beamelement),3*length(obj.beam(j).beamelement));
                   for i=1:length(obj.beam(j).beamelement)
                       %aerodynamic load            
                           inertial_primary=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).el_m_p/obj.beam(j).beamelement(i).le;
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_primary/ld,2);
                           fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_primary/ld,2);
                           nq=nq+1;
                           % inertial fuel load
                           
                               payload=obj.beam(j).beamelement(i).az*obj.beam(j).beamelement(i).el_m_s/obj.beam(j).beamelement(i).le;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_primary/ld,4);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_primary/ld,4);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_primary/ld+payload*0.2/ld,4);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_primary/ld+payload*0.2/ld,4);
                               nq=nq+1;
                               fuel_load=0;
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_primary/ld+payload*0.2/ld,5);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_primary/ld+payload*0.2/ld,5);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3)+inertial_primary/ld+payload*0.2/ld+payload*0.8/ld,5);
                               fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3)+inertial_primary/ld+payload*0.2/ld+payload*0.8/ld,5);
                               nq=nq+1;
                   end
                   

                   k=1;
                   for i =1:nq
                       fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
                       k=k+4;
                   end
                                 
               end
               
            
            %   end
           end
           
%            for j=1:1:length(obj.beam)
%                for i=1:(length(obj.beam(j).beamelement))
%                    fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4*length(obj.beam(j).beamelement),length(obj.beam(j).beamelement));
%                    for i =1:length(obj.beam(j).beamelement)
%                        acc_y=obj.beam(j).beamelement(i).ay;
%                        
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),2);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),2);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2)+acc_y/10,obj.beam(j).node_coords(i+1,3),2);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2)+acc_y/10,obj.beam(j).node_coords(i,3),2);
%                    end
%                    k=1;
%                    for i=1:length(obj.beam(j).beamelement)
%                        fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
%                        k=k+4;
%                    end
%                end
%            end
%            
%            for j=1:1:length(obj.beam)
%                for i=1:(length(obj.beam(j).beamelement))
%                    fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',4*length(obj.beam(j).beamelement),length(obj.beam(j).beamelement));
%                    for i=1:length(obj.beam(j).beamelement)
%                        acc_x=obj.beam(j).beamelement(i).ax;
%                        
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1),obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),3);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1),obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),3);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i+1,1)+acc_x/10,obj.beam(j).node_coords(i+1,2),obj.beam(j).node_coords(i+1,3),3);
%                        fprintf(fileID,'%f %f %f %f\n',obj.beam(j).node_coords(i,1)+acc_x/10,obj.beam(j).node_coords(i,2),obj.beam(j).node_coords(i,3),3);
%                    end
%                    k=1;
%                    for i =1:length(obj.beam(j).beamelement)
%                        fprintf(fileID,'%i %i %i %i \n',k,k+1,k+2,k+3);
%                        k=k+4;
%                    end
%                end
%            end
           
        end
        
        function obj=plot_externalforces(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces)
            
            qx=0;
            qy=0;
            qz=0;
            q_max=0;
            for j=1:1:length(obj.beam)
                for i=1:1:length(obj.beam(j).beamelement)
                    qx(i)=cell2mat({obj.beam(j).beamelement(i).qx});
                    qy(i)=cell2mat({obj.beam(j).beamelement(i).qy});
                    qz(i)=cell2mat({obj.beam(j).beamelement(i).qz});
                end
                qt=max([abs(qx),abs(qy),abs(qz)]);
                q_max=max([q_max qt])
            end
            
            halfspan=sum(cell2mat({obj.beam(1).beamelement(:).le}));
            
            %max_q=max([abs(qx),abs(qy),abs(qz)]);
            
            scale=2*q_max/halfspan;
            
            for i=1:length(obj.beam)
                if isa(obj.beam(i),'class_wing')
                    obj.beam(i).plot_externalforces(add_eigenmass,add_fuelmass,add_engineforces,add_gearforces,scale);
                end
            end
        end
        
        function obj=f_structural_layout(obj,overwrite)            
            for i=1:1:length(obj.beam)
                if obj.beam(i).isExternalFEM==0
                    obj.beam(i)=obj.beam(i).f_structural_layout(overwrite);
                end
            end
        end
        
        function obj=f_load_based_self_design(obj,weights,varargin)
            
            if nargin==3
                overwrite=varargin{1};
                iterate_mass_flag=1;
            elseif nargin==4
                overwrite=varargin{1};
                iterate_mass_flag=varargin{2};
            else
                overwrite=1;
                iterate_mass_flag=1;
            end
            
  
            mem=obj.settings.nonlinear;
            obj.settings.nonlinear=0;
            % solve for initial internal forces
            obj=obj.f_solve();
            obj.settings.nonlinear=mem;
            % first perform initial layout
            obj=obj.f_structural_layout(overwrite);
            % calculate mass
            obj=obj.f_calc_mass(weights);
            
            prev_val=obj.m_total;
            err=100;
            its=1;
            %% Structural Layout Iteration
            if iterate_mass_flag==1
                while(err>1)
                    %% perform structural layout (determine wing box size from loading)
                    % resolve system with new mass
                    obj=obj.f_solve();
                    
                    %figure
                    %obj.plot_internalforces;
                    obj=obj.f_structural_layout(overwrite);
                    % calculate mass
                    obj=obj.f_calc_mass(weights);
                    
                    err=100*abs((prev_val-obj.m_total)/abs(obj.m_total));
                    
                    disp(['            self design iteration: ' num2str(its) '      error: ' num2str(err) '% mass' num2str(obj.m_total)]);
                    prev_val=obj.m_total;
                    its=its+1;
                end
            end
        end
        
        function obj=f_load_based_self_design_unrestrained(obj,weights,varargin)
            
            if nargin==3
                overwrite=varargin{1};
                iterate_mass_flag=1;
            elseif nargin==4
                overwrite=varargin{1};
                iterate_mass_flag=varargin{2};
            else
                overwrite=1;
                iterate_mass_flag=1;
            end
            
            mem=obj.settings.nonlinear;
            obj.settings.nonlinear=0;
            % solve for initial internal forces
            obj=obj.f_solve_unrestrained();
            obj.settings.nonlinear=mem;
            % first perform initial layout
            obj=obj.f_structural_layout(overwrite);
            % calculate mass
            obj=obj.f_calc_mass(weights);
            
            prev_val=obj.m_total;
            err=100;
            its=1;
            %% Structural Layout Iteration
            if iterate_mass_flag==1
                while(err>1)
                    %% perform structural layout (determine wing box size from loading)
                    % resolve system with new mass
                    obj=obj.f_solve_unrestrained();
                    
                    %figure
                    %obj.plot_internalforces;
                    obj=obj.f_structural_layout(overwrite);
                    % calculate mass
                    obj=obj.f_calc_mass(weights);
                    
                    err=100*abs((prev_val-obj.m_total)/abs(obj.m_total));
                    
                    disp(['            self design iteration: ' num2str(its) '      error: ' num2str(err) '% mass' num2str(obj.m_total)]);
                    prev_val=obj.m_total;
                    its=its+1;
                end
            end
        end
        
        function obj=f_calc_mass(obj,weights)
            obj.m_total=0;
            k=0;
            for i=1:1:length(obj.beam)
                if  isa(obj.beam(i),'class_wing')
                    if obj.beam(i).isExternalFEM==0
                        wingweights.WingSystemsEstimate=weights.WingSystemsEstimate(i);
                        wingweights.WingSkinEstimate=weights.WingSkinEstimate(i);
                        wingweights.WingInitialGuess=weights.WingInitialGuess(i);
                        obj.beam(i)=obj.beam(i).f_calc_mass(wingweights);

                        obj.m_total=obj.m_total+obj.beam(i).m_total;
                    elseif obj.beam(i).isExternalFEM==1
                        obj.m_total = obj.m_total + obj.beam(i).f_compute_totalMass;
                    end
                    k=k+1;
                elseif isa(obj.beam(i),'class_fuselage') && ~obj.beam(i).is_rigid
                    if obj.beam(i).isExternalFEM==0
                        fuselageweights.FuselageSystemsEstimate=weights.FuselageSystemsEstimate(i-k);
                        fuselageweights.FuselageNonStructuralEstimate=weights.FuselageNonStructuralEstimate(i-k);
                        fuselageweights.NumberPAX=weights.FuselagePAX(i-k);
                        obj.beam(i)=obj.beam(i).f_calc_mass(fuselageweights);
                        obj.m_total=obj.m_total+obj.beam(i).m_total;
                    elseif obj.beam(i).isExternalFEM==1
                        obj.m_total = obj.m_total + obj.beam(i).f_compute_totalMass;
                    end   
                elseif isa(obj.beam(i),'class_pylon') && ~obj.beam(i).is_rigid
                    % fuselageweights.FuselageSystemsEstimate=weights.FuselageSystemsEstimate(i-k);
                    % fuselageweights.FuselageNonStructuralEstimate=weights.FuselageNonStructuralEstimate(i-k);
                    obj.beam(i)=obj.beam(i).f_calc_mass(fuselageweights);
                    obj.m_total=obj.m_total+obj.beam(i).m_total;
                end
                if obj.beam(i).isExternalFEM==0
                    obj.m_total=obj.m_total+sum(obj.beam(i).nodal_masses(:,1));
                end
            end
        end
        
        function [T,db]=assemble_force_transformation_matrix(obj,aircraft)
            
            for i=1:length(aircraft.wings)
                len=length(obj.beam(i).Q);
                if ~isempty(obj.beam(i).sort_vec)
                    obj.beam(i).T_sort=[aircraft.wings(i).T obj.beam(i).sort_vec(:,1)];
                    obj.beam(i).T_sort=sortrows(obj.beam(i).T_sort,len+1);
                    obj.beam(i).T_sort=aircraft.wings(i).T_sort(1:len,1:len)';
                    obj.beam(i).T_sort=[obj.beam(i).T_sort obj.beam(i).sort_vec(:,1)];
                    obj.beam(i).T_sort=sortrows(obj.beam(i).T_sort,len+1);
                    obj.beam(i).T_sort=obj.beam(i).T_sort(1:len,1:len)';
                    k=size(obj.beam(i).Kff,2)-len+1;
                    obj.beam(i).T_sort=aircraft.wings(i).T_sort(k:end,k:end);
                    
                    obj.beam(i).db_sort=[aircraft.wings(i).db obj.beam(i).sort_vec(:,1)];
                    obj.beam(i).db_sort=sortrows(obj.beam(i).db_sort,len+1);
                    obj.beam(i).db_sort=aircraft.wings(i).db_sort(1:len,1:len)';
                    obj.beam(i).db_sort=[obj.beam(i).db_sort obj.beam(i).sort_vec(:,1)];
                    obj.beam(i).db_sort=sortrows(obj.beam(i).db_sort,len+1);
                    obj.beam(i).db_sort=obj.beam(i).db_sort(1:len,1:len)';
                    k=size(obj.beam(i).Kff,2)-len+1;
                    obj.beam(i).db_sort=aircraft.wings(i).db_sort(k:end,k:end);
                else
                    obj.beam(i).T_sort=aircraft.wings(i).T(1:end,1:end);
                    obj.beam(i).db_sort=aircraft.wings(i).db(1:end,1:end);
                end
            end
            
            if length(aircraft.wings)<length(obj.beam)
                for i=(length(aircraft.wings)+1):1:length(obj.beam)
                    obj.beam(i).T_sort=zeros(size(obj.beam(i).Kff,2),size(aircraft.wings(1).T,2));
                    obj.beam(i).db_sort=zeros(size(obj.beam(i).Kff,2),size(aircraft.wings(1).db,2));
                end
            end
            
            % calculate overall system DOFs
            for i=1:1:length(obj.beam)
                sys_ndof(i)=length(obj.beam(i).Kff);
            end
            % preallocate system global aerostructural transformation
            % matrix
            obj.T_glob=zeros(sum(sys_ndof),size(aircraft.wings(1).T,2));
            obj.db_glob=zeros(sum(sys_ndof),size(aircraft.wings(1).db,2));
            
            obj.sort_vec=zeros(sum(sys_ndof),2);
            obj.dof_node_beam=zeros(sum(sys_ndof),3);
            
            sys_pos=cumsum(sys_ndof);
            sys_pos=[0 sys_pos(1:end-1)];
            
            for i=1:1:length(obj.beam)
                obj.T_glob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),:)=obj.T_glob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),:)+obj.beam(i).T_sort
                obj.db_glob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),:)=obj.db_glob(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),:)+obj.beam(i).db_sort;
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),3)=i;
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),2)=obj.beam(i).Kff_node_dof_info(:,1);
                obj.dof_node_beam(sys_pos(i)+1:sys_pos(i)+sys_ndof(i),1)=obj.beam(i).Kff_node_dof_info(:,2);
            end
            
            el_ndof=6;
            
            ctr=1;
            for i=1:1:length(obj.coupling_condition)
                for j=1:1:el_ndof
                    if(obj.coupling_condition(i).dof(j)==1)
                        for k=1:1:length(obj.coupling_condition(i).beam_idx)
                            idx=find(obj.beam(obj.coupling_condition(i).beam_idx(k)).Kff_node_dof_info(:,1)==obj.coupling_condition(i).node_idx(k),1,'first');
                            sort_idx=sys_pos(obj.coupling_condition(i).beam_idx(k))+idx+j-1;
                            obj.sort_vec(sort_idx,1)=ctr;
                            obj.sort_vec(sort_idx,2)=sort_idx;
                        end
                        ctr=ctr+1;
                    end
                end
            end
            
            j=ctr;
            for i=1:sum(sys_ndof)
                if obj.sort_vec(i,1)==0
                    obj.sort_vec(i,1)=j;
                    obj.sort_vec(i,2)=i;
                    j=j+1;
                end
            end
            
            len=sum(sys_ndof);
            col_ndof=max(obj.sort_vec(:,1));
            tot_size=sum(sys_ndof);
            % T Matrix
            Tsort=[obj.T_glob obj.sort_vec(:,1)];
            Tsort=sortrows(Tsort,2);
            T=zeros(col_ndof,size(aircraft.wings(1).T,2));
            for i=1:1:max(obj.sort_vec(:,1))
                idx=find(Tsort(:,end)==i);
                if length(idx)>1
                    T(i,:)=sum(Tsort(idx,1:end-1));
                else
                    T(i,:)=Tsort(idx,1:end-1);
                end
            end
            % db Matrix
            
            dbsort=[obj.db_glob obj.sort_vec(:,1)];
            dbsort=sortrows(dbsort,2);
            db=zeros(col_ndof,size(aircraft.wings(1).db,2));
            for i=1:1:max(obj.sort_vec(:,1))
                idx=find(dbsort(:,end)==i);
                if length(idx)>1
                    db(i,:)=sum(dbsort(idx,1:end-1));
                else
                    db(i,:)=dbsort(idx,1:end-1);
                end
            end
        end

        function obj=f_solve_unrestrained(obj)
            
            eval_nonlin=obj.settings.nonlinear;
            obj.err=100;
            
%             if eval_nonlin==1
%                 for i=1:1:length(obj.beam)
%                     obj.beam(i)=obj.beam(i).f_reset();
%                 end
%                 if ~isempty(obj.nodal_deflections)
%                     obj.nodal_deflections=obj.nodal_deflections*0;
%                 end
%             end
            shape=obj.modefrequencies*0;
            loadstep=1;
            if eval_nonlin==1
                loadstep=1;
            end
            
            while loadstep<=1
                obj.err=100;

                if eval_nonlin==1
                    disp(['                 nonlinear loadstep: ' num2str(loadstep)]);
                end
                itctr=1;
                while obj.err>0.1
                    %obj=obj.f_assemble_free(loadstep,eval_nonlin);
                    obj=obj.f_solve_free_modes(1,0);
                    shape=obj.modefrequencies*0;
                    %Mdiag=obj.modeshapes'*obj.Mff*obj.modeshapes;
                    Kdiag=obj.modeshapes'*obj.Kff*obj.modeshapes;
                    Qdiag=obj.modeshapes'*obj.Ftest;
                    x_now=Kdiag(7:end,7:end)^-1*Qdiag(7:end);
                    for ss=1:length(obj.modefrequencies)-6
                        shape=shape+obj.modeshapes(:,6+ss)*x_now(ss);
                    end
                    obj.nodal_deflections=shape;
                    obj.modal_deflections=x_now;
%                     if eval_nonlin==1
%                         R=obj.Ptest-obj.Ftest;
%                         err=sum(abs(R))/sum(abs(obj.Ftest));%sum(abs(R(2:6:end)))/sum(abs(obj.Ftest(2:6:end)))+sum(abs(R(3:6:end)))/sum(abs(obj.Ftest(3:6:end)))+sum(abs(R(4:6:end)))/sum(abs(obj.Ftest(4:6:end)))+sum(abs(R(5:6:end)))/sum(abs(obj.Ftest(5:6:end)));
%                         if isnan(err)
%                             err=10E7;
%                         end
%                         if eval_nonlin==1
%                             disp(['                    nonlinear iteration: ' num2str(itctr) '    error: ' num2str(err)]);
%                             itctr=itctr+1;
%                             if itctr>1000
%                                 break;
%                             end
%                         end
%                         obj.err=err;
%                         if isempty(obj.nodal_deflections)
%                             obj.nodal_deflections=linsolve(obj.Kff,-R);
%                         else
%                             obj.nodal_deflections=obj.nodal_deflections+linsolve(obj.Kff,-R);
%                         end
%                     else
%                         obj.nodal_deflections=linsolve(obj.Kff,obj.Ftest);
%                     end
                    obj=obj.f_postprocess();
                end
                loadstep=loadstep+0.2;
            end
        end
                
            
        function obj=f_solve(obj)
            
            eval_nonlin=obj.settings.nonlinear;
            obj.err=100;
            
            if eval_nonlin==1
                for i=1:1:length(obj.beam)
                    obj.beam(i)=obj.beam(i).f_reset();
                end
                if ~isempty(obj.nodal_deflections)
                    obj.nodal_deflections=obj.nodal_deflections*0;
                end
            end
            
            loadstep=1;
            if eval_nonlin==1
                loadstep=1;
            end
            
            while loadstep<=1
                obj.err=100;
                
                if eval_nonlin==1
                    disp(['                 nonlinear loadstep: ' num2str(loadstep)]);
                end
                itctr=1;
                while obj.err>0.1
                    obj=obj.f_assemble(loadstep,eval_nonlin);
                    
                    if eval_nonlin==1
                        R=obj.Ptest-obj.Ftest;
                        err=sum(abs(R))/sum(abs(obj.Ftest));%sum(abs(R(2:6:end)))/sum(abs(obj.Ftest(2:6:end)))+sum(abs(R(3:6:end)))/sum(abs(obj.Ftest(3:6:end)))+sum(abs(R(4:6:end)))/sum(abs(obj.Ftest(4:6:end)))+sum(abs(R(5:6:end)))/sum(abs(obj.Ftest(5:6:end)));
                        if isnan(err)
                            err=10E7;
                        end
                        if eval_nonlin==1
                            disp(['                    nonlinear iteration: ' num2str(itctr) '    error: ' num2str(err)]);
                            itctr=itctr+1;
                            if itctr>1000
                                break;
                            end
                        end
                        obj.err=err;
                        if isempty(obj.nodal_deflections)
                            obj.nodal_deflections=linsolve(obj.Kff,-R);
                        else
                            obj.nodal_deflections=obj.nodal_deflections+linsolve(obj.Kff,-R);
                        end
                    else
                        obj.nodal_deflections=linsolve(obj.Kff,obj.Ftest);
                    end
                    
                    obj=obj.f_postprocess();
                end
                loadstep=loadstep+0.2;
            end
        end
        
        function obj=f_postprocess(obj)
            eval_nonlin=obj.settings.nonlinear;
            for i=1:obj.n
                xbeam(i).nodal_deflections(1:obj.beam(i).ndof)=zeros(obj.beam(i).ndof,1);
            end
            tot_size=obj.n_dof_sys;
            for i=1:1:tot_size
                dof=obj.dof_node_beam(i,1);
                node=obj.dof_node_beam(i,2);
                beam_idx=obj.dof_node_beam(i,3);
                idx=(node-1)*6+dof;
                xbeam(beam_idx).nodal_deflections(idx)=obj.nodal_deflections(obj.sort_vec(i,1));
            end
            
            for i =1:obj.n
                obj.beam(i).nodal_deflections(:)=xbeam(i).nodal_deflections(:);
            end
            
            % only for now
            % obj.nodal_deflections=xbeam(1).nodal_deflections(:);
            
            for i=1:1:length(obj.beam)
                %if obj.beam(i).isExternalFEM==0
                    if eval_nonlin==1
                        obj.beam(i)=obj.beam(i).f_postprocess_nonlin();
                    else
                        obj.beam(i)=obj.beam(i).f_postprocess();
                        obj.err=0;
                    end
                %end
            end
        end
        
        
        
        function obj=f_solve_free_modes(obj,varargin)
            
            
            if(nargin==4)
                p_ref=varargin{1};
                nonlin=varargin{2};
                Kprv=varargin{3};
            else
                nonlin=0;
                Kprv=obj.Kff*0;
                p_ref=[0 0 0];
            end
            
            obj=obj.f_assemble_free(1,0);
            K=obj.Kff;
            M=obj.Mff;
            if(sum(sum(M)))==0 || isnan(sum(sum(M)))
                [obj.modeshapes, omega2]= eig(K*1E-8,K);
            else
                [obj.modeshapes, omega2]= eig(M,K);
            end
            
            obj.modefrequencies = (1./sqrt(diag(omega2)))/(2*pi); % Since inverse iteration is used in determining matrix A
            [obj.modefrequencies,index]=sort(obj.modefrequencies);
            obj.modeshapes(:,:)=obj.modeshapes(:,index);
            
            %replace the rigid body modes by unit deflections about body
            %axis
%             obj.modeshapes(:,1:6)=0.0;
%             obj.modeshapes(1:6:end,1)=0.1;
%             obj.modeshapes(2:6:end,2)=0.1;
%             obj.modeshapes(3:6:end,3)=0.1;
%             n_dof=length(obj.nodal_deflections);
%             I_Hat=[];
%             for i=1:n_dof/6
%                 I_Hat=[I_Hat; eye(3,3);zeros(3,3)];
%             end
%             Euler=[0.01 0 0];
%             M_BI=[  cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
%                 cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
%                 cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
%             k=1;
%             for j=1:1:size((obj.node_coords),1)
%                 s_inertial(k:k+5)=[M_BI'*(obj.node_coords(j,1:3)-p_ref)' [0.01 0 0]'];
%                 s_inertial(k:k+2)=s_inertial(k:k+2)-(obj.node_coords(j,1:3)-p_ref);
%                 k=k+6;
%             end
%             obj.modeshapes(:,4)= s_inertial;
%             Euler=[0 0.01 0];
%             M_BI=[  cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
%                 cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
%                 cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
%             k=1;
%             for j=1:1:size((obj.node_coords),1)
%                 s_inertial(k:k+5)=[M_BI'*(obj.node_coords(j,1:3)-p_ref)' [0 0.01 0]'];
%                 s_inertial(k:k+2)=s_inertial(k:k+2)-(obj.node_coords(j,1:3)-p_ref);
%                 k=k+6;
%             end
%             obj.modeshapes(:,5)= s_inertial;
%             Euler=[0 0 0.01];
%             M_BI=[  cos(Euler(3))*cos(Euler(2))                                             sin(Euler(3))*cos(Euler(2))                                                     -sin(Euler(2));
%                 cos(Euler(3))*sin(Euler(2))*sin(Euler(1))-sin(Euler(3))*cos(Euler(1))       sin(Euler(3))*sin(Euler(2))*sin(Euler(1))+cos(Euler(3))*cos(Euler(1))           cos(Euler(2))*sin(Euler(1));
%                 cos(Euler(3))*sin(Euler(2))*cos(Euler(1))+sin(Euler(3))*sin(Euler(1))       sin(Euler(3))*sin(Euler(2))*cos(Euler(1))-cos(Euler(3))*sin(Euler(1))           cos(Euler(2))*cos(Euler(1));];
%             k=1;
%             for j=1:1:size((obj.node_coords),1)
%                 s_inertial(k:k+5)=[M_BI'*(obj.node_coords(j,1:3)-p_ref)' [0 0 -0.01]'];
%                 s_inertial(k:k+2)=s_inertial(k:k+2)-(obj.node_coords(j,1:3)-p_ref);
%                 k=k+6;
%             end
%             obj.modeshapes(:,6)= s_inertial;
            
            %A=linsolve(K,M); % inverse iteration method (http://en.wikipedia.org/wiki/Modal_analysis_using_FEM )
            try
            [obj.modeshapes_lumped, omega2_lumped]= eig(obj.Mff_lumped,K);
                        obj.modefrequencies_lumped = (1./sqrt(diag(omega2_lumped)))/(2*pi);
                                    [obj.modefrequencies_lumped,index_lumped]=sort(obj.modefrequencies_lumped);
            obj.modeshapes_lumped(:,:)=obj.modeshapes_lumped(:,index_lumped);
                        Msx_lumped=real(obj.modeshapes_lumped)'*obj.Mff_lumped*real(obj.modeshapes_lumped);
            Ksx_lumped=real(obj.modeshapes_lumped)'*obj.Kff*real(obj.modeshapes_lumped);
            norm_factor_lumped=diag(Msx_lumped);
            for i=1:length(norm_factor_lumped)
                obj.modeshapes_lumped(:,i)=real(obj.modeshapes_lumped(:,i)/sqrt(norm_factor_lumped(i)));
            end
            end
            % Script to extract the relevant information from the eigen modes matrix
            %             j=1;
            %             for i = 1:size(obj.Mff,2)
            %                 [ EMA(j), idx(j)] = max(abs(obj.modeshapes(:,i))); % returns the maximum of each column vector from the V matrix
            %                 r(j) = mod(idx(j),6); % output 0 means the 6th DOF
            %                 j=j+1;
            %             end


            % %omega = sqrt(diag(omega2));
            


            
            
            %% norm modes for reduced mass=1
            Msx=real(obj.modeshapes)'*obj.Mff*real(obj.modeshapes);
            Ksx=real(obj.modeshapes)'*obj.Kff*real(obj.modeshapes);
            norm_factor=diag(Msx);
            for i=1:length(norm_factor)
                obj.modeshapes(:,i)=real(obj.modeshapes(:,i)/sqrt(norm_factor(i)));
            end
            

        end
        
        function obj=f_solve_modes(obj)
            K=obj.Kff;
            M=obj.Mff;
            %A=K\M; % inverse iteration method (http://en.wikipedia.org/wiki/Modal_analysis_using_FEM )
            % check inverse iteration (where is the iteration?) matlab eig
            % delivers better results
            [obj.modeshapes, omega2]= eig(M,K);
            [obj.modeshapes_lumped, omega2_lumped]= eig(obj.Mff_lumped,K);
            
            % Script to extract the relevant information from the eigen modes matrix
            
            j=1;
            for i = 1:size(obj.Mff,2)
                [ EMA(j), idx(j)] = max(abs(obj.modeshapes(:,i))); % returns the maximum of each column vector from the V matrix
                r(j) = mod(idx(j),6); % output 0 means the 6th DOF
                j=j+1;
            end
            
            obj.modefrequencies = (1./sqrt(diag(omega2)))/(2*pi); % Since inverse iteration is used in determining matrix A
            obj.modefrequencies_lumped = (1./sqrt(diag(omega2_lumped)))/(2*pi);
            % %omega = sqrt(diag(omega2));
            
            [obj.modefrequencies,index]=sort(obj.modefrequencies);
            obj.modeshapes(:,:)=obj.modeshapes(:,index);
            
            [obj.modefrequencies_lumped,index_lumped]=sort(obj.modefrequencies_lumped);
            obj.modeshapes_lumped(:,:)=obj.modeshapes_lumped(:,index_lumped);
            
            
            %% norm modes for reduced mass=1
            Msx=real(obj.modeshapes)'*obj.Mff*real(obj.modeshapes);
            Ksx=real(obj.modeshapes)'*obj.Kff*real(obj.modeshapes);
            norm_factor=diag(Msx);
            for i=1:length(norm_factor)
                obj.modeshapes(:,i)=real(obj.modeshapes(:,i)/sqrt(norm_factor(i)));
            end
            
            Msx_lumped=real(obj.modeshapes_lumped)'*obj.Mff_lumped*real(obj.modeshapes_lumped);
            Ksx_lumped=real(obj.modeshapes)'*obj.Kff*real(obj.modeshapes);
            norm_factor_lumped=diag(Msx_lumped);
            for i=1:length(norm_factor_lumped)
                obj.modeshapes_lumped(:,i)=real(obj.modeshapes_lumped(:,i)/sqrt(norm_factor_lumped(i)));
            end
            
            %             for i=1:length(obj.modeshapes(1,:))
            %                 obj.modeshapes(:,i)=obj.modeshapes(:,i)/(max(abs(obj.modeshapes(:,i))));
            %             end
            %
            %             kk=1;
            %             tol=0.1;
            %             offset=0;
            %             while kk<length(obj.modefrequencies)-1
            %                if(abs(obj.modefrequencies(kk)-obj.modefrequencies(kk+1))<tol)
            %                    obj.modefrequencies(kk)=(obj.modefrequencies(kk)+obj.modefrequencies(kk+1))/2;
            %                    obj.modefrequencies(kk+1:end-1)=obj.modefrequencies(kk+2:end);
            %                    obj.modeshapes(:,kk)=(obj.modeshapes(:,kk)+obj.modeshapes(:,kk+1));
            %                    obj.modeshapes(:,kk+1:end-1)=obj.modeshapes(:,kk+2:end);
            %                end
            %                kk=kk+1;
            %             end
            
            %             for i = 1:size(r,2)
            %                 if isequal(r(i),2)
            %                     y=i;break;
            %                 end
            %             end
            
            %             %script to filter out the 2nd mode from the matrix
            %             k = 1;
            %             for i = 1:size(obj.Mff,2)
            %                 if isequal(r(i),2)
            %
            %                 else
            %                     r1(k)=r(i);          % Eigen modes
            %                     EMa(k)= EMA(i);       % Max Eigen Amplitude
            %                     omegamod(k)=omega(i); % Corresponding Eigen Frequency
            %                     k=k+1;
            %                 end
            %             end
            %
            %             for i = 1:size(r1,2)
            %                 if isequal(r1(i),4)
            %                     r1(i)=1;
            %                 elseif isequal(r1(i),0)
            %                     r1(i)=3;
            %                 else
            %
            %                 end
            %             end
            %
            %             %identify the first position of occurance of mode 5
            %             for i = 1:size(r1,2)
            %                 if isequal(r1(i),5)
            %                     q=i;break;
            %                 end
            %             end
        end
    end
end

