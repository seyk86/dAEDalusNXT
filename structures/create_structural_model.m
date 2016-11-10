%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,aircraft_structure]=create_structural_model(aircraft)

for nwings=1:length(aircraft.wings_structural_properties)
    wing_settings=class_structural_settings_wing;
    %TODO: different material for each segment
    material=class_material(aircraft.wings_structural_properties(nwings).material(1));
    wing_settings=wing_settings.f_set_material(material);
    aircraft=aircraft.compute_wingbox_coords();
    wingstructure=class_wing((size(aircraft.wings(nwings).wingbox_coords,2)-1)*(aircraft.wings(nwings).symmetric+1),'class_crosssection_wingbox',aircraft.wings(nwings).name);
    wingstructure=wingstructure.f_init_structure(aircraft.wings(nwings),aircraft.wings_structural_properties(nwings).is_fueled);
    if wingstructure.isExternalFEM==0
        wingstructure=wingstructure.f_init_material_properties(wing_settings);
    elseif wingstructure.isExternalFEM==1
        % TODO: do material properties and allowed stresses need to be given?
        
        % Load external FEM mass and stiffness matrices
        iWing = [];
        nameTemp = wingstructure.identifier;
        for j=1:length(aircraft.wings)
            if strcmp(nameTemp, aircraft.wings(j).name)
                iWing = j;
            end
        end
        
        wingstructure.externalFEM.Mext = importdata(aircraft.wings(iWing).pathMassMatrix);
        wingstructure.externalFEM.Kext = importdata(aircraft.wings(iWing).pathStiffnessMatrix);
        
        wingstructure.update_M=1;
        wingstructure.update_K=1;
        wingstructure.update_Q=1;
        
        for j=1:length(wingstructure.beamelement)
            wingstructure.beamelement(j).m = 0;
        end
    end
    for nengines=1:length(aircraft.engines)
        if strcmp(aircraft.engines(nengines).mounting,aircraft.wings(nwings).name)
            wingstructure=wingstructure.f_add_engine(aircraft.engines(nengines));
            wingstructure.wingmountedengines=1;
        end
    end 
    if nwings==1
        aircraft_structure=class_beam_collection(wingstructure);
    else
        aircraft_structure=aircraft_structure.f_add_beam(wingstructure);
    end
end

for nfuse=1:length(aircraft.fuselages)
    fuselage_settings=class_structural_settings_fuselage;  
    fuselage_settings.delta_pressure=8.9632e+04;
    fuselage_settings=fuselage_settings.f_set_material(class_material('aluminum'));
    aircraft=aircraft.compute_shell_coords();
    aircraft.fuselages(nfuse)=aircraft.fuselages(nfuse).compute_shell_coords(aircraft.grid_settings.y_max_grid_size);
    fuselage_structure=class_fuselage(size(aircraft.fuselages(nfuse).center_coords,2)-1,'class_crosssection_fuselage',aircraft.fuselages(nfuse).name);
    
    fuselage_structure=fuselage_structure.f_init_structure(aircraft.fuselages(nfuse));
    for i=1:size(fuselage_structure.beamelement,2)
        fuselage_structure.beamelement(i).crosssection.delta_pressure=fuselage_settings.delta_pressure;
    end
    
    if fuselage_structure.isExternalFEM==0
        fuselage_structure=fuselage_structure.f_init_material_properties(fuselage_settings);
    elseif fuselage_structure.isExternalFEM==1
        % TODO: do material properties and allowed stresses need to be given?
        
        % Load external FEM mass and stiffness matrices
        iFuselage = [];
        nameTemp = fuselage_structure.identifier;
        for j=1:length(aircraft.fuselages)
            if strcmp(nameTemp, aircraft.fuselages(j).name)
                iFuselage = j;
            end
        end
        
        fuselage_structure.externalFEM.Mext = importdata(aircraft.fuselages(iFuselage).pathMassMatrix);
        fuselage_structure.externalFEM.Kext = importdata(aircraft.fuselages(iFuselage).pathStiffnessMatrix);
        
        fuselage_structure.update_M=1;
        fuselage_structure.update_K=1;
        fuselage_structure.update_Q=1;
        
        for j=1:length(fuselage_structure.beamelement)
            fuselage_structure.beamelement(j).m = 0;
        end
    end
    
   % fuselage_structure=fuselage_structure.f_add_boundary_condition(class_boundary_condition(1,[1 1 1 1 1 1],[0 0 0 0 0 0]));
   %TODO use xml specified constrainment instead of the following line in which the fuselage is constrained close to the leading edge of the main wing (wing(1))
   
%    fuselage_structure=fuselage_structure.f_add_boundary_condition(class_boundary_condition(ceil((abs(aircraft.fuselages(1).fuselage_segments(1,1).pos(1,1))+aircraft.wings(1,1).wing_segments(1,1).pos(1,1))/aircraft.grid_settings.dy_max_struct_grid),[1 1 1 1 1 1],[0 0 0 0 0 0]));
   wingRootCoordX = aircraft_structure.beam(1).node_coords(ceil(aircraft_structure.beam(1).nel/2)+1,1);
   distanceToWing = fuselage_structure.node_coords(:,1)-wingRootCoordX;
   [~,boundaryIdx] = min(abs(distanceToWing));
   fuselage_structure = fuselage_structure.f_add_boundary_condition(class_boundary_condition(boundaryIdx,[1 1 1 1 1 1],[0 0 0 0 0 0]));

   aircraft_structure=aircraft_structure.f_add_beam(fuselage_structure);
end

n_cc=length(aircraft.boundary_conditions);
for ne=1:length(aircraft.engines)
    %find closest node:
    for nw=1:length(aircraft.wings)
        if(strcmp(aircraft.engines(ne).mounting,aircraft.wings(nw).name))
            dist=zeros(length(aircraft_structure.beam(nw).node_coords(:,2)),1);
            for i=1:length(aircraft_structure.beam(nw).node_coords(:,2))
                dist(i)=sqrt(sum((aircraft_structure.beam(nw).node_coords(i,2)-aircraft.engines(ne).cg_pos(2)').^2));
            end
            [Y,node_idx_beam] = min(dist);
            %add forces and moments
            pyl_coords=[aircraft_structure.beam(nw).node_coords(node_idx_beam,:)' aircraft.engines(ne).cg_pos'];
            pylon=class_pylon(2,'none',['Engine' num2str(ne) 'Pylon' ]);
            geof.center_coords=pyl_coords;
            pylon=pylon.f_init_structure(geof);
            pylon.beamelement(1).el_m_sys=0.001;
            pylon.beamelement(2).el_m_sys=aircraft.engines(ne).m-0.001;
            connection{1}=aircraft.engines(ne).mounting;
            connection{2}=['Engine' num2str(ne) 'Pylon' ];
            aircraft_structure=aircraft_structure.f_add_beam(pylon);
            aircraft.boundary_conditions{n_cc+ne}=connection;
        end
    end
    for nf=1:length(aircraft.fuselages)
        if(strcmp(aircraft.engines(ne).mounting,aircraft.fuselages(nf).name))
            dist=zeros(length(aircraft_structure.beam(nf+nw).node_coords(:,1)),1);
            for i=1:length(aircraft_structure.beam(nf+nw).node_coords(:,1))
                dist(i)=sqrt(sum((aircraft_structure.beam(nf+nw).node_coords(i,1)-aircraft.engines(ne).cg_pos(1)').^2));
            end
            [Y,node_idx_beam] = min(dist);
            %add forces and moments
            pyl_coords=[aircraft_structure.beam(nf+nw).node_coords(node_idx_beam,:)' aircraft.engines(ne).cg_pos'];
            pylon=class_pylon(2,'none',['Engine' num2str(ne) 'Pylon' ]);
            geof.center_coords=pyl_coords;
            pylon=pylon.f_init_structure(geof);
            pylon.beamelement(1).el_m_sys=0.001;
            pylon.beamelement(2).el_m_sys=aircraft.engines(ne).m-0.001;
            connection{1}=aircraft.engines(ne).mounting;
            connection{2}=['Engine' num2str(ne) 'Pylon' ];
            aircraft_structure=aircraft_structure.f_add_beam(pylon);
            aircraft.boundary_conditions{n_cc+ne}=connection;
        end
    end
end

for nbc=1:length(aircraft.boundary_conditions)
    balljoint=0;
    hingejoint=0;
    zjoint=0;
    for nbeam=1:length(aircraft_structure.beam)
        if strcmp(aircraft.boundary_conditions{nbc}{1},aircraft_structure.beam(nbeam).identifier)
           idx_1=nbeam; 
        end
        if strcmp(aircraft.boundary_conditions{nbc}{2},aircraft_structure.beam(nbeam).identifier)
           idx_2=nbeam; 
        end
        if length(aircraft.boundary_conditions{nbc})==3
            if strcmp(aircraft.boundary_conditions{nbc}{3},'BallJoint')
                 balljoint=1;
            elseif strcmp(aircraft.boundary_conditions{nbc}{3},'ZJoint')
                 zjoint=1;
            elseif strcmp(aircraft.boundary_conditions{nbc}{3},'KardanJoint')
                 hingejoint=1;
            end
        end
    end
    dis=1E6;
    dis_min=dis;
    node_idx_1=0;
    node_idx_2=0;
    for nc_ctr1=1:length(aircraft_structure.beam(idx_1).node_coords)
         for nc_ctr2=1:length(aircraft_structure.beam(idx_2).node_coords)
            dis=norm(aircraft_structure.beam(idx_1).node_coords(nc_ctr1,:)'-aircraft_structure.beam(idx_2).node_coords(nc_ctr2,:)');
            if dis<dis_min
               dis_min=dis;
               node_idx_1=nc_ctr1;
               node_idx_2=nc_ctr2;
            elseif abs(dis-dis_min)<1E-5
               node_idx_1=[node_idx_1 nc_ctr1];
               node_idx_2=[node_idx_2 nc_ctr2];  
            end
         end
    end
    for nnodes=1:length(node_idx_1)
        if balljoint==0 && zjoint==0
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 1 1 1]);
        elseif balljoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 0 1 1]); 
        elseif zjoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[0 0 1 0 0 0]); 
        elseif hingejoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 0 1 1]); 
        end
    end
end

% pre-merge coupling conditions
del_idx=[];
for i=1:length(aircraft_structure.coupling_condition)-1
   for j=i+1:length(aircraft_structure.coupling_condition)
       for kk=1:min([length(aircraft_structure.coupling_condition(i).beam_idx) length(aircraft_structure.coupling_condition(j).beam_idx)])
            if (aircraft_structure.coupling_condition(i).beam_idx(kk)==aircraft_structure.coupling_condition(j).beam_idx(kk))&&(aircraft_structure.coupling_condition(i).node_idx(kk)==aircraft_structure.coupling_condition(j).node_idx(kk))
                if kk==1
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(2)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(2)];
                 %   aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                elseif kk==2
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(1)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(1)]; 
                   % aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                end
                del_idx=[del_idx j];
            end
       end
   end
end

del_idx=unique(del_idx);
if ~isempty(del_idx)
    temp_cc=[];
    jj=1;
    for i=1:length(aircraft_structure.coupling_condition)
        if i==del_idx(jj)
            if jj<length(del_idx)
                jj=jj+1;
            end
        else
            temp_cc=[temp_cc aircraft_structure.coupling_condition(i)];
        end
    end

    aircraft_structure.coupling_condition=temp_cc;
end
aircraft_structure=aircraft_structure.f_calc_mass(aircraft.weights);
aircraft_structure.identifier=aircraft.name;

% pre-merge coupling conditions
del_idx=[];
for i=1:length(aircraft_structure.coupling_condition)-1
   for j=i+1:length(aircraft_structure.coupling_condition)
       for kk=1:min([length(aircraft_structure.coupling_condition(i).beam_idx) length(aircraft_structure.coupling_condition(j).beam_idx)])
            if (aircraft_structure.coupling_condition(i).beam_idx(kk)==aircraft_structure.coupling_condition(j).beam_idx(kk))&&(aircraft_structure.coupling_condition(i).node_idx(kk)==aircraft_structure.coupling_condition(j).node_idx(kk))
                if kk==1
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(2)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(2)];
                 %   aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                elseif kk==2
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(1)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(1)]; 
                   % aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                end
                del_idx=[del_idx j];
            end
       end
   end
end

del_idx=unique(del_idx);
if ~isempty(del_idx)
    temp_cc=[];
    jj=1;
    for i=1:length(aircraft_structure.coupling_condition)
        if i==del_idx(jj)
            if jj<length(del_idx)
                jj=jj+1;
            end
        else
            temp_cc=[temp_cc aircraft_structure.coupling_condition(i)];
        end
    end

    aircraft_structure.coupling_condition=temp_cc;
end
aircraft_structure=aircraft_structure.f_calc_mass(aircraft.weights);
