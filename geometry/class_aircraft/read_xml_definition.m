function obj=read_xml_definition(obj,filename)
wings_structural_properties=[];
xmlstruct=xmltools(filename);
if strcmp(xmlstruct.child(1).tag,'DAEDALUS')
    xml_aircraft=xmlstruct.child(1).child;
    obj.name= xml_aircraft.attribs(1).value;
    for i=1:length(xml_aircraft.child)
        if strcmp(xml_aircraft.child(i).tag,'WINGS')
            for j=1:length(xml_aircraft.child(i).child)
                frontspar=[];
                rearspar=[];
                is_fueled=[];
                material={};
                if j==1
                    obj.wings=class_aerosurface(xml_aircraft.child(i).child(j));
                else
                    obj.wings(j)=class_aerosurface(xml_aircraft.child(i).child(j));
                end
                for seg=1:length(xml_aircraft.child(i).child(j).child)
                    factor=1;
                    for wb=1:length(xml_aircraft.child(i).child(j).child(seg).child)
                        if strcmp(xml_aircraft.child(i).child(j).child(seg).child(wb).tag,'CONTROL_SURFACE')
                            if ~(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value)==0)
                                factor=[factor 1];
                            end
                            if ~(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value)==1)
                                factor=[factor 1];
                            end
                        end
                        if strcmp(xml_aircraft.child(i).child(j).child(seg).child(wb).tag,'WINGBOX')
                            if length(factor)>1
                                root_fs=str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(1).value);
                                tip_fs=str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value);
                                root_rs=str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value);
                                tip_rs=str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(4).value);
                                
                                if ~(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value)==0)
                                    inner_fs=(tip_fs)*(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value))+(root_fs)*(1-str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value));
                                    inner_rs=(tip_rs)*(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value))+(root_rs)*(1-str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value));
                                else
                                    inner_fs=[];
                                    inner_rs=[];
                                end
                                if ~(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value)==1)
                                    outer_fs=(tip_fs)*(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value))+(root_fs)*(1-str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value));
                                    outer_rs=(tip_rs)*(str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value))+(root_rs)*(1-str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value));
                                else
                                    outer_fs=[];
                                    outer_rs=[];
                                end
                                frontspar=[frontspar str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(1).value) inner_fs inner_fs  outer_fs outer_fs str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value)];
                                rearspar=[rearspar str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value) inner_rs inner_rs outer_rs outer_rs str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(4).value)];
                            else
                                frontspar=[frontspar str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(1).value) str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(2).value)];
                                rearspar=[rearspar str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(3).value) str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(4).value)];
                            end
                            is_fueled=[is_fueled factor*str2double(xml_aircraft.child(i).child(j).child(seg).child(wb).child(5).value)];
                            material{seg}= xml_aircraft.child(i).child(j).child(seg).child(wb).child(6).attribs.value;
                        end
                    end
                end
                
                wings_structural_properties(j).frontspar=frontspar;
                wings_structural_properties(j).rearspar=rearspar;
                wings_structural_properties(j).is_fueled=is_fueled;
                wings_structural_properties(j).material=material;
            end
        end
        
        if strcmp(xml_aircraft.child(i).tag,'FUSELAGES')
            for j=1:length(xml_aircraft.child(i).child)
                if j==1
                    obj.fuselages=class_fuselage_geometry(xml_aircraft.child(i).child(j));
                else
                    obj.fuselages(j)=class_fuselage_geometry(xml_aircraft.child(i).child(j));
                end
            end
        end
        
        if strcmp(xml_aircraft.child(i).tag,'NACELLES')
            for j=1:length(xml_aircraft.child(i).child)
                if j==1
                    obj.nacelles=class_fuselage_geometry(xml_aircraft.child(i).child(j));
                else
                    obj.nacelles(j)=class_fuselage_geometry(xml_aircraft.child(i).child(j));
                end
            end
        end
        
        if strcmp(xml_aircraft.child(i).tag,'ENGINES')
            for j=1:length(xml_aircraft.child(i).child)
                engine=class_engine;
                engine.cg_pos=[str2double(xml_aircraft.child(i).child(j).child(1).child(1).value) str2double(xml_aircraft.child(i).child(j).child(1).child(2).value) str2double(xml_aircraft.child(i).child(j).child(1).child(3).value)];
                engine.mass_engine=str2double(xml_aircraft.child(i).child(j).child(2).value);
                engine.mass_pylon=str2double(xml_aircraft.child(i).child(j).child(3).value);
                engine.m=engine.mass_engine+engine.mass_pylon;
                engine.thrust=str2double(xml_aircraft.child(i).child(j).child(4).value);
                engine.mounting=xml_aircraft.child(i).child(j).attribs(2).value;
                
                if j==1
                    obj.engines=engine;
                else
                    obj.engines(j)=engine;
                end
            end
        end
        
        if strcmp(xml_aircraft.child(i).tag,'WEIGHTS')
            weights=class_weights;
            weights.WingSystemsEstimate=[];
            weights.WingSkinEstimate=[];
            weights.FuselageNonStructuralEstimate=[];
            weights.FuselageSystemsEstimate=[];
            weights.WingInitialGuess=[];
            for j=1:length(xml_aircraft.child(i).child)
                if strcmp(xml_aircraft.child(i).child(j).tag,'GENERAL_WEIGHTS')
                    weights.MTOW=str2num(xml_aircraft.child(i).child(j).child(1).value);
                    weights.MZFW=str2num(xml_aircraft.child(i).child(j).child(2).value);
                    weights.MLW=str2num(xml_aircraft.child(i).child(j).child(3).value);
                    weights.OWE=str2num(xml_aircraft.child(i).child(j).child(4).value);
                end
                if strcmp(xml_aircraft.child(i).child(j).tag,'WING_WEIGHTS')
                    for jkj=1:length(xml_aircraft.child(i).child(j).child)
                        weights.WingSystemsEstimate=[weights.WingSystemsEstimate str2num(xml_aircraft.child(i).child(j).child(jkj).child(1).value)];
                        weights.WingSkinEstimate=[weights.WingSkinEstimate str2num(xml_aircraft.child(i).child(j).child(jkj).child(2).value)];
                        weights.WingInitialGuess=[weights.WingInitialGuess 0];
                    end
                end
                if strcmp(xml_aircraft.child(i).child(j).tag,'FUSELAGE_WEIGHTS')
                    for jkj=1:length(xml_aircraft.child(i).child(j).child)
                        weights.FuselageNonStructuralEstimate=[weights.FuselageNonStructuralEstimate str2num(xml_aircraft.child(i).child(j).child(jkj).child(1).value)];
                        weights.FuselageSystemsEstimate=[weights.FuselageSystemsEstimate str2num(xml_aircraft.child(i).child(j).child(jkj).child(2).value)];
                        weights.FuselagePAX=[weights.FuselagePAX str2num(xml_aircraft.child(i).child(j).child(jkj).child(3).value)];
                    end
                end
            end
            obj.weights=weights;
        end
        
        
        if strcmp(xml_aircraft.child(i).tag,'BOUNDARY_CONDITIONS')
            for j=1:length(xml_aircraft.child(i).child)
                if strcmp(xml_aircraft.child(i).child(j).tag,'CONNECTION')
                    connection{1}=xml_aircraft.child(i).child(j).child(1).attribs(1).value;
                    connection{2}=xml_aircraft.child(i).child(j).child(2).attribs(1).value;
                    if length(xml_aircraft.child(i).child(j).child)==3
                        connection{3}=xml_aircraft.child(i).child(j).child(3).attribs(1).value;
                    end
                    obj.boundary_conditions{j}=connection;
                end
            end
        end
        
        if strcmp(xml_aircraft.child(i).tag,'REFERENCE')
            obj.reference=class_reference(str2double(xml_aircraft.child(i).child(1).value),str2double(xml_aircraft.child(i).child(2).value),str2double(xml_aircraft.child(i).child(3).value),...
                [str2double(xml_aircraft.child(i).child(4).child(1).value),str2double(xml_aircraft.child(i).child(4).child(2).value),str2double(xml_aircraft.child(i).child(4).child(3).value)]);
        end
        

        
        %                   if strcmp(xml_aircraft.child(i).tag,'PROFILES')
        %                       for j=1:length(xml_aircraft.child(i).child)
        %                           obj.wings=class_aerosurface(xml_aircraft.child(i).child(j));
        %                       end
        %                   end
        %%TODO Fuselage and so on!
    end
else
    
    fprintf('Unknown Data Format');
end
obj.wings_structural_properties=wings_structural_properties;
k=1;
for i=1:length(obj.wings)
    for j=1:length(obj.wings(i).wing_segments)
        if ~isempty(obj.wings(i).wing_segments(j).te_device)
            
            if(obj.wings(i).wing_segments(j).te_device.is_sym_defl==1)
                obj.control_surfaces{k}=obj.wings(i).wing_segments(j).te_device.name;
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).te_device.delta;
            else
                obj.control_surfaces{k}=[obj.wings(i).wing_segments(j).te_device.name '_left'];
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).te_device.delta;
                k=k+1;
                obj.control_surfaces{k}=[obj.wings(i).wing_segments(j).te_device.name '_right'];
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).te_device.delta;
            end
            k=k+1;
        end
        if ~isempty(obj.wings(i).wing_segments(j).le_device)
            if(obj.wings(i).wing_segments(j).le_device.is_sym_defl==1)
                obj.control_surfaces{k}=obj.wings(i).wing_segments(j).le_device.name;
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).le_device.delta;
            else
                obj.control_surfaces{k}=[obj.wings(i).wing_segments(j).le_device.name '_left'];
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).le_device.delta;
                k=k+1;
                obj.control_surfaces{k}=[obj.wings(i).wing_segments(j).le_device.name '_right'];
                obj.control_deflections{k}=obj.wings(i).wing_segments(j).le_device.delta;
            end
            k=k+1;
        end
    end
end
obj.control.trim_surface_idx=[];

obj.control.lat_names=[];
obj.control.lat_ctrl_idx=[];
obj.control.lon_names=[];
n_eng=length(obj.engines);
obj.control.lon_ctrl_idx=1:n_eng;
for eng_ctr=1:n_eng
    obj.control.lon_names=[ obj.control.lon_names {['engine' num2str(eng_ctr)]}];
end


obj.control.control_allocation_matrix=zeros(length(obj.control_surfaces),4);
for i=1:length(xml_aircraft.child)
       if strcmp(xml_aircraft.child(i).tag,'CONTROL')
           for j=1:length(xml_aircraft.child(i).child)
                if strcmp(xml_aircraft.child(i).child(j).tag,'CONTROL_ALLOCATION')
                    for k=1:length(xml_aircraft.child(i).child(j).child)
                        ts_ctr=1;

                        if strcmp(xml_aircraft.child(i).child(j).child(k).tag,'TRIM_SURFACE')
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp(xml_aircraft.child(i).child(j).child(k).attribs(1).value,obj.control_surfaces{cs_ctr})
                                     obj.control.trim_surfaces{ts_ctr}=xml_aircraft.child(i).child(j).child(k).attribs(1).value;
                                     obj.control.control_allocation_matrix(cs_ctr,4)=1;
                                     ts_ctr=ts_ctr+1;
                                      obj.control.trim_surface_idx=[obj.control.trim_surface_idx cs_ctr];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_left'],obj.control_surfaces{cs_ctr})
                                     obj.control.trim_surfaces{ts_ctr}=[xml_aircraft.child(i).child(j).child(k).attribs(1).value  '_left'];
                                     obj.control.control_allocation_matrix(cs_ctr,4)=-1;
                                     ts_ctr=ts_ctr+1;
                                     obj.control.trim_surface_idx=[obj.control.trim_surface_idx cs_ctr];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_right'],obj.control_surfaces{cs_ctr})
                                     obj.control.trim_surfaces{ts_ctr}=[xml_aircraft.child(i).child(j).child(k).attribs(1).value '_right'];
                                     obj.control.control_allocation_matrix(cs_ctr,4)=1;
                                     ts_ctr=ts_ctr+1;
                                     obj.control.trim_surface_idx=[obj.control.trim_surface_idx cs_ctr];
                                end
                            end
                        end
                     
                        if strcmp(xml_aircraft.child(i).child(j).child(k).tag,'PITCH_CMD')
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp(xml_aircraft.child(i).child(j).child(k).attribs(1).value,obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,1)=1;
                                    obj.control.lon_ctrl_idx=[obj.control.lon_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lon_names=[ obj.control.lon_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_left'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,1)=-1;
                                    obj.control.lon_ctrl_idx=[obj.control.lon_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lon_names=[ obj.control.lon_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_right'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,1)=1;
                                    obj.control.lon_ctrl_idx=[obj.control.lon_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lon_names=[ obj.control.lon_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                        end
                        
                        if strcmp(xml_aircraft.child(i).child(j).child(k).tag,'ROLL_CMD')
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp(xml_aircraft.child(i).child(j).child(k).attribs(1).value,obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,2)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_left'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,2)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_right'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,2)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            
                        end
                        
                        if strcmp(xml_aircraft.child(i).child(j).child(k).tag,'YAW_CMD')
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp(xml_aircraft.child(i).child(j).child(k).attribs(1).value,obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,3)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_left'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,3)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            for cs_ctr=1:length(obj.control_surfaces)
                                if strcmp([xml_aircraft.child(i).child(j).child(k).attribs(1).value '_right'],obj.control_surfaces{cs_ctr})
                                    obj.control.control_allocation_matrix(cs_ctr,3)=1;
                                    obj.control.lat_ctrl_idx=[obj.control.lat_ctrl_idx cs_ctr+n_eng];
                                    obj.control.lat_names=[ obj.control.lat_names {obj.control_surfaces{cs_ctr}}];
                                end
                            end
                            
                        end    
                    end
                end
           end
       end
        
end
obj.control.control_allocation_matrix_lon=obj.control.control_allocation_matrix(obj.control.lon_ctrl_idx(n_eng+1:end)-n_eng,:);
obj.control.control_allocation_matrix_lon(:,[2 3])=0;
obj.control.control_allocation_matrix_lat=obj.control.control_allocation_matrix(obj.control.lat_ctrl_idx(1:end)-n_eng,:);
obj.control.control_allocation_matrix_lat(:,[1 4])=0;

end