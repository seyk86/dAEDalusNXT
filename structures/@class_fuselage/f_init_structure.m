%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function[fuselage_struct]=f_init_structure(fuselage_struct,fuselage_geo)

fuselage_struct.node_coords(:,1)=fuselage_geo.center_coords(1,:)';
fuselage_struct.node_coords(:,2)=fuselage_geo.center_coords(2,:)';
fuselage_struct.node_coords(:,3)=fuselage_geo.center_coords(3,:)';

fuselage_struct.r_Ref(1)=fuselage_geo.center_coords(1,1);
fuselage_struct.r_Ref(2)=fuselage_geo.center_coords(2,1);
fuselage_struct.r_Ref(3)=fuselage_geo.center_coords(3,1);

%% initialize fuselage element structure
crosssection=class_crosssection_fuselage;

% Identification of which fuselage segment the element belongs
seg_id=[];
for i=1:size(fuselage_geo.fuselage_segments,2)
    n=size(fuselage_geo.fuselage_segments(i).center_coords,2);
    seg_id=[seg_id i*ones(1,n-1)];
end

for i=1:size(fuselage_geo.center_coords,2)-1
    w=0.5*(fuselage_geo.shell_width(i)+fuselage_geo.shell_width(i+1));
    h=0.5*(fuselage_geo.shell_height(i)+fuselage_geo.shell_height(i+1));
    % Equivalent radius for a circle with the same area of an
    % ellipse defined by w and h
    r=(w*h)^0.5;
    fuselage_struct.beamelement(i).crosssection=crosssection.setGeometry(r);
    
    elem_vec=fuselage_geo.center_coords(:,i+1)-fuselage_geo.center_coords(:,i);      
        
    le=sqrt(elem_vec(1)*elem_vec(1)+elem_vec(2)*elem_vec(2)+elem_vec(3)*elem_vec(3));
    phi=-asin(elem_vec(1)/le);
    nu=0;
    sweep=fuselage_geo.fuselage_segments(seg_id(i)).sweep*pi/180;

    fuselage_struct.beamelement(i)=fuselage_struct.beamelement(i).setElementGeometry(le,phi,nu,sweep);
    

        epsilon(i)=sweep;
        nnu(i)=0;
        nphi(i)=phi;

    %.setElementGeometry(le,phi,nu,twist);
    fuselage_struct.Vfuselage=fuselage_struct.Vfuselage+pi*r^2*le;
end
        epsilon(i+1)=sweep;
        nnu(i+1)=0;
        nphi(i+1)=phi;
% 
 fuselage_struct.epsilon=epsilon;
 fuselage_struct.nu=nnu;
 fuselage_struct.phi=nphi;
 
 fuselage_struct.isExternalFEM = fuselage_geo.isExternalFEM;
end
