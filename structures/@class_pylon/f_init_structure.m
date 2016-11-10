%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function[pylon_struct]=f_init_structure(pylon_struct,pylon_geo)

pylon_struct.node_coords(1,1)=pylon_geo.center_coords(1,1)';
pylon_struct.node_coords(1,2)=pylon_geo.center_coords(2,1)';
pylon_struct.node_coords(1,3)=pylon_geo.center_coords(3,1)';

pylon_struct.node_coords(2,1)=pylon_geo.center_coords(1,2)'-0.1;%0.1*pylon_geo.center_coords(1,1)'+0.9*pylon_geo.center_coords(1,2)';
pylon_struct.node_coords(2,2)=pylon_geo.center_coords(2,2)';%0.1*pylon_geo.center_coords(2,1)'+0.9*pylon_geo.center_coords(2,2)';
pylon_struct.node_coords(2,3)=pylon_geo.center_coords(3,2)';%0.1*pylon_geo.center_coords(3,1)'+0.9*pylon_geo.center_coords(3,2)';

pylon_struct.node_coords(3,1)=pylon_geo.center_coords(1,2)';
pylon_struct.node_coords(3,2)=pylon_geo.center_coords(2,2)';
pylon_struct.node_coords(3,3)=pylon_geo.center_coords(3,2)';

pylon_struct.r_Ref(1)=pylon_geo.center_coords(1,1);
pylon_struct.r_Ref(2)=pylon_geo.center_coords(2,1);
pylon_struct.r_Ref(3)=pylon_geo.center_coords(3,1);

%% initialize fuselage element structure
crosssection='na';%class_crosssection_fuselage;

for i=1:2
    elem_vec=pylon_struct.node_coords(i+1,:)-pylon_struct.node_coords(i,:);      
        
    le=sqrt(elem_vec(1)*elem_vec(1)+elem_vec(2)*elem_vec(2)+elem_vec(3)*elem_vec(3));
    if elem_vec(2)==0
        phi=-pi/2;
    end
     phi=pi/2;
    %phi=sign(elem_vec(2))*pi/2;
   % phi=-acos(elem_vec(2)/le);
    nu=0;%-asin(elem_vec(3)/le);
    epsi=0;

    pylon_struct.beamelement(i)=pylon_struct.beamelement(i).setElementGeometry(le,phi,nu,epsi);
    

    epsilon(i)=0;
    nnu(i)=nu;
    nphi(i)=phi;
end
epsilon(i+1)=0;
nnu(i+1)=0;
nphi(i+1)=phi;
% 
 pylon_struct.epsilon=epsilon;
 pylon_struct.nu=nnu;
 pylon_struct.phi=nphi;
end
