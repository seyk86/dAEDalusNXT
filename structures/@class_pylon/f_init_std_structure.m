%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% =================================================================
%> @brief initialize beam geometry manually
%>
%> @return initialized instance of the class_beam
% =================================================================

function [obj] = f_init_std_structure(obj,le,phi,nu,epsilon,r,t)

    for k=1:1:obj.nel
         obj.beamelement(k)=obj.beamelement(k).setElementGeometry(le,phi(k),nu(k),epsilon(k));
         obj.beamelement(k).crosssection=obj.beamelement(k).crosssection.setGeometry(r(k),t(k));
         obj.beamelement(k)=obj.beamelement(k).f_calcCrossProp();    
    end

    %calculate coordinates of nodes in xyz system
    obj.node_coords=zeros(length(obj.beamelement),3);

    obj.node_coords(1,1)=obj.r_Ref(1);
    obj.node_coords(1,2)=obj.r_Ref(2);
    obj.node_coords(1,3)=obj.r_Ref(3);

    for i=1:1:length(obj.beamelement)
        obj.node_coords(i+1,1)=obj.node_coords(i,1)+obj.beamelement(i).le*sin(obj.beamelement(i).phi);
        obj.node_coords(i+1,2)=obj.node_coords(i,2)+obj.beamelement(i).le*cos(obj.beamelement(i).nu)*cos(obj.beamelement(i).phi);
        obj.node_coords(i+1,3)=obj.node_coords(i,3)+obj.beamelement(i).le*sin(obj.beamelement(i).nu);
    end
end
