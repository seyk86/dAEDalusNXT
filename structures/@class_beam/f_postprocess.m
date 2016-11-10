%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [obj] = f_postprocess(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %compute reaction forces for each node
    obj.reaction_forces=obj.K*obj.nodal_deflections-obj.Q;
    obj.node_loadings=obj.reaction_forces;
    obj.node_loadings_loc=zeros(length(obj.reaction_forces),1);
    %obj.node_loadings(1:6)=obj.reaction_forces(1:6);
    el_ndof=obj.el_ndof;
    k=1;%+el_ndof;
	%k=7;
    %el_ndof=6;
    obj.node_loadings_loc(1:el_ndof)=obj.beamelement(1).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(1:el_ndof);
    if obj.isExternalFEM==0
        for i=1:obj.nel
            el_load=obj.beamelement(i).elKglobal*obj.nodal_deflections(k:k+2*el_ndof-1);
            obj.node_loadings(k:k+el_ndof-1)=el_load(1:el_ndof);
            if(el_ndof==6)
                el_load_loc=obj.beamelement(i).T*el_load;
                obj.node_loadings_loc(k:k+el_ndof-1)=el_load_loc(1:el_ndof);
            end
            k=k+el_ndof;
        end
    end
    obj.node_loadings(k:k+el_ndof-1)=-obj.reaction_forces(k:k+el_ndof-1);
    obj.node_loadings_loc(k:k+el_ndof-1)=-obj.beamelement(end).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(k:k+el_ndof-1);
    
    %add nodal load later!
    %obj.node_loadings(k:k+el_ndof-1)=[0,0,0,0,0,0];
end

