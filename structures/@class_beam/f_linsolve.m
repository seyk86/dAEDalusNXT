%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_linsolve :  file is part of nlFEM class_wing
%               solve linear FEM problem, compute deflections and reaction
%               forces
% Author:       Klaus Seywald
%               klaus.seywald@mytum.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_linsolve(obj)

% isboxwing=obj.isboxwing;
ndof=obj.ndof;
el_ndof=obj.el_ndof;
% Solve equation system


%def=obj.Ks^-1*obj.Qs;

% def=linsolve(obj.Ks,obj.Qs);

obj.wf=linsolve(obj.Kff,(obj.Ff-obj.Kfp*obj.wp));

obj.Frem=obj.Kpp*obj.wp+obj.Kfp'*obj.wf;


% Reaction loads are calculated
% Create result vector containing deflections, rotations and twist

obj.nodal_deflections=zeros(ndof,1);
sort_vec=sortrows(obj.sort_vec,1);
obj.w=[obj.wp' obj.wf']';
obj.w=[obj.w sort_vec(:,2)];
obj.w=sortrows(obj.w,2);
obj.w=obj.w(:,1);

%obj.w=obj.w(:,obj.sort_vec);
% if isboxwing==1
%     obj.nodal_deflections((1+el_ndof):ndof-el_ndof)=def;
% else
%     obj.nodal_deflections((1+el_ndof):ndof)=def;
% end 

obj.nodal_deflections=obj.w;
%compute reaction forces for each node
% obj.reaction_forces=obj.K*obj.nodal_deflections-obj.Q;
% obj.node_loadings=obj.reaction_forces;
% obj.node_loadings_loc=zeros(length(obj.reaction_forces),1);
% %obj.node_loadings(1:6)=obj.reaction_forces(1:6);
% 
% k=1+el_ndof;
% %k=7;
% %el_ndof=6;
% obj.node_loadings_loc(1:el_ndof)=obj.beamelement(1).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(1:el_ndof);
% for i=2:obj.nel
%    el_load=obj.beamelement(i).elKglobal*obj.nodal_deflections(k:k+2*el_ndof-1);
%    obj.node_loadings(k:k+el_ndof-1)=el_load(1:el_ndof);
%    if(el_ndof==6)
%        el_load_loc=obj.beamelement(i).T*el_load;
%        obj.node_loadings_loc(k:k+el_ndof-1)=el_load_loc(1:el_ndof);
%    end
%    k=k+el_ndof;
% end
% obj.node_loadings(k:k+el_ndof-1)=-obj.reaction_forces(k:k+el_ndof-1);
% obj.node_loadings_loc(k:k+el_ndof-1)=-obj.beamelement(end).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(k:k+el_ndof-1);
% %add nodal load later!
% %obj.node_loadings(k:k+el_ndof-1)=[0,0,0,0,0,0];

end

