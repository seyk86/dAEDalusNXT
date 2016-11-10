%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus structures
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se
%  f_set_BC :       set boundary conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_set_no_BC(obj)
% Apply boundary conditions
% Remove locked dofs at x=0

len=length(obj.Q);

node_dof_vec=zeros(len,2);
node_dof_vec(:,1)=reshape(ones(obj.el_ndof,1)*(1:1:len/obj.el_ndof),1,obj.el_ndof*(obj.nel+1));

for i=1:obj.el_ndof:len
    node_dof_vec(i:i+obj.el_ndof-1,2)=1:obj.el_ndof;
end

    obj.Kff=obj.K;
    obj.Mff=obj.M;
    obj.Mff_lumped=obj.M_lumped;
        
    obj.Ff=obj.Q;
    obj.Pf=obj.P;
    obj.Kff_node_dof_info=node_dof_vec;

end

