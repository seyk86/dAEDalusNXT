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

function obj = f_set_BC(obj, varargin)
% Apply boundary conditions
% Remove locked dofs at x=0

len=length(obj.Q);

bcs=obj.boundary_condition;

sort_vec=zeros(len,2);
wp=0;

node_dof_vec=zeros(len,2);
node_dof_vec(:,1)=reshape(ones(obj.el_ndof,1)*(1:1:len/obj.el_ndof),1,obj.el_ndof*(obj.nel+1));

for i=1:obj.el_ndof:len
    node_dof_vec(i:i+obj.el_ndof-1,2)=1:obj.el_ndof;
end

if (obj.n_bc==0)
    obj.Kff=obj.K;
    obj.Mff=obj.M;
    obj.Mff_lumped=obj.M_lumped;
    
    obj.Ff=obj.Q;
    obj.Pf=obj.P;
    obj.Kff_node_dof_info=node_dof_vec;
    return
end

k=1;

for n=1:1:obj.n_bc
    for i=1:obj.el_ndof
        if bcs(n).dof(i)==1
            sort_idx=obj.el_ndof*(bcs(n).node_idx-1)+i;
            sort_vec(sort_idx,1)=k;
            sort_vec(sort_idx,2)=sort_idx;
            wp(k)=bcs(n).prescribed_val(i);
            k=k+1;
        end
    end
end
j=k;

for i=1:len
    if sort_vec(i,1)==0
        sort_vec(i,1)=j;
        sort_vec(i,2)=i;
        j=j+1;
    end
end

obj.sort_vec=sort_vec;
if isempty(varargin)
    Ksort=[obj.K sort_vec(:,1)];
    Ksort=sortrows(Ksort,len+1);
    Ksort=Ksort(1:len,1:len)';
    Ksort=[Ksort sort_vec(:,1)];
    Ksort=sortrows(Ksort,len+1);
    Ksort=Ksort(1:len,1:len)';

    Msort=[obj.M sort_vec(:,1)];  % Mass Matrix sorting  
    Msort=sortrows(Msort,len+1);
    Msort=Msort(1:len,1:len)';
    Msort=[Msort sort_vec(:,1)];
    Msort=sortrows(Msort,len+1);
    Msort=Msort(1:len,1:len)';

    Msort_lumped=[obj.M_lumped sort_vec(:,1)];  % Mass Matrix sorting  
    Msort_lumped=sortrows(Msort_lumped,len+1);
    Msort_lumped=Msort_lumped(1:len,1:len)';
    Msort_lumped=[Msort_lumped sort_vec(:,1)];
    Msort_lumped=sortrows(Msort_lumped,len+1);
    Msort_lumped=Msort_lumped(1:len,1:len)';


    node_dof_info=sortrows([node_dof_vec obj.sort_vec(:,1)],3);
    obj.Kff_node_dof_info=node_dof_info(k:end,1:2);

    obj.Mpp=Msort(1:k-1,1:k-1);  % Mass part
    obj.Kpp=Ksort(1:k-1,1:k-1);

    obj.Mff=Msort(k:end,k:end); % Mass Part
    obj.Mff_lumped=Msort_lumped(k:end,k:end);
    obj.Kff=Ksort(k:end,k:end);

 
    obj.Mfp=Msort(1:k-1,k:end)';
    obj.Kfp=Ksort(1:k-1,k:end)';

    obj.wp=wp';
end


Qsort=[obj.Q sort_vec(:,1)];
Qsort=sortrows(Qsort,2);
Qsort=Qsort(:,1);
if isempty(obj.wf)
    obj.wf=Qsort(k:end)*0;
end
Psort=[obj.P sort_vec(:,1)];
Psort=sortrows(Psort,2);
Psort=Psort(:,1);

obj.Fp=Qsort(1:k-1);
obj.Ff=Qsort(k:end);
obj.Pf=Psort(k:end);
end

