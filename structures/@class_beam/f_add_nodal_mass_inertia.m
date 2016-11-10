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
%     Author:   	Klaus Seywald/Simon Binder
%                   klaus.seywald@mytum.de
%                   seywald@kth.se

function obj = f_add_nodal_mass_inertia(obj, add_eigenmass)

%for every nodal mass which has nonzero weight do
if ~isempty(obj.nodal_masses)
    conm=find(obj.nodal_masses(:,1));
    if add_eigenmass
        for i=1:length(conm)
            m=obj.nodal_masses(conm(i),1);
            dx=obj.nodal_masses(conm(i),2);
            dy=obj.nodal_masses(conm(i),3);
            dz=obj.nodal_masses(conm(i),4);
            if conm(i)==1
                acc_x=obj.beamelement(conm(i)).ax;
                acc_y=obj.beamelement(conm(i)).ay;
                acc_z=obj.beamelement(conm(i)).az;
            elseif conm(i)==obj.nel+1
                acc_x=obj.beamelement(conm(i)-1).ax;
                acc_y=obj.beamelement(conm(i)-1).ay;
                acc_z=obj.beamelement(conm(i)-1).az;
            else
                acc_x=(obj.beamelement(conm(i)).ax+obj.beamelement(conm(i)-1).ax)/2;
                acc_y=(obj.beamelement(conm(i)).ay+obj.beamelement(conm(i)-1).ay)/2;
                acc_z=(obj.beamelement(conm(i)).az+obj.beamelement(conm(i)-1).az)/2;
            end
            Fx=acc_x*m;
            Fy=acc_y*m;
            Fz=acc_z*m;
            Mx=Fz*dy-Fy*dz;
            My=Fx*dz-Fz*dx;
            Mz=Fy*dx-Fx*dy;
            obj.Q(1+6*(conm(i)-1):1+6*(conm(i)-1)+5)=obj.Q(1+6*(conm(i)-1):1+6*(conm(i)-1)+5)+[Fx,Fy,Fz,Mx,My,Mz]';
        end
    end
end
end

