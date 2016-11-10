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
% f_assemble:       Assemble system stiffness matrix and load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_assemble(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces)
    %% initialize variables
    
    el_ndof=obj.el_ndof;
    k=1; %counter
    
    if obj.update_K==1
        if ~obj.isExternalFEM
            obj.K=obj.K*0;
            obj=obj.f_calc_element_stiffness();
            
            for i=1:obj.nel  
                %Assemble System Stiffness Matrix
                obj.K(k:k+2*el_ndof-1,k:k+2*el_ndof-1)=...
                    obj.K(k:k+2*el_ndof-1,k:k+2*el_ndof-1)...
                    +obj.beamelement(i).elKglobal;
                k=k+el_ndof;
            end 
        else
            obj.K = obj.externalFEM.Kext;
            obj.Kff = obj.externalFEM.Kext;
        end
        obj.update_K=0;
    end
    
    k=1; %counter
    if obj.update_M==1
        if ~obj.isExternalFEM
            obj.M=zeros(size(obj.M,1),size(obj.M,2));
            obj=obj.f_calc_element_mass();
            
            for i=1:obj.nel
                %Assemble System mass Matrix
                obj.M(k:k+2*el_ndof-1,k:k+2*el_ndof-1)=...
                    obj.M(k:k+2*el_ndof-1,k:k+2*el_ndof-1)...
                    +obj.beamelement(i).elMglobal;
                k=k+el_ndof;
            end
        
            k=1; %counter
            obj.M_lumped=obj.M_lumped*0;
            for i=1:obj.nel
                %Assemble System mass Matrix
                obj.M_lumped(k:k+2*el_ndof-1,k:k+2*el_ndof-1)=...
                    obj.M_lumped(k:k+2*el_ndof-1,k:k+2*el_ndof-1)...
                    +obj.beamelement(i).elM_lumped_global;
                k=k+el_ndof;
            end
        else
            obj.M = obj.externalFEM.Mext;
            obj.Mff = obj.externalFEM.Mext;
            obj.Mff_lumped = obj.externalFEM.Mext;
            obj.M_lumped = obj.externalFEM.Mext;
        end
%         if any(any(obj.nodal_masses))
%            obj=obj.f_add_nodal_masses();
%         end
        obj.update_M=0;
    end
    
    k=1; %counter
    if obj.update_Q==1
         obj=obj.f_calc_element_loadvec(add_eigenmass,add_fuelmass);
         obj.update_Q=0;
         obj.Q=obj.Q*0;
        %% assemble system matrix and load vector
        for i=1:obj.nel  
            %Assemble System Load Vector
            obj.Q(k:k+2*el_ndof-1)=obj.Q(k:k+2*el_ndof-1)...
                +obj.beamelement(i).elqglobal;
            k=k+el_ndof;
        end 
        if obj.isExternalFEM == 0
%             obj=obj.f_add_inertial_forces(add_engineforces,add_gearforces);
            obj=obj.f_add_nodal_mass_inertia(add_eigenmass);
        elseif obj.isExternalFEM == 1
            g = [obj.beamelement(1).ax;obj.beamelement(1).ay;obj.beamelement(1).az;0;0;0];
            a = repmat(g,(length(obj.M)/6),1);
            obj.Q = obj.Q + obj.M*a;
        end
    
        obj=obj.f_add_nodal_loads();
    end
    

    
end
