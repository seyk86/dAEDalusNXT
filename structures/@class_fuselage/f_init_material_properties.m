%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj = f_init_material_properties(obj,structure)

%%set Material properties
    for i=1:obj.nel
        %% todo now all has Uskin.E and Uskin.G
        obj.beamelement(i).E=structure.materials.E;
        obj.beamelement(i).G=structure.materials.G;
        %Equivalent Skin
        E_sk_eq = structure.materials.E;
        G_sk_eq = structure.materials.G;
        rho_sk_eq = structure.materials.rho;
        %Frames
        E_fr_eq = structure.materials.E;
        G_fr_eq = structure.materials.G;
        rho_fr_eq = structure.materials.rho;
        %Allowables
        sigma_allowable_sk = structure.materials.sigma_allowable;
        sigma_allowable_fr=structure.materials.sigma_allowable;      
        t_min_sk=structure.t_min_sk;
        
        obj.fuel_density=structure.fuel_density;
        
        obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(E_sk_eq,G_sk_eq,rho_sk_eq,... 
            E_fr_eq,G_fr_eq,rho_fr_eq,sigma_allowable_sk,sigma_allowable_fr,t_min_sk);
        
    end
end
