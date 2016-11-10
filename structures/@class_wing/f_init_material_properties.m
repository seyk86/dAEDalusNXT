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
% f_initMaterialProperties :   file is part of nlFEM class_wing
%               initialize FEM model with required material properties
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_init_material_properties(obj,structure)

%%set Material properties
    for i=1:obj.nel
        %% todo now all has Uskin.E and Uskin.G
        obj.beamelement(i).E=structure.materials.E;
        obj.beamelement(i).G=structure.materials.G;
        %Upper Skin
        E_sk_up = structure.materials.E;
        G_sk_up = structure.materials.G;
        rho_sk_up = structure.materials.rho;
        %Lower Skin
        E_sk_lo = structure.materials.E;
        G_sk_lo = structure.materials.G;
        rho_sk_lo = structure.materials.rho;
        %Spars
        rho_sp = structure.materials.rho;
        sigma_allowable_sp = structure.materials.sigma_allowable;
        sigma_allowable_sk_u=structure.materials.sigma_allowable;
        sigma_allowable_sk_l=structure.materials.sigma_allowable;
        
        
        t_min_sp=structure.t_min_sp;
        t_min_sk=structure.t_min_sk;
        
        obj.fuel_density=structure.fuel_density;
        
        obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(E_sk_up,G_sk_up,rho_sk_up,... 
            E_sk_lo,G_sk_lo,rho_sk_lo,rho_sp,sigma_allowable_sp,sigma_allowable_sk_u,sigma_allowable_sk_l,t_min_sp,t_min_sk);
        
    end
end
