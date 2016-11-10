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
%     elM_6dof:     class method of class_beamelement
%                   calculates 6DOF mass matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj=elM_lumped_6dof(obj)
   
   rho=obj.m/obj.A;
   if isa(obj.crosssection, 'class_crosssection_wingbox')
       rho_t=obj.crosssection.rho_sp;
   elseif isa(obj.crosssection, 'class_crosssection_fuselage')
       rho_t=obj.crosssection.rho_sk_eq;
   else
       rho_t=obj.m/obj.A;
   end
   massless=0;
   if massless
       obj.m=0;
       rho_t=obj.m/obj.A;
   end
   
    L=obj.le;
 %   m=obj.m*obj.le;
    J=obj.J;
    
    
    if obj.is_fueled
        rho_fuel=obj.ref_parent_beam.fuel_density;
        obj.el_fuel_vol=obj.A_enclosed*obj.le;
        m=obj.m+obj.A_enclosed*rho_fuel;
    else
        m=obj.m;
    end

    m=m*obj.le;
    %COMMENT SIMON:
    %in the following the modified density rho (which includes 
    %nonstructural masses) is used for calculation of the torsional inertia. 
    %when comparing with nastran change rho*J*L to "rho_t*J*L" in order to get
    %the same torsional mode frequencies
    % Calculate Lumped Element Mass Matrices in Local Coordinates
    Q1=[1,      0,      0,      0,      0,  0;
        0,      1,      0,      0,      0,  0;
        0,      0,      1,      0,      0,  0;
        0,      0,      0,      0,      0,  0;
        0,      0,      0,      0,      (rho_t*J*L)*(1/2)*2/m,  0;
        0,      0,      0,      0,      0,  0;];
    
    Q2=[0,     0,      0,      0,      0,  0;
        0,     0,      0,      0,      0,  0;
        0,     0,      0,      0,      0,  0;
        0,     0,      0,      0,      0,  0;
        0,     0,      0,      0,      0,  0;
        0,     0,      0,      0,      0,  0;];
    
    Q4=[1,      0,      0,      0,      0,  0;
        0,      1       0,      0,      0,  0;
        0,      0,      1,      0,      0,  0;
        0,      0,      0,      0,      0,  0;
        0,      0,      0,      0,      (rho_t*J*L)*(1/2)*2/m,  0;
        0,      0,      0,      0,      0,  0;];
        
    obj.elM_lumped=m/2*[Q1,Q2;Q2,Q4];
   if massless
        obj.elM_lumped(isnan(obj.elM_lumped))=0;
    end
    % Calculate Element Mass Matrices in Global Coordinates
    obj.elM_lumped_global=obj.T'*obj.elM_lumped*obj.T;

end
