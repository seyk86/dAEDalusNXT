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

function obj=elM_6dof(obj)

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
    
    Q1=[156,    0,      0,      0,      0,  -22*L;
        0,      140,    0,      0,      0,  0;
        0,      0,      156,    22*L,   0,  0;
        0,      0,      22*L,   4*L^2,  0,  0;
        0,      0,      0,      0,      (rho_t*J*L)*(2/6)*420/m,  0;
        -22*L,   0,      0,      0,      0,  4*L^2;];
    
    Q2=[54,    0,      0,      0,      0,  13*L;
        0,     70,     0,      0,      0,  0;
        0,     0,      54,     -13*L,   0,  0;
        0,     0,      13*L,  -3*L^2, 0,  0;
        0,     0,      0,      0,      (rho_t*J*L)*(1/6)*420/m,  0;
        -13*L,  0,      0,      0,      0,  -3*L^2;];
    
    Q3=[54,     0,      0,      0,      0,  -13*L;
        0,      70,     0,      0,      0,  0;
        0,      0,      54,     13*L,   0,  0;
        0,      0,      -13*L,  -3*L^2, 0,  0;
        0,      0,      0,      0,      (rho_t*J*L)*(1/6)*420/m,  0;
        13*L,  0,      0,      0,      0,  -3*L^2;];
    
    Q4=[156,    0,      0,      0,       0,  22*L;
        0,      140,    0,      0,       0,  0;
        0,      0,      156,    -22*L,   0,  0;
        0,      0,      -22*L,   4*L^2,  0,  0;
        0,      0,      0,      0,      (rho_t*J*L)*(2/6)*420/m,  0;
        22*L,   0,      0,      0,      0,  4*L^2;];
    obj.elM=m/420*[Q1,Q2;Q3,Q4];
    if massless
        obj.elM(isnan(obj.elM))=0;
    end
   % obj.elM=m/420*diag([-1 1 1 -1 1 1 -1 1 1 -1 1 1])'*[Q1,Q2;Q3,Q4]*diag([-1 1 1 -1 1 1 -1 1 1 -1 1 1]);
    % Calculate Element Mass Matrices in Local Coordinates
%     obj.elM= m.*[156/420,0,54/420,22*L/420,0,-13*L/420,0,0,0,0,0,0;...      %ux1
%         0,2/6,0,0,0,0,0,1/6,0,0,0,0;...                                     %uy1
%         54/420,0,156/420,13*L/420,0,-22*L/420,0,0,0,0,0,0;...               %uz1
%         22*L/420,0,13*L/420,(4*L^2)/420,0,(-3*L^2)/420,0,0,0,0,0,0;...      %theta1
%         0,0,0,0,(rho*J*L/m)*(2/6),0,0,0,0,0,(rho*J*L/m)*(1/6),0;...         %phi1
%         -13*L/420,0,-22*L/420,(-3*L^2)/420,0,(4*L^2)/420,0,0,0,0,0,0;...    %theta1
%         0,0,0,0,0,0,156/420,0,54/420,22*L/420,0,-13*L/420;...               %ux2
%         0,1/6,0,0,0,0,0,2/6,0,0,0,0;...                                     %uy2        
%         0,0,0,0,0,0,54/420,0,156/420,13*L/420,0,-22*L/420;...               %uz2
%         0,0,0,0,0,0,22*L/420,0,13*L/420,(4*L^2)/420,0,(-3*L^2)/420;...      %theta2
%         0,0,0,0,(rho*J*L/m)*(1/6),0,0,0,0,0,(rho*J*L/m)*(2/6),0;...         %phi2
%         0,0,0,0,0,0, -13*L/420,0,-22*L/420,(-3*L^2)/420,0,(4*L^2)/420;      %theta2
%         ];
    % Calculate Element Mass Matrices in Global Coordinates
    obj.elMglobal=obj.T'*obj.elM*obj.T;

end
