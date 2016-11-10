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
% elq_6dof:         is a method of class beamelement
%                   assemble element load vector for 6DOF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = elq_6dof(obj,add_eigenmass,add_fuelmass)

if isa(obj.crosssection,'class_crosssection_C_shapeV2');
%     mt=obj.mt+obj.crosssection.e*obj.qz;
    mt=obj.mt+(obj.crosssection.e+obj.crosssection.w/2)*obj.qz;  %center of cross section
else
    mt=obj.mt;
end
    mx=obj.mx*0;
    mz=obj.mz*0;
    le=obj.le;
    
    %perform load rotation (follower load) for aerodynamic forces
    %(nonlinear solution only)
%% Uncommented for performance improvement: REQUIRED FOR NONLINEAR COMPUTATION  
%     if sum(obj.nodal_deflections_loc)~=0.0
%         T=obj.f_rotVec(-obj.nodal_deflections_loc(10),0,0);  
%         
%         %obj.nodal_deflections_loc(10)*180/pi;
%         qaerrot=T*[qx;qy;qz];
%         
%         qx=qaerrot(1);
%         qy=qaerrot(2);
%         qz=qaerrot(3);
%     end

    %% adding distributed loads due to eigenmass and fuelweight
%% Uncommented for performance improvement: REQUIRED FOR NONLINEAR COMPUTATION      
    %if add_eigenmass
%         m=obj.m;
%         qm=-m*g;                        % distributed load due to eigenmass
%         obj.qm=qm;
%         T=obj.f_rotVec(0,gamma,0);      % transformations matrix from NED to aerodynamic system
%         qmrot=T*[0;0;qm];               % distributed loading due to eigenmass in x,y,z coordinates
%         qx=qx+qmrot(1);
%         qy=qy+qmrot(2);
%         qz=qz+qmrot(3);
    %end
    
    if obj.is_fueled
        rho_fuel=obj.ref_parent_beam.fuel_density;
        obj.el_fuel_vol=obj.A_enclosed*obj.le;
        obj.el_m_fuel=obj.el_fuel_vol*rho_fuel;
        m_tot=obj.m+obj.el_fuel_vol*rho_fuel/obj.le;
    else
        m_tot=obj.m;
    end
    
    q_inertia=obj.T(1:3,1:3)*[obj.ax*m_tot;obj.ay*m_tot;obj.az*m_tot];
    qloc=[obj.qx;obj.qy;obj.qz;];
    qx=qloc(1);
    qy=qloc(2);
    qz=qloc(3);
    
    qx=qx+q_inertia(1);
    qy=qy+q_inertia(2);
    qz=qz+q_inertia(3);
    
%     if add_fuelmass 
%         if obj.is_fueled
%             rho_fuel=obj.ref_parent_beam.fuel_density;
%             obj.el_fuel_vol=obj.A_enclosed*obj.le;
%             qf=-obj.A_enclosed*obj.le*rho_fuel/obj.le*g;                     % distributed load due to eigenmass
%             obj.qf=qf;
%             T=obj.f_rotVec(0,gamma,0);  % transformations matrix from NED to aerodynamic system
%             qfrot=T*[0;0;qf];               % distributed loading due to eigenmass in x,y,z coordinates
%            
%             qx=qx+qfrot(1);
%             qy=qy+qfrot(2);
%             qz=qz+qfrot(3);
%         end
%     end
    
    % add eigenmass to external distributed loads
    fract=0.03333333333333333333333333333;
    %% assuming linear shape functions and piecewise constant distributed loading! (possible extension)
    obj.elq=([qx*le/2+0.15*obj.dqx*le^2,qy*le/2+0.15*obj.dqy*le^2,qz*le/2+0.15*obj.dqz*le^2,...   
        mx*le/2+0.15*obj.dmx*le^2+qz*le^2/12+fract*obj.dqz*le^3,mt*le/2+0.15*obj.dmt*le^2,mz*le/2+0.15*obj.dmz*le^2-qx*le^2/12-fract*obj.dqx*le^3,...    
        (qx+obj.dqx)*le/2-0.15*obj.dqx*le^2,(qy+obj.dqy)*le/2-0.15*obj.dqy*le^2,(qz+obj.dqz)*le/2-0.15*obj.dqz*le^2,...     
        (mx+obj.dmx)*le/2-0.15*obj.dmx*le^2-(qz+obj.dqz)*le^2/12-fract*obj.dqz*le^3,(mt+obj.dmt)*le/2-0.15*obj.dmt*le^2,(mz+obj.dmz)*le/2+0.15*obj.dmz*le^2+(qx+obj.dqx)*le^2/12+fract*obj.dqx*le^3])';
%     obj.elq=[qx*le/2,qy*le/2,qz*le/2,...
%         mx*le/2+qz*le^2/12,mt*le/2,mz*le/2+qx*le^2/12,...
%         qx*le/2,qy*le/2,qz*le/2,...
% mx*le/2-qz*le^2/12,(mt)*le/2,mz*le/2-qx*le^2/12]';
    % T=obj.f_rotM_6dof(obj.nu,0,-obj.phi);
    
    %obj.elqglobal=obj.T^(-1)*obj.elq;
    obj.elqglobal=obj.T'*obj.elq;
    %obj.elqglobal=obj.elq;
end
