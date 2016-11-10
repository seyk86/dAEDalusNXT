%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_beamelement.m
%> @brief File containing the class for a finite element beam collection
%>        This file is part of dAEDalus structures, Copyright (C) 2011,
%>        Klaus Seywald
%>   dAEDalus is published under the terms of the GNU General Public
%>   License by the Free Software Foundation;
%>   either version 2 or any later version.
%>
%>   You should have received a copy of the GNU General Public
%>   License along with dAEDalus; see the file GNU GENERAL 
%>   PUBLIC LICENSE.TXT.  If not, write to the Free Software 
%>   Foundation, 59 Temple Place -Suite 330, Boston, MA
%>   02111-1307, USA.
%>
% ======================================================================
%> @brief class for a finite beam element
%> This class is the base class for a finite element beam, the follwing
%> coordinate system is used
%> Global Coordinate Definition (according to Airbus Wing Geometry Definitions ):
%>   y= along starboard wing
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
%> Local Coordinate Definition (Cartesian Coordinate System )
%>   y= beam axis
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
% ======================================================================

classdef class_beamelement
    properties
        %% general element information 
        %> element length                                [m]
        le=0.0;    
        %> element rotation around z axis (sweep)        [rad]
        phi=0.0; 
        %> element rotation around x axis (dihedral)     [rad] 
        nu=0.0;     
        %> element rotation around y axis (twist)        [rad]
        epsilon=0.0;
        %> Youngs modulus                                [N/m2]
        E=21E+10;
        %> Shear modulus                                 [N/m2]
        G=7.93E+10; 
        %> element cross-section area                    [m2]
        A=1;
        %> Moment of inertia about x-axis                [m4]
        Ix=1; 
        %> Moment of inertia about z-axis                [m4]
        Iz=1; 
        %> element polar moment of inertia               [m4]
        Ip=1;   
        %> zx Moment of inertia about
        Izx=0; 
        %> excentrictiy
        e=0;
           
        J=1; %torsional moment 
        
        %% Forces and Moments (in element local coordinates)
        %> element transverse pressure load              [N/m]
        qz=0.0; 
        %> 
        dqz=0.0;
        
        %> element transverse pressure load              [N/m]
        qy=0.0;
        
        dqy=0.0;
        %> element transverse pressure load              [N/m]
        qx=0.0; 
        dqx=0.0;
        %> element distributed torsion load              [Nm/m]
        mt=0.0; 
        dmt=0.0;
        %> element distributed moment load               [Nm/m]
        mx=0.0;
        dmx=0.0;
        %> element distributed moment load               [Nm/m]
        mz=0.0; 
        dmz=0.0;
        %> element distributed load due to eigenmass     [Nm/m]
        qm=0.0;     
        
        %% accelerations defined in global coordinate system
        %> element distributed load due to x acceleration     [Nm/m]
        ax=0.0; 
        %> element distributed load due to y acceleration       [Nm/m]
        ay=0.0; 
        %> element distributed load due to z acceleration      [Nm/m]
        az=0.0; 
        
        
        %> element distrubuted load due to fuel          [Nm/m]
        qf=0.0;    

        %% additions for non linear element
        elp=0.0;
        elpglobal=0.0;              % element internal force vector
        nodal_deflections_loc;      
        
        %% Masses
        %>  element structural mass per unit length       [kg/m]
        m=1;          
        
        %> element fueled y/n                     []
        is_fueled=0;
        %> fuel volume per element                [m3]
        el_fuel_vol=0.0 
        %> enclosed area
        A_enclosed=0;         

        sigma_b=0.0;
        %> element degrees of freedom
        eldof;          
        %> element stiffness matrix in local coordinates
        elK;            
        %> element stiffness matrix in global system coordinates
        elKglobal;    
        %> element mass matrix
        elM;
        %> element mass matrix in global coordinates
        elMglobal;
        %> element mass matrix
        elM_lumped;
        %> element mass matrix in global coordinates
        elM_lumped_global;
        %> element load vector in local coordinates
        elq;           
        %> element load vector in global coordinates
        elqglobal;
        %> element mass
        el_m_s=0.0;
        %> element system mass
        el_m_sys=0.0;   
        el_m_p=0.0;
        
        %> element fuel mass
        el_m_fuel=0.0;

        %> transformation matrix from global system to local system
        T;
        
        Nm=0;
        Vm;
        Mm;
        kappa;
        
        crosssection;
        %> reference to beam to which the element belongs to
        ref_parent_beam;  
        

    end
    
    methods
        function obj =class_beamelement(crosssection,parent_ref)
            obj.crosssection=crosssection;
            obj.ref_parent_beam=parent_ref;
        end
       
       
        %compute element stiffness matrix for 3 DOF element(outdated, remove in
        %release)
        obj=lin_elK_3dof(obj);
        
        %compute element stiffness matrix for 6 DOF element 
        obj=lin_elK_6dof(obj);
        
        % Compute element mass matrix for 6DOF element
         obj=elM_6dof(obj);
         
        obj=elM_lumped_6dof(obj);
        
        %compute element load vector for 3 DOF element (outdated, remove in
        %release)
        obj=elq_3dof(obj);
        %compute element load vector for 6 DOF element
        obj=elq_6dof(obj,add_eigenmass,add_fuelmass);
        
        obj = f_calcCrossProp(obj)
        
        
        %set geometry of element
        function obj=setElementGeometry(obj,le,phi,nu,twist)
           obj.le=le;
           obj.phi=phi;
           obj.nu=nu; 
           obj.epsilon=twist;
           obj.T=obj.f_rotM_6dof(real(obj.nu),real(obj.epsilon),real(obj.phi));
           %obj.T=obj.f_rotM_3dof(0,0,-obj.phi);
        end
    end    
end





