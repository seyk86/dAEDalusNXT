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

function [obj] = f_postprocess(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %compute reaction forces for each node
    obj.reaction_forces=obj.K*obj.nodal_deflections-obj.Q;
    obj.node_loadings=obj.reaction_forces;
    obj.node_loadings_loc=zeros(length(obj.reaction_forces),1);
    %obj.node_loadings(1:6)=obj.reaction_forces(1:6);
    el_ndof=obj.el_ndof;
    k=1;
    
    %obj.node_loadings_loc(1:el_ndof)=obj.beamelement(1).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(1:el_ndof);
    if obj.isExternalFEM==0
    	for i=1:obj.nel-1
        	el_load=obj.beamelement(i).elKglobal*obj.nodal_deflections(k:k+2*el_ndof-1);
        	
        	el_loadnxt=obj.beamelement(i+1).elKglobal*obj.nodal_deflections(k+el_ndof:k+3*el_ndof-1);
        	obj.node_loadings(k:k+el_ndof-1)=el_load(1:el_ndof);
        	if(el_ndof==6)
            	el_load_loc=obj.beamelement(i).T*el_load;
            	el_load_loc2=obj.beamelement(i+1).T*el_loadnxt;
            	if i==1
            	    obj.node_loadings_loc(k:k+el_ndof-1)=el_load_loc(1:el_ndof);%+el_load_loc2(1:el_ndof)*0.5;
            	    obj.node_loadings_loc(k+el_ndof:k+2*el_ndof-1)=-el_load_loc(el_ndof+1:2*el_ndof)*0.5+el_load_loc2(1:el_ndof)*0.5;
            	elseif i==obj.nel-1
            	    obj.node_loadings_loc(k+2*el_ndof:k+3*el_ndof-1)=-el_load_loc2(el_ndof+1:2*el_ndof);%+el_load_loc2(1:el_ndof)*0.5;
            	    obj.node_loadings_loc(k+el_ndof:k+2*el_ndof-1)=-el_load_loc(el_ndof+1:2*el_ndof)*0.5+el_load_loc2(1:el_ndof)*0.5;%+el_load_loc2(1:el_ndof)*0.5;
            	else
            	    obj.node_loadings_loc(k+el_ndof:k+2*el_ndof-1)=-el_load_loc(el_ndof+1:2*el_ndof)*0.5+el_load_loc2(1:el_ndof)*0.5;%+el_load_loc2(1:el_ndof)*0.5;
            	end
        	end
        	k=k+el_ndof;
    	end
    end
   % obj.node_loadings(k:k+el_ndof-1)=-obj.reaction_forces(k:k+el_ndof-1);
  %  obj.node_loadings_loc(k:k+el_ndof-1)=-obj.beamelement(end).T(1:el_ndof,1:el_ndof)*obj.reaction_forces(k:k+el_ndof-1);
    %%TODO: check algorithm
    % now calculate nodal deflections at quarter chordline --> is input for
    % deflection for aerodynamic mesh
    for k=1:obj.nel+1
        if k==1
            T=obj.beamelement(k).T(1:6,1:6);
            obj.nodal_deflections_loc((k-1)*el_ndof+1:(k-1)*el_ndof+6)=T*obj.nodal_deflections((k-1)*el_ndof+1:(k-1)*el_ndof+6);
            rot=obj.nodal_deflections((k-1)*el_ndof+4:(k-1)*el_ndof+6);
            Theta_re=T(1:3,1:3)*[0 ;obj.epsilon(k)*0; 0]; 
        elseif k<=obj.nel
            T=obj.beamelement(k).T(1:6,1:6);
            obj.nodal_deflections_loc((k-1)*el_ndof+1:(k-1)*el_ndof+6)=T*obj.nodal_deflections((k-1)*el_ndof+1:(k-1)*el_ndof+6);
            rot=obj.nodal_deflections((k-1)*el_ndof+4:(k-1)*el_ndof+6);
            Theta_re=T(1:3,1:3)*[0 ;obj.epsilon(k)*0; 0];
        else
            rot=obj.nodal_deflections((k-1)*el_ndof+4:(k-1)*el_ndof+6);
            Theta_re=obj.beamelement(end).T(1:3,1:3)*[0; obj.epsilon(k); 0];
        end
        obj.nodal_deflections_c4((k-1)*el_ndof+4:(k-1)*el_ndof+6)=obj.nodal_deflections((k-1)*el_ndof+4:(k-1)*el_ndof+6);
        %distance in global x direction between node coordinate and c4 line
        %of wing (wingbox_c4-wingbox_coords)
        dist_x=abs(obj.dist_c4_sc(k));
        sgnx=sign(obj.dist_c4_sc(k));
        
        dTheta=rot(2);
        Theta=Theta_re(2)*0;
        dx2=sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
        %dy=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta))*0;
        dz2=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));

        dTheta=rot(3);
        Theta=Theta_re(3)*0;
        dx3=sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
        dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
        
        delta_twist=[dx2+dx3;dy3;dz2];
        obj.nodal_deflections_c4((k-1)*el_ndof+1:(k-1)*el_ndof+3)=obj.nodal_deflections((k-1)*el_ndof+1:(k-1)*el_ndof+3)+delta_twist;
        % obj.nodal_deflections_c4((k-1)*el_ndof+1:(k-1)*el_ndof+3)=obj.nodal_deflections((k-1)*el_ndof+1:(k-1)*el_ndof+3)+(M*diff'-diff');
    end
    %add nodal load later!
    %obj.node_loadings(k:k+el_ndof-1)=[0,0,0,0,0,0];
end

