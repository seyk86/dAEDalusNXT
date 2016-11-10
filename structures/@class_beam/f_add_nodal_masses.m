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

function obj = f_add_nodal_masses(obj)

  %for every nodal mass which has nonzero weight do
  conm=find(obj.nodal_masses(:,1));
  for i=conm'
    %calculate distance to connected node
    %dist=abs(obj.node_coords(i,:)-obj.nodal_masses(i,2:4));
    
    %better way: in nodal_masses the position is already defined as offset
    %from the respective node
    distabs=abs(obj.nodal_masses(i,2:4));
    dist=obj.nodal_masses(i,2:4);
    %create additional mass matrix
    m=obj.nodal_masses(i,1);
    mass_matrix1=[m 0 0; 0 m 0; 0 0 m];
    
    mass_matrix2= [  0          m*dist(3)	-m*dist(2)	;
                     -m*dist(3) 0         	m*dist(1)	;
                     m*dist(2)  -m*dist(1)	0           ];
                 
    mass_matrix3= -mass_matrix2;
    
    mass_matrix4=[   obj.nodal_masses(i,5)+(distabs(2)^2+distabs(3)^2)*m          -dist(1)*dist(2)*m                              -dist(1)*dist(3)*m;
                            -dist(1)*dist(2)*m                              obj.nodal_masses(i,6)+(distabs(1)^2+distabs(3)^2)*m   -dist(2)*dist(3)*m;
                            -dist(1)*dist(3)*m                              -dist(2)*dist(3)*m                              obj.nodal_masses(i,7)+(distabs(1)^2+distabs(2)^2)*m];
                   
   mass_matrix5=[   obj.nodal_masses(i,5)+(distabs(2)^2+distabs(3)^2)*m          0                             0;
                            0                             obj.nodal_masses(i,6)+(distabs(1)^2+distabs(3)^2)*m   0;
                           0                              0                             obj.nodal_masses(i,7)+(distabs(1)^2+distabs(2)^2)*m];
    mass_matrix=[mass_matrix1 mass_matrix2; mass_matrix3 mass_matrix4];
%     lumped_mass_matrix=[mass_matrix1 zeros(3,3); zeros(3,3) mass_matrix5];%not sure if the lumped mass matrix is built up correctly!!!!
         
                    
%       lumped_mass_matrix=   [m,     0,      0,      0,                                                  m*dist(3),                                      -m*dist(2);
%                             0,      m,      0,      -m*dist(3),                                         0,                                              m*dist(1);
%                             0,      0,      m,      m*dist(2),                                          -m*dist(1),                                     0;
%                             0,   -m*dist(3),m*dist(2),      obj.nodal_masses(i,5)+(dist(2)^2+dist(3)^2)*m,      -dist(1)*dist(2)*m,                             -dist(1)*dist(3)*m;
%                            m*dist(3),0, -m*dist(1),       -dist(1)*dist(2)*m,                                 obj.nodal_masses(i,6)+(dist(1)^2+dist(3)^2)*m,  -dist(2)*dist(3)*m;
%                             -m*dist(2),m*dist(1),      0,      -dist(1)*dist(3)*m,                                 -dist(2)*dist(3)*m,                             obj.nodal_masses(i,7)+(dist(1)^2+dist(2)^2)*m;];
%                 
    lumped_mass_matrix=mass_matrix;
    %add mass matrix to beam.M and beam.M_lumped
    obj.M(1+6*(i-1):6*i,1+6*(i-1):6*i)=obj.M(1+6*(i-1):6*i,1+6*(i-1):6*i)+mass_matrix;
                
    obj.M_lumped(1+6*(i-1):6*i,1+6*(i-1):6*i)=obj.M_lumped(1+6*(i-1):6*i,1+6*(i-1):6*i)+lumped_mass_matrix;
                
  end
    



% 
%   %add engine forces at nearest node
%   if add_engines==1     
%         if beam.wingmountedengines == 1
%             for ne=1:length(beam.engine)
%                 %find closest node:
%                 dist=zeros(beam.nel,1);
%                 for i=1:beam.nel
%                     dist(i)=sqrt(sum((beam.node_coords(i,2)-engine(ne).cg_pos(2)').^2));
%                 end
%                 [Y,I] = min(dist);
%                 %add forces and moments
% 
%                 engine_lumped_mass=[1,      0,      0,      0,      0,  0;
%                                     0,      1,      0,      0,      0,  0;
%                                     0,      0,      1,      0,      0,  0;
%                                     0,      0,      0,      0,      0,  0;
%                                     0,      0,      0,      0,      norm(-engine(ne).cg_pos(1)+beam.node_coords(I,1))^2*2,  0;
%                                     0,      0,      0,      0,      0,  0;]*engine(ne).m/2;
%                 
%                 
%                 beam.M(1+6*(I-1):1+6*(I-1)+5,1+6*(I-1):1+6*(I-1)+5)=beam.M(1+6*(I-1):1+6*(I-1)+5,1+6*(I-1):1+6*(I-1)+5)+engine_lumped_mass;
%                 
%                 beam.M_lumped(1+6*(I-1):1+6*(I-1)+5,1+6*(I-1):1+6*(I-1)+5)=beam.M_lumped(1+6*(I-1):1+6*(I-1)+5,1+6*(I-1):1+6*(I-1)+5)+engine_lumped_mass;
%                 
%             end
%         end
%   end  

end

