%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function fuselage = f_add_inertial_forces(fuselage,add_engines,add_gears)
  
%add engine forces at nearest node
  engine=fuselage.engine;
  gear=fuselage.gear;
  
  if add_engines==1     
        if fuselage.fuselagemountedengines == 1
            for ne=1:length(fuselage.engine)
                %find closest node:
                dist=zeros(fuselage.nel,1);
                for i=1:fuselage.nel
                    dist(i)=sqrt(sum((fuselage.node_coords(i,1)-engine(ne).cg_pos(1)').^2));
                end
                [Y,I] = min(dist);
                %add forces and moments
                
                Fx=engine(ne).delta_t.*engine(ne).thrust_vec(1)*engine(ne).thrust;
                Fy=engine(ne).delta_t.*engine(ne).thrust_vec(2)*engine(ne).thrust;
                Fz=engine(ne).delta_t.*engine(ne).thrust_vec(3)*engine(ne).thrust-fuselage.load_factor*fuselage.g*engine(ne).m;
                
                My=Fz*(engine(ne).cg_pos(2)-fuselage.node_coords(I,2))-Fx*(engine(ne).cg_pos(3)-fuselage.node_coords(I,3));
                    
                fuselage.Q(1+6*(I-1):1+6*(I-1)+5)=fuselage.Q(1+6*(I-1):1+6*(I-1)+5)+[Fx,Fy,Fz,0,My,0]';
            end
        end
  end
  
  %add gear forces at neares node
  if add_gears==1
        if fuselage.fuselagemountedgears == 1    
        %find closest node:
            for ne=1:length(fuselage.gear)
                dist=zeros(fuselage.nel,1);
                for i=1:fuselage.nel
                    dist(i)=sqrt(sum((fuselage.node_coords(i,1)-gear(ne).pos(1)').^2));
                end
                [Y,I] = min(dist);
                %add forces and moments
                if fuselage.landingimpact==1
                fuselage.Q(1+6*(I-1):1+6*(I-1)+5)=fuselage.Q(1+6*(I-1):1+6*(I-1)+5)+[gear(ne).Fx,gear(ne).Fy,gear(ne).Fz,0,0,0]';
                end
            end
        end
  end
  
end
