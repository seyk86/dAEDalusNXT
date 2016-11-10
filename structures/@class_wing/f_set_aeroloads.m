%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_initaeroloads:  file is part of nlFEM class_wing 
%                   initializes FEM model with aerodynamic loads from 
%                   dAEDalus-VLM
%   Author:         Klaus Seywald
%                   klaus.seywald@mytum.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wingstr]= f_set_aeroloads(wingstr,wingaero)

        for i=1:1:length(wingstr.beamelement)      
            qloc=(wingaero.beam_forces_structmesh(:,i));
           % qloc=wingaero.c4_forces_structmesh(:,i);
           % dqloc=fac*cos(wingstr.beamelement(i).phi)*wingstr.beamelement(i).T(1:3,1:3)*(wingaero.c4_forces_structmesh(:,i+1)-wingaero.c4_forces_structmesh(:,i));
            mloc=(wingaero.beam_moments_structmesh(:,i));
           % mloc=(wingaero.c4_moments_structmesh(:,i));

            wingstr.beamelement(i).qx=qloc(1);
            wingstr.beamelement(i).qy=qloc(2);
            wingstr.beamelement(i).qz=qloc(3);
            
%             wingstr.beamelement(i).dqx=dqloc(1);
%             wingstr.beamelement(i).dqy=dqloc(2);
%             wingstr.beamelement(i).dqz=dqloc(3);
%             wingstr.beamelement(i).mx=mloc(1);
%             wingstr.beamelement(i).mz=mloc(3);
%             wingstr.beamelement(i).dmx=dmloc(1);
%             wingstr.beamelement(i).dmz=dmloc(3);
                      
            wingstr.beamelement(i).mt=mloc(2);
           % wingstr.beamelement(i).mt=mloc(2)+qloc(3)*(0.5*wingstr.dist_c4_sc(i+1)+0.5*wingstr.dist_c4_sc(i));%*cos(wingstr.nodal_deflections(5+6*(i-1))*0.5+wingstr.nodal_deflections(5+6*(i))*0.5);
%           wingstr.beamelement(i).dmt=dmloc(2)+dqloc(3)*(0.5*wingstr.dist_c4_sc(i+1)+0.5*wingstr.dist_c4_sc(i));
        end
        wingstr.update_Q=1;
end
