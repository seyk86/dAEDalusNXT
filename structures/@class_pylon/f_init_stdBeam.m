%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_init_stdBeam :   file is part of nlFEM class_wing
%                    initalize entire structe with unified values for
%                    A, E, G, Ix, Iz, Ip
%                   for boxwing set: isboxwing=1;
% Author:           Klaus Seywald
%                   klaus.seywald@mytum.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % =================================================================
        %> @brief initialize beam with unified A E G Iy Iz Ip
        %>
        %> @param A vector containing crosssectional area of each element
        %> @param E vector containing young's modulus of each element
        %> @param G vector containing shear modulus of each element
        %> @param Iy vector containing second moment of area about y of each element
        %> @param Iz vector containing second moment of area about z of each element
        %> @param Ip vector containing torsional constant of each element
        %> @param m containing mass of each element
        %>
        %> @return initialized instance of the class_beam
        % =================================================================
function obj =f_init_stdBeam(obj,E,G,m)

    for i=1:1:obj.nel
        obj.beamelement(i).E=E;
        obj.beamelement(i).G=G;
        obj.beamelement(i).m=m;
    end
    
end
