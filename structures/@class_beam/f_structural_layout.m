%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_structural_layout:  class file is part of nlFEM, class class_wing
%                       performs structural layout a'la Rafik 
%     Author:           Klaus Seywald
%                       klaus.seywald@mytum.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % ================================================================
         %> @brief perform element wise structural layout of crosssection
         %>
         % ================================================================

function obj = f_structural_layout(obj,overwrite)
    disp('            performing structural layout calculation');
    Mbx=abs(obj.node_loadings_loc(4:6:end));%+abs(wing.node_loadings(6:6:end)); %bending moment about x (extend.. 
    Mbz=abs(obj.node_loadings_loc(6:6:end));
    Mt=abs(obj.node_loadings_loc(5:6:end)); %torsional moment
    Qz=abs(obj.node_loadings_loc(3:6:end));
    Qx=abs(obj.node_loadings_loc(1:6:end));
   
    
%     wingsection={obj.beamelement(:).wingsection};
%     beambeam={obj.beamelement(:)};

    for i=1:1:obj.nel
        % self design crosssections
        
        
        if obj.is_sym
            %Simon: Mittelung der Lasten rechts und Links
%             MBX=Mbx(i)*0.25+Mbx(i+1)*0.25+Mbx(end-i+1)*0.25+Mbx(end-i)*0.25;
%             MBZ=Mbz(i)*0.25+Mbz(i+1)*0.25+Mbz(end-i+1)*0.25+Mbz(end-i)*0.25;
%             MT=Mt(i)*0.25+Mt(i+1)*0.25+Mt(end-i+1)*0.25+Mt(end-i)*0.25;
%             QX=Qx(i)*0.25+Qx(i+1)*0.25+Qx(end-i+1)*0.25+Qx(end-i)*0.25;
%             QZ=Qz(i)*0.25+Qz(i+1)*0.25+Qz(end-i+1)*0.25+Qz(end-i)*0.25;

            %Simon: besser: jeweils h�here Lasten, entweder re oder li
            MBX=max(Mbx(i)*0.5+Mbx(i+1)*0.5,Mbx(end-i+1)*0.5+Mbx(end-i)*0.5);
            MBZ=max(Mbz(i)*0.5+Mbz(i+1)*0.5,Mbz(end-i+1)*0.5+Mbz(end-i)*0.5);
            MT=max(Mt(i)*0.5+Mt(i+1)*0.5,Mt(end-i+1)*0.5+Mt(end-i)*0.5);
            QX=max(Qx(i)*0.5+Qx(i+1)*0.5,Qx(end-i+1)*0.5+Qx(end-i)*0.5);
            QZ=max(Qz(i)*0.5+Qz(i+1)*0.5,Qz(end-i+1)*0.5+Qz(end-i)*0.5);
        else
            MBX=Mbx(i)*0.5+Mbx(i+1)*0.5;
            MBZ=Mbz(i)*0.5+Mbz(i+1)*0.5;
            MT=Mt(i)*0.5+Mt(i+1)*0.5;
            QX=Qx(i)*0.5+Qx(i+1)*0.5;
            QZ=Qz(i)*0.5+Qz(i+1)*0.5;
        end
        if isa(obj,'class_fuselage')
            if ~obj.is_rigid
                obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.f_self_design_crosssection(MBX,MBZ,MT,QX,QZ,obj.loadcase_index);
                obj.beamelement(i)=obj.beamelement(i).f_calcCrossProp();
            else
                obj.beamelement(i).Iz=1E10;
                obj.beamelement(i).Ix=1E10;
                obj.beamelement(i).Ip=1E10;
                obj.beamelement(i).J=1E10;
                obj.beamelement(i).A=1E10;
                obj.beamelement(i).A_enclosed=1E10;
            end
        elseif isa(obj,'class_pylon')
            obj.beamelement(i).A=0.01;
            obj.beamelement(i).Ix=5e-6;
            obj.beamelement(i).Iz=5e-6;
            obj.beamelement(i).Ip=10e-6;
            obj.beamelement(i).J=10e-6;
            obj.beamelement(i).m=0;
            obj.beamelement(i).A_enclosed=0.01;   
        else
            obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.f_self_design_crosssection(MBX,MBZ,MT,QX,QZ,obj.loadcase_index,overwrite);
            obj.beamelement(i)=obj.beamelement(i).f_calcCrossProp();
        end
        %wingsection{i}.f_self_design_crosssection(Mb(i),Mt(i),Q(i)); 
        % compute fuel volume per element
%         if(beambeam{1}(i).is_fueled)
%             beambeam{1}(i).el_fuel_vol=wingsection{i}.A_fuel*beambeam{1}(i).le;
%         end
    end
    obj.update_K=1;
    obj.update_M=1;
%     for i=1:1:obj.nel
%         beambeam{1}(i).wingsection=wingsection{i};
%         beambeam{1}(i)=beambeam{1}(i).f_calcCrossProp();
%         obj.beamelement(i)=beambeam{1}(i);
%     end
end

