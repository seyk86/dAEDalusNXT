%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ac,aircraft_structure,wingaero] = critical_case_layout(aircraft,aircraft_structure,state,weights,aeroelastic_solver_settings, varargin)


if nargin==6
    overwrite=varargin{1};
else
    overwrite=1;
end

    disp(['performing critical case layout for: ' aircraft_structure.f_get_name()]);
    
    
%     for i=1:length(wingstructure.beamelement)
%             wingstructure.beamelement(i).wingsection.t_sp_fr_lc_idx=1;
%             wingstructure.beamelement(i).wingsection.t_sp_re_lc_idx=1;
%             wingstructure.beamelement(i).wingsection.t_sk_up_lc_idx=1;
%             wingstructure.beamelement(i).crosssection.t_sk_lo_lc_idx=1;
%     end


for i=1:length(state)
   disp(['   critical case ' num2str(i) ' out of ' num2str(length(state))]);
    wingstructure_prv=aircraft_structure.f_copy; %deep copy wingstructure
    [ac,aircraft_structure,wingaero]=structural_sizing_loop(aircraft,aircraft_structure,state(i),weights,aeroelastic_solver_settings,overwrite);
    
    for nb=1:length(aircraft_structure.beam)
        if (isa(aircraft_structure.beam(nb),'class_wing') && (aircraft_structure.beam(nb).is_sym==1))
            for ne=1:floor(length(aircraft_structure.beam(nb).beamelement)/2)
                aircraft_structure.beam(nb).beamelement(ne).Ix=max([aircraft_structure.beam(nb).beamelement(ne).Ix aircraft_structure.beam(nb).beamelement(end-ne+1).Ix]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).Ix=aircraft_structure.beam(nb).beamelement(ne).Ix;
                
                aircraft_structure.beam(nb).beamelement(ne).Ip=max([aircraft_structure.beam(nb).beamelement(end-ne+1).Ip aircraft_structure.beam(nb).beamelement(ne).Ip]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).Ip=aircraft_structure.beam(nb).beamelement(ne).Ip;
                
                aircraft_structure.beam(nb).beamelement(ne).Iz=max([aircraft_structure.beam(nb).beamelement(ne).Iz aircraft_structure.beam(nb).beamelement(end-ne+1).Iz]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).Iz=aircraft_structure.beam(nb).beamelement(ne).Iz;
                
                aircraft_structure.beam(nb).beamelement(ne).Izx=max([aircraft_structure.beam(nb).beamelement(ne).Izx aircraft_structure.beam(nb).beamelement(end-ne+1).Izx]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).Izx=aircraft_structure.beam(nb).beamelement(ne).Izx;

                aircraft_structure.beam(nb).beamelement(ne).J=max([aircraft_structure.beam(nb).beamelement(ne).J aircraft_structure.beam(nb).beamelement(end-ne+1).J]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).J=aircraft_structure.beam(nb).beamelement(ne).J;
                
                aircraft_structure.beam(nb).beamelement(ne).A=max([aircraft_structure.beam(nb).beamelement(ne).A aircraft_structure.beam(nb).beamelement(end-ne+1).A]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).A=aircraft_structure.beam(nb).beamelement(ne).A;
                
                aircraft_structure.beam(nb).beamelement(ne).A_enclosed=max([aircraft_structure.beam(nb).beamelement(ne).A_enclosed aircraft_structure.beam(nb).beamelement(end-ne+1).A_enclosed]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).A_enclosed=aircraft_structure.beam(nb).beamelement(ne).A_enclosed;
                
                aircraft_structure.beam(nb).beamelement(ne).m=max([aircraft_structure.beam(nb).beamelement(ne).m aircraft_structure.beam(nb).beamelement(end-ne+1).m]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).m=aircraft_structure.beam(nb).beamelement(ne).m;

                aircraft_structure.beam(nb).beamelement(ne).el_m_s=max([aircraft_structure.beam(nb).beamelement(ne).el_m_s aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_s]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_s=aircraft_structure.beam(nb).beamelement(ne).el_m_s;

                aircraft_structure.beam(nb).beamelement(ne).el_m_sys=max([aircraft_structure.beam(nb).beamelement(ne).el_m_sys aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_sys]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_sys=aircraft_structure.beam(nb).beamelement(ne).el_m_sys;
                 
                aircraft_structure.beam(nb).beamelement(ne).el_m_p=max([aircraft_structure.beam(nb).beamelement(ne).el_m_p aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_p]);
                aircraft_structure.beam(nb).beamelement(end-ne+1).el_m_p=aircraft_structure.beam(nb).beamelement(ne).el_m_p;
                [~,I]=max([aircraft_structure.beam(nb).beamelement(end-ne+1).Iz aircraft_structure.beam(nb).beamelement(ne).Iz]);
                if I==1
                    aircraft_structure.beam(nb).beamelement(ne).crosssection=aircraft_structure.beam(nb).beamelement(end-ne+1).crosssection;
                else
                    aircraft_structure.beam(nb).beamelement(end-ne+1).crosssection=aircraft_structure.beam(nb).beamelement(ne).crosssection;
                end
            end
        end
        
       if i>=2
           aircraft_structure=combine_structures(wingstructure_prv,aircraft_structure,state(i).loadcase_index);
           aircraft_structure=aircraft_structure.f_calc_mass(weights);
       end
    end
end

