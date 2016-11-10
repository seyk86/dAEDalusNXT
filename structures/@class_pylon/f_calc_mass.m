%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% ================================================================
%> @brief initialize beam geometry manually
%>
%> @param weights class of class_weights
%>
%> @return object of class_beam with calculated mass
% ================================================================

function pylon = f_calc_mass(pylon,weights)
    
    disp('            performing eigenmass calculation');
   m_total=0;
    for i=1:1:pylon.nel
         le=pylon.beamelement(i).le;
        el_m_sum=pylon.beamelement(i).el_m_sys+pylon.beamelement(i).el_m_s+pylon.beamelement(i).el_m_p;
        pylon.beamelement(i).m=(el_m_sum)/le;

        m_total=m_total+el_m_sum;    
    end

    pylon.m_total=m_total;
    pylon.update_Q=1;   
end
