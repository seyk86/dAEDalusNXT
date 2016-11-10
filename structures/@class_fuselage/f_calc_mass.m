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

function fuselage = f_calc_mass(fuselage,weights)
    
    m_total=0.0;
    m_shell_total=0.0;
    Vfuselage=fuselage.Vfuselage;
    
    % Raymer coefficients accounting for cut-outs in the shell for doors,
    % windows and landing gears
%     k_doors=1.25;
%     k_lg=1.12;
    k_doors=1;
    k_lg=1;
    
    for i=1:length(weights.FuselageNonStructuralEstimate)
       if (weights.FuselageNonStructuralEstimate(i)<weights.NumberPAX(i)*(80+30+10)) %hardcoded get from class weights later
           weights.FuselageNonStructuralEstimate=weights.NumberPAX(i)*(80+30+10); %passenger + baggage =80+30+10
              disp(['           warning nonstructural mass estimate lower than passenger weight for' num2str(weights.NumberPAX(i)) ' Passengers with 80kg and 30kg luggage and 10kgs for seat']);
              disp('           estimate overwritten!!');
       end    
    end
    
    fuselage.n_PAX=weights.NumberPAX;
    
    disp('            performing eigenmass calculation');

    for i=1:1:fuselage.nel
        
        [r]=fuselage.beamelement(i).crosssection.get_dimensions();            
        le=fuselage.beamelement(i).le;
        Vel=pi*le*r^2;
        
        % calculates distributed mass from non structural mass
        if i<(fuselage.nel)
            el_m_nonstruct=(weights.FuselageNonStructuralEstimate-(1-fuselage.beamelement(i).crosssection.filling_ratio)*80*weights.NumberPAX)/Vfuselage*Vel;
        else
            el_m_nonstruct=0;
        end
        % calculates mass from fuselage systems
        el_m_sys=weights.FuselageSystemsEstimate/Vfuselage*Vel;
        % calculates distributed fuselage stiffened shell mass 
        el_m_shell=k_doors*k_lg*fuselage.beamelement(i).crosssection.f_calc_dm*le;
%         el_m_shell=fuselage.beamelement(i).crosssection.f_calc_dm*le;
        
        fuselage.beamelement(i).el_m_s=el_m_nonstruct;
        fuselage.beamelement(i).el_m_sys=el_m_sys;
        fuselage.beamelement(i).el_m_p=el_m_shell;
        
        el_m_sum=el_m_nonstruct+el_m_sys+el_m_shell;
        fuselage.beamelement(i).m=(el_m_sum)/le;
        
        m_shell_total=m_shell_total+el_m_shell;
        m_total=m_total+el_m_sum;
        
    end

    fuselage.m_total=m_total;

    fuselage.m_shell_total=m_shell_total;
    
    fuselage.update_Q=1;
    
end
