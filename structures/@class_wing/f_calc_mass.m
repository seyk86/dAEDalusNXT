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

function wing = f_calc_mass(wing,weights)
    m_total=0.0;
    m_wingbox_total=0.0;
    Awing=wing.Awing;
    V_wingbox=0;
    
    disp('            performing eigenmass calculation');
    %wingsection={wing.beamelement(:).crosssection};
    
    m_fuel_total=0;
    fuel_volume_total=0;
    for i=1:1:wing.nel
        %[h w c]=wing.beamelement(i).crosssection.get_dimensions();             
        le=wing.beamelement(i).le;
        phi=wing.beamelement(i).phi;
        
        c=wing.beamelement(i).crosssection.c;
        h=wing.beamelement(i).crosssection.h;
        w=wing.beamelement(i).crosssection.w;
        
        % calculate distributed mass from TSW (skin??)
        el_m_s=weights.WingSkinEstimate/Awing*le*c;%*cos(phi);
        % calculate mass from wing systems
%         el_m_sys=weights.WingSystemsEstimate/Awing*le*c;%*cos(phi);
        el_m_sys=weights.WingSystemsEstimate/wing.V_wingbox*h*w*le;
        % calculate wing box mass from 
        el_m_wingbox=wing.beamelement(i).crosssection.f_calc_dm*le;
       
        wing.beamelement(i).el_m_s=el_m_s;
        wing.beamelement(i).el_m_sys=el_m_sys;
        wing.beamelement(i).el_m_p=el_m_wingbox;
        
        wing.beamelement(i).el_m_fuel=wing.beamelement(i).el_fuel_vol*wing.fuel_density*wing.beamelement(i).crosssection.fueling_factor;
        
        el_m_sum=el_m_s+el_m_sys+el_m_wingbox;%+el_ribweight;
        wing.beamelement(i).m=(el_m_sum)/le;
        
        m_wingbox_total=m_wingbox_total+el_m_wingbox;
        m_total=m_total+el_m_sum;
        m_fuel_total=m_fuel_total+wing.beamelement(i).el_m_fuel;
        fuel_volume_total=fuel_volume_total+wing.beamelement(i).el_fuel_vol*wing.beamelement(i).is_fueled;
    end

    kc=1;

    wing.m_ribs=0;%kc*5*10^(-4)*2800*wing.Awing*((wing.beamelement(1).crosssection.w+wing.beamelement(end).crosssection.w)/2+1);
    wing.m_total=m_total+wing.m_ribs+m_fuel_total;
    wing.m_wingbox_total=m_wingbox_total;
    wing.m_fuel_total=m_fuel_total;
    wing.fuel_volume_total=fuel_volume_total;

    wing.update_Q=1;
end

