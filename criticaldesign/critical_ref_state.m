%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ref_state] =critical_ref_state(aircraft,MCruise,StructuralDesignAltitude)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%based on: http://adg.stanford.edu/aa241/structures/placard.html


%aircraft_state=class_aircraft_state(aircraft,[0.083333333*aircraft.reference.c_ref 0 0]);

aircraft_state=class_aircraft_state(aircraft,aircraft.reference.p_ref);
[rho_air,a,T,P,mu]=stdatmo(StructuralDesignAltitude);

V=[MCruise*a 0 0];

ref_state=class_flight_state(StructuralDesignAltitude,V,0,0,aircraft_state);

Mc=1.06*MCruise; %maximum operating mach number, kroo
MD=1.07*Mc; %design dive mach number, kroo

if(ref_state.h<11000)
    VDTAS=MD*a;%*sqrt((288.15-0.0065*ref_state.Altitude)/288.15); % design dive TAS
    VCTAS=MCruise*a;%*sqrt((288.15-0.0065*ref_state.Altitude)/288.15); % design dive TAS
else
    VDTAS=MD*a;%*sqrt((288.15-0.0065*11000)/288.15); % design dive TAS
    VCTAS=MCruise*a;%*sqrt((288.15-0.0065*11000)/288.15); % design dive TAS
end
VCEAS=VCTAS*sqrt(rho_air/1.225);  % design dive EAS
VDEAS=VDTAS*sqrt(rho_air/1.225);  % design dive EAS
ref_state=ref_state.f_set_V_A(norm(V));

ref_state.VD=VDTAS;
ref_state.VC=VCTAS;
ref_state.VC_EAS=VCEAS;
ref_state.VD_EAS=VDEAS;

end

