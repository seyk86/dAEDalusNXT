%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%

% THIS FUNCTION IS THE MODIFIED VERSION OF THE BONIN IMPLEMENTATION!
%Bonin: generate_garrick_data

r=k_UVLM(kk);

%% Lift deficiency Factor
[ck_imag,ck_real,~,~]=Theodorsen_parameters2(r);
c_s=ck_real+1i*ck_imag;

%% Omega*time
ndt=Uinf*wingaero.t_vec/(wingaero.reference.c_ref/2);

%% Induced drag for Heave
Cd_heave=-2*pi*r^2*hh(ii)^2/(wingaero.reference.c_ref/2)^2*(imag(c_s)*cos(r*ndt)+real(c_s)*sin(r*ndt)).^2;
heave=-hh(ii)*cos(r*ndt);
alpha_eff=atan(r*2*norm(Uinf)/wingaero.reference.c_ref*hh(ii)/norm(Uinf))*sin(r*ndt);

%% Induced drag for Pitch
alpha_pitch=ah*sin(r*ndt);
Y1=2*(real(c_s)-r*imag(c_s)*(1/2-a));
Y2=2*(imag(c_s)+r*real(c_s)*(1/2-a))-r;
Cd_suck=pi*ah^2/2*(Y1*sin(r*ndt)+Y2*cos(r*ndt)).^2;
Cd_pitch=-Cd_suck+alpha_pitch*pi*ah.*(r*cos(r*ndt)+a*r^2*sin(r*ndt)+...
    2*real(c_s)*(sin(r*ndt)+(1/2-a)*r*cos(r*ndt))+...
    2*imag(c_s)*(cos(r*ndt)-(1/2-a)*r*sin(r*ndt)));
