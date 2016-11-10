%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ck_imag,ck_real,ck_magnit,ck_phase]=Theodorsen_parameters2(k)

% function to calculate the value of the lift deficiency factor using
% Hankel and Bessel functions

ss=k*exp(1i*pi*0.5);
Z=-1i*ss;
H_0=besselh(0,2,Z);
H_1=besselh(1,2,Z);
c_s=(H_1)/(H_1+1i*H_0);

ck_real=real(c_s);
ck_imag=imag(c_s);
ck_magnit=abs(c_s);
ck_phase=angle(c_s);

end
    
    
    
    
    
    
