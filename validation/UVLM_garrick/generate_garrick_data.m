%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%------------------------------------------------------------------------
% Name:
% FlightSim_Aeroelasticity_Bonin.m
%
% Description:
% ***
%
%-------------------------------------------------------------------------
%
% Author: J. D. Bonin                           Airbus Group Innovations
% Email : jannis.bonin@rwth-aachen.de
%
%-------------------------------------------------------------------------
%
% CHANGE_LOG:
% 13-Jan-14 , ***
%
%
% NOTES :
%  To set up phy_con working you need one k value very small!
% 
%%------------------------------------------------------------------------
%function [ nothing] =plot_fun(theta,rho)

 theta=0.5;

%% Data

 rho=1.225;

b=1;
a=-0.5;
c=0.4;


%% Aero

T1=-(1/3)*(2+c^2)*sqrt(1-c^2)+c*(cos(c))^-1;
T3=-(1/8)*(1-c^2)*(5*(c^2)+4)+(1/4)*c*(7+2*c^2)*sqrt(1-c^2)*((cos(c))^-1)-((c^2)+1/8)*((cos(c))^-1)^2;
T4=c*sqrt(1-c^2)-(cos(c))^-1;
T5=-(1-c^2)-(((cos(c))^-1)^2)+2*c*sqrt(1-c^2)*(cos(c))^-1;
T7=(1/8)*c*(7+2*c^2)*sqrt(1-c^2)-((c^2)+1/8)*((cos(c))^-1);
T8=-(1/3)*(1+2*c^2)*sqrt(1-c^2)+c*((cos(c))^-1);
T9=(1/2)*((1/3)*((1-c^2)^(3/2))+a*T4);
T10=sqrt(1-c^2)+(cos(c))^-1;
T11=(2-c)*sqrt(1-c^2)+(1-2*c)*(cos(c))^-1;
T12=sqrt(1-c^2)*(2+c)-(2*c+1)*(cos(c))^-1;
T13=-(1/2)*(T7+(c-a)*T1);
T15=T4+T10;
T16=T1-T8-(c-a)*T4+(1/2)*T11;
T17=-2*T9-T1+(a-(1/2))*T4;
T18=T5-T4*T10;
T19=(-1/2)*T4*T11;

M_nc=[-pi pi*a T1;pi*a -pi*(a^2+1/8) -T13;T1 -T13 T3/pi];
B_nc=[0 -pi T4;0 pi*(a-1/2) -T16;0 -T17 -T19/pi];
K_nc=[0 0 0;0 0 -T15;0 0 -T18/pi];

R_1=[-2*pi;2*pi*(a+1/2);-T12];
S_1=[0 1 T10/pi];
S_2=[1 0.5-a T11/(2*pi)];


r=k_UVLM(kk);

%     Z=(-j*k(kk));
    ss=r*exp(1i*pi*theta);
    Z=-1i*ss;
    H_0=besselh(0,2,Z);
    H_1=besselh(1,2,Z);
    c_s=(H_1)/(H_1+1i*H_0);
    
 %   c_s_approx=1-0.165/(1-0.0455/r*1i)-0.335/(1-0.3/r*1i);
  %  c_s=c_s_approx;
   % ck(kk,1)=c_s;
  %  ck(kk,2)=c_s_approx;
   % c_s=c_s_approx;
    Q(:,:,kk)=2*(b^2)*(M_nc*((ss)^2)+(B_nc+c_s*R_1*S_2)*(ss)+K_nc+c_s*R_1*S_1); 

    ndt=Uinf*wingaero.t_vec/(wingaero.reference.c_ref/2);
    
    Cd_heave=-2*pi*r^2*hh(kk)^2/(wingaero.reference.c_ref/2)^2*(imag(c_s)*cos(r*ndt)+real(c_s)*sin(r*ndt)).^2;
    heave=-hh(kk)*cos(r*ndt);
    alpha_eff=atan(r*2*norm(Uinf)/wingaero.reference.c_ref*hh(kk)/norm(Uinf))*sin(r*ndt);
    
    alpha_pitch=ah*sin(r*ndt);
    Y1=2*(real(c_s)-r*imag(c_s)*(1/2-aa));
    Y2=2*(imag(c_s)+r*real(c_s)*(1/2-aa))-r;
    Cd_suck=pi*ah^2/2*(Y1*sin(r*ndt)+Y2*cos(r*ndt)).^2;
    Cd_pitch=-Cd_suck+alpha_pitch*pi*ah.*(r*cos(r*ndt)+aa*r^2*sin(r*ndt)+...
        2*real(c_s)*(sin(r*ndt)+(1/2-aa)*r*cos(r*ndt))+...
        2*imag(c_s)*(cos(r*ndt)-(1/2-aa)*r*sin(r*ndt)));
  
