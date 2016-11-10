%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Q]=generate_ts_data_mod_function(a,b,c,k_UVLM)

% THIS FUNCTION IS THE MODIFIED VERSION OF THE BONIN IMPLEMENTATION!
% MODIFICATION: script to function adaptation
% Bonin: generate_ts_data
% REFERENCE ftp://161.24.15.247/Gil/AE-249/Computational%20Aids%20in%20Aeroservoelastic%20Analysis%20Using%20MATLAB.pdf

% Q matrix -  Lift(row1), Pitching Moment(row2), column1(Heave Motion),
% column2(Pitch Motion)

%% Constants
 theta=0.5; % constant to calculate the lift deficiency factor

%% Geometric Theodorsen constants

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

%% Mass, damping and Stiffness matrices
M_nc=[-pi pi*a T1;pi*a -pi*(a^2+1/8) -T13;T1 -T13 T3/pi]; %% check for correctness see  control of spacecraft and aircraft book
B_nc=[0 -pi T4;0 pi*(a-1/2) -T16;0 -T17 -T19/pi];
K_nc=[0 0 0;0 0 -T15;0 0 -T18/pi];

%% Matrices for the circulatory part of the damping and stiffness matrices
R_1=[-2*pi;2*pi*(a+1/2);-T12];
S_1=[0 1 T10/pi];
S_2=[1 0.5-a T11/(2*pi)];

%%

r=k_UVLM;
for ij=1:length(r);
    
%% Lift deficiency Factor
[ck_imag,ck_real,~,~]=Theodorsen_parameters2(r);
c_s=ck_real+1i*ck_imag;
   
   %% % Q matrix -  
   %Lift(row1), Pitching Moment(row2), column1(Heave Motion),
   % column2(Pitch Motion)
    ss=r*exp(1i*pi*theta);
    Q(:,:,1)=2*(b^2)*(M_nc*((ss)^2)+(B_nc+c_s*R_1*S_2)*(ss)+K_nc+c_s*R_1*S_1); 
    
end

end
