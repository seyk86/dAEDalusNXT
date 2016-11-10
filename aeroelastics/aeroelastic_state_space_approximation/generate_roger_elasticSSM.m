%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [AeroelasticSSM]=generate_roger_elasticSSM(A0,A1,A2,Arest,Ms,Ks,aircraft,gamma,V,rho,AeroelasticSSM)
% UNTITLED2 Creates an Aeroelastic State Space Model from Roger's
% Approximation
% Detailed explanation goes here
 %fullstructure=fullstructure.f_assemble_free(1,0);
bb=0.5*aircraft.reference.c_ref;
%nQ=size(A0,1)+6;
nQ=size(A0,1);
nL=length(gamma);
nR=AeroelasticSSM.nR;
nE=length(Ms);
nG=AeroelasticSSM.nG;
nC=AeroelasticSSM.nC;

A0elastic=A0(nR+1:nR+nE,nR+1:nR+nE);
A1elastic=A1(nR+1:nR+nE,nR+1:nR+nE);
A2elastic=A2(nR+1:nR+nE,nR+1:nR+nE);

A0elastic_rigid=A0(1:nR,nR+1:nR+nE);
A1elastic_rigid=A1(1:nR,nR+1:nR+nE);
A2elastic_rigid=A2(1:nR,nR+1:nR+nE);

for rac_c=1:nL
    Arest_elastic(:,1+(rac_c-1)*nE:nE*rac_c)=Arest(nR+1:nR+nE,1+nR+nQ*(rac_c-1):(rac_c-1)*nQ+nR+nE);
    Arest_elastic_rigid(:,1+(rac_c-1)*nE:nE*rac_c)=Arest(1:nR,1+nR+nQ*(rac_c-1):(rac_c-1)*nQ+nR+nE);
    Arest_rigid(:,1+(rac_c-1)*nR:nR*rac_c)=Arest(nR+1:nR+nE,1+nQ*(rac_c-1):(rac_c-1)*nQ+nR);
    Arest_control(:,1+(rac_c-1)*nC:nC*rac_c)=Arest(nR+1:nR+nE,1+nR+nE+nQ*(rac_c-1):(rac_c-1)*nQ+nR+nE+nC); 
end

A0rigid=A0(nR+1:nR+nE,1:nR);
A1rigid=A1(nR+1:nR+nE,1:nR);
A2rigid=A2(nR+1:nR+nE,1:nR);

A0control=A0(nR+1:nR+nE,nR+1+nE:nR+nE+nC);
A1control=A1(nR+1:nR+nE,nR+1+nE:nR+nE+nC);
A2control=A2(nR+1:nR+nE,nR+1+nE:nR+nE+nC);

E=-0.5*aircraft.reference.S_ref*rho*(V^2)*A0elastic+Ks;
D=-0.5*aircraft.reference.S_ref*rho*V*bb*A1elastic;
A=-0.5*aircraft.reference.S_ref*rho*bb^2*A2elastic*0+Ms;
Alag=0.5*aircraft.reference.S_ref*rho*V^2*Arest_elastic;

Kar=0.5*aircraft.reference.S_ref*rho*(V^2)*A0rigid;
Car=0.5*aircraft.reference.S_ref*rho*V*bb*A1rigid;
Mar=0.5*aircraft.reference.S_ref*rho*bb^2*A2rigid;
Alag_rigid=0.5*aircraft.reference.S_ref*rho*V^2*Arest_rigid;

Ker=0.5*aircraft.reference.S_ref*rho*(V^2)*A0elastic_rigid;
Cer=0.5*aircraft.reference.S_ref*rho*V*bb*A1elastic_rigid;
Mer=0.5*aircraft.reference.S_ref*rho*bb^2*A2elastic_rigid;
Aer=0.5*aircraft.reference.S_ref*rho*(V^2)*Arest_elastic_rigid;

Kac=0.5*aircraft.reference.S_ref*rho*(V^2)*A0control;
Cac=0.5*aircraft.reference.S_ref*rho*V*bb*A1control;
Mac=0.5*aircraft.reference.S_ref*rho*bb^2*A2control;
Alag_control=0.5*aircraft.reference.S_ref*rho*V^2*Arest_control;


AS=zeros((nE)*(nL+2)+nR*nL+nC*nL,(nE)*(nL+2)+nR*nL+nC*nL);

AS=[zeros(nE,nE),eye(nE,nE) ,A^(-1)*Alag*0,A^(-1)*Alag_rigid*0,A^(-1)*Alag_control*0;...
    -A^(-1)*E,-A^(-1)*D ,A^(-1)*Alag,A^(-1)*Alag_rigid,A^(-1)*Alag_control;];
for i=1:nL
    AS=[AS;...
        zeros(nE,nE),eye(nE,nE),A^(-1)*Alag*0,A^(-1)*Alag_rigid*0,A^(-1)*Alag_control*0;];
    AS(2*nE+(i-1)*nE+1:2*nE+i*nE,2*nE+(i-1)*nE+1:2*nE+i*nE)=AS(2*nE+(i-1)*nE+1:2*nE+i*nE,2*nE+(i-1)*nE+1:2*nE+i*nE)-V/bb*gamma(i)*eye(nE,nE);
end

start_idx=nE*(nL+2);
for i=1:nL
    AS=[AS;...
        zeros(nR,(nE)*(nL+2)+nR*nL+nC*nL)];
    AS(start_idx+1+(i-1)*nR:start_idx+i*nR,start_idx+1+(i-1)*nR:start_idx+i*nR)=AS(start_idx+1+(i-1)*nR:start_idx+i*nR,start_idx+1+(i-1)*nR:start_idx+i*nR)-V/bb*gamma(i)*eye(nR,nR);
end

start_idx=(nE)*(nL+2)+nR*nL;
for i=1:nL
    AS=[AS;...
        zeros(nC,(nE)*(nL+2)+nR*nL+nC*nL)];
    AS(start_idx+1+(i-1)*nC:start_idx+i*nC,start_idx+1+(i-1)*nC:start_idx+i*nC)=AS(start_idx+1+(i-1)*nC:start_idx+i*nC,start_idx+1+(i-1)*nC:start_idx+i*nC)-V/bb*gamma(i)*eye(nC,nC);
end

B_rigid=[zeros(nE,nR),zeros(nE,nR) ,A^(-1)*Mar*0;...
    A^(-1)*Kar,A^(-1)*Car ,A^(-1)*Mar;
    zeros(nE*nL,nR*3)];
for i=1:nL
    B_rigid=[B_rigid;...
        zeros(nR,nR),eye(nR,nR),zeros(nR,nR)];
end
    B_rigid=[B_rigid;...
        zeros(nC*nL,nR),zeros(nC*nL,nR),zeros(nC*nL,nR)];
    
    
B_control=[zeros(nE,nC),zeros(nE,nC) ,A^(-1)*Mac*0;...
    A^(-1)*Kac,A^(-1)*Cac ,A^(-1)*Mac;
    zeros(nE*nL+nR*nL,nC*3)];
for i=1:nL
    B_control=[B_control;...
        zeros(nC,nC),eye(nC,nC),zeros(nC,nC)];
end


C_force=[Ker Cer Aer];

AeroelasticSSM.ASee=AS;
AeroelasticSSM.ASce=B_control;
AeroelasticSSM.ASre=B_rigid;
AeroelasticSSM.Cer=C_force;
end

