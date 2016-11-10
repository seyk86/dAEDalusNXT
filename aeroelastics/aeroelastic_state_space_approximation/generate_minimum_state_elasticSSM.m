%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [AeroelasticSSM]=generate_minimum_state_elasticSSM(A0,A1,A2,E,D,Ms,Ks,aircraft,gamma,V,rho,AeroelasticSSM)
% UNTITLED2 Creates an Aeroelastic State Space Model from Roger's
% Approximation
% Detailed explanation goes here
 %fullstructure=fullstructure.f_assemble_free(1,0);
nG=AeroelasticSSM.nG;
nC=AeroelasticSSM.nC;
bb=0.5*aircraft.reference.c_ref;
nQ=size(A0,1);
nL=length(gamma);
nR=AeroelasticSSM.nR;
nE=length(Ms);

A0elastic=A0(nR+1:nR+nE,nR+1:nR+nE);
AeroelasticSSM.A0ee=A0elastic;
A1elastic=A1(nR+1:nR+nE,nR+1:nR+nE);
AeroelasticSSM.A1ee=A1elastic;
A2elastic=A2(nR+1:nR+nE,nR+1:nR+nE);
AeroelasticSSM.A2ee=A2elastic;
Eelastic=E(:,nR+1:nR+nE);
AeroelasticSSM.Eelastic=Eelastic;
Delastic=D(nR+1:nR+nE,:);
AeroelasticSSM.Delastic=Delastic;


A0rigid=A0(nR+1:nR+nE,1:nR);
AeroelasticSSM.A0re=A0rigid;
A1rigid=A1(nR+1:nR+nE,1:nR);
AeroelasticSSM.A1re=A1rigid;
A2rigid=A2(nR+1:nR+nE,1:nR);
AeroelasticSSM.A2re=A2rigid;
Erigid=E(:,1:nR);
AeroelasticSSM.Erigid=Erigid;
Drigid=D(1:nR,:);
AeroelasticSSM.Drigid=Drigid;

A0elastic_rigid=A0(1:nR,nR+1:nR+nE);
AeroelasticSSM.A0er=A0elastic_rigid;
A1elastic_rigid=A1(1:nR,nR+1:nR+nE);
AeroelasticSSM.A1er=A1elastic_rigid;
A2elastic_rigid=A2(1:nR,nR+1:nR+nE);
AeroelasticSSM.A2er=A2elastic_rigid;

A0control=A0(nR+1:nR+nE,nR+1+nE:nR+nE+nC);
AeroelasticSSM.A0ce=A0control;
A1control=A1(nR+1:nR+nE,nR+1+nE:nR+nE+nC);
AeroelasticSSM.A1ce=A1control;
A2control=A2(nR+1:nR+nE,nR+1+nE:nR+nE+nC);
AeroelasticSSM.A2ce=A2control;
Econtrol=E(:,nR+1+nE:nR+nE+nC);
AeroelasticSSM.Econtrol=Econtrol;
Dcontrol=D(nR+1+nE:nR+nE+nC,:);

Kaee=-0.5*aircraft.reference.S_ref*rho*(V^2)*A0elastic+Ks;
AeroelasticSSM.Ks=Ks;
Caee=-0.5*aircraft.reference.S_ref*rho*V*bb*A1elastic;
Maee=-0.5*aircraft.reference.S_ref*rho*bb^2*A2elastic+Ms;
AeroelasticSSM.Ms=Ms;
Kar=0.5*aircraft.reference.S_ref*rho*(V^2)*A0rigid;
Car=0.5*aircraft.reference.S_ref*rho*V*bb*A1rigid;
Mar=0.5*aircraft.reference.S_ref*rho*bb^2*A2rigid;

Kac=0.5*aircraft.reference.S_ref*rho*(V^2)*A0control;
Cac=0.5*aircraft.reference.S_ref*rho*V*bb*A1control;
Mac=0.5*aircraft.reference.S_ref*rho*bb^2*A2control;

Ker=0.5*aircraft.reference.S_ref*rho*(V^2)*A0elastic_rigid;
Cer=0.5*aircraft.reference.S_ref*rho*V*bb*A1elastic_rigid;
Mer=0.5*aircraft.reference.S_ref*rho*bb^2*A2elastic_rigid;
Aer=Drigid*0.5*aircraft.reference.S_ref*rho*(V^2);


A=[zeros(nE,nE), eye(nE,nE),zeros(nE,nL),zeros(nE,nL),zeros(nE,nL)
    -Maee^(-1)*Kaee , -Maee^(-1)*Caee,Maee^(-1)*Delastic*0.5*aircraft.reference.S_ref*rho*(V^2),Maee^(-1)*Delastic*0.5*aircraft.reference.S_ref*rho*(V^2),Maee^(-1)*Delastic*0.5*aircraft.reference.S_ref*rho*(V^2);
    zeros(nL,nE),Eelastic,-diag(gamma)*V/bb,zeros(nL,nL),zeros(nL,nL);
    zeros(nL,nE),Eelastic*0,zeros(nL,nL),-diag(gamma)*V/bb,zeros(nL,nL);
    zeros(nL,nE),Eelastic*0,zeros(nL,nL),zeros(nL,nL),-diag(gamma)*V/bb;];

B_rigid=[zeros(nE,nR), zeros(nE,nR),zeros(nE,nR);
    Maee^(-1)*Kar , Maee^(-1)*Car,Maee^(-1)*Mar;...
    zeros(nL,3*nR);
    zeros(nL,nR),Erigid,zeros(nL,nR);
    zeros(nL,3*nR);];

B_control=[zeros(nE,nC), zeros(nE,nC),zeros(nE,nC);
    Maee^(-1)*Kac ,Maee^(-1)*Cac,Maee^(-1)*Mac;...
            zeros(nL,3*nC);
            zeros(nL,3*nC);
    zeros(nL,nC),Econtrol,zeros(nL,nC);];

C_force=[Ker Cer Aer];

AeroelasticSSM.msASee=A;
AeroelasticSSM.msASce=B_control;
AeroelasticSSM.msASre=B_rigid;
AeroelasticSSM.msCer=C_force;
AeroelasticSSM.Kaee=Kaee;
end

