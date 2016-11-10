%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%close all
clear all
addpath(genpath('../../structures'));

%number of beam elements
nel=10;
totalLength=10; % [m]
phi=-20*pi/180; %backward swept
nu=20*pi/180;   %
twist=40*pi/180; %in rad!

%material properties
E=6.8900e+10;
G=2.6000e+10;

%circular Beam Crosssection Properties
d=0.01; %[m]
r=d/2;
A=pi*r^2;
Ix=A/4*r^2;  %main bending
Iz=Ix/2;     %inplane
J = (pi*r^3)/4; % torsional moment for a circular bar
Ip=J; 
%torsional constant
%polar moment of inertia

% calculate node coords
nodeCoordsLocal=zeros(3,nel+1);
nodeCoordsLocal(2,:)=linspace(0,totalLength,nel+1);

a=nu; % around x
b=twist; % around y
c=phi; % around z

Lx=[1       0       0
    0   cos(a)  sin(a)
    0   -sin(a) cos(a)];%defines negative rotation around x

Ly=[cos(b) 0 -sin(b)
    0      1    0
    sin(b)  0   cos(b)];%defines negative rotation around y

Lz=[cos(c) sin(c)   0
    -sin(c) cos(c)  0
    0           0   1]; %defines negative rotation around z


T=Lz*Ly*Lx;
nodeCoords=T'*nodeCoordsLocal %backward trafo "from local to global"

zVectorGlobal=[T'*[0;0;1]]';



%initialize beam 
beam=class_beam(nel,'class_crosssection_wingbox');
beam.node_coords=nodeCoords;
for i=1:nel
    beam.beamelement(i)=beam.beamelement(i).setElementGeometry(totalLength/nel,phi,nu,twist);
    beam.beamelement(i).Ix=Ix;
    beam.beamelement(i).Ip=Ip;
    beam.beamelement(i).Iz=Iz;
    beam.beamelement(i).J=J;
    beam.beamelement(i).A=A;
    beam.beamelement(i).E=E;
    beam.beamelement(i).G=G;
    beam.beamelement(i).A_enclosed=0;
end
beam.update_K=1;
beam.update_Q=1;
beam.update_M=1;

%build beam_collection
structure=class_beam_collection(beam);

%set loads
%set distributed load (defined in local coordinate system)
localload=T*[1;1;1];
for i=1:nel
    structure.beam.beamelement(i).qx=localload(1); %[N/m]
    structure.beam.beamelement(i).qy=localload(2); %[N/m]
    structure.beam.beamelement(i).qz=localload(3); %[N/m]
    structure.beam.beamelement(i).mt=0.1; %[Nm/m] 
end
%set force
%structure.beam(1).nodal_loads(end-3)=10;
structure.beam(1).update_Q=1;

%set boundary condition
structure.beam(1)=structure.beam(1).f_add_boundary_condition(class_boundary_condition(1,[1 1 1 1 1 1],[0 0 0 0 0 0]));

%initialize structural solversettings
structure.settings=class_wingstructure_solver_settings;

%solve
structure=structure.f_solve;

%read results
nasData=read_pch('C:\Daten\Nastran_Workdir\beam_valid\beamvalid.pch');

dAEDdata=reshape(structure.nodal_deflections,6,10)';
dAEDdata=[zeros(1,6); dAEDdata];

err=(dAEDdata-nasData)./nasData;

%plot
figure;
subplot(2,3,1)
grid on
hold on
plot(dAEDdata(:,1),'-o');
plot(nasData(:,1),'-+');
subplot(2,3,2)
grid on
hold on
plot(dAEDdata(:,2),'-o');
plot(nasData(:,2),'-+');
subplot(2,3,3)
grid on
hold on
plot(dAEDdata(:,3),'-o');
plot(nasData(:,3),'-+');
subplot(2,3,4)
grid on
hold on
plot(dAEDdata(:,4),'-o');
plot(nasData(:,4),'-+');
subplot(2,3,5)
grid on
hold on
plot(dAEDdata(:,5),'-o');
plot(nasData(:,5),'-+');
subplot(2,3,6)
grid on
hold on
plot(dAEDdata(:,6),'-o');
plot(nasData(:,6),'-+');


figure;
subplot(2,3,1)
plot(err(:,1),'-o');
grid on
hold on
subplot(2,3,2)
plot(err(:,2),'-o');
grid on
hold on
subplot(2,3,3)
plot(err(:,3),'-o');
grid on
hold on
subplot(2,3,4)
plot(err(:,4),'-o');
grid on
hold on
subplot(2,3,5)
plot(err(:,5),'-o');
grid on
hold on
subplot(2,3,6)
plot(err(:,6),'-o');
grid on
hold on


