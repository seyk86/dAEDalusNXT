%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%AeroelasticSSM
clear all
load('input_data_for_test_dynamics.mat')
n_modes=16;
[AeroelasticSSM]=generate_SSM_from_GAF(Q,k,aircraft_structure_ow,aircraft_ow,n_modes,50,1.225);
%% INIT
V_max=50;
Aircraft.AeroelasticSSM=AeroelasticSSM;
Aircraft.Reference=aircraft_ow.reference;
Aircraft.rho=1.225;
Aircraft.V=50;
t_sim=5;
aeroflag=0;
%Gust def
H=0.1;
Uds=10;
t_start=0.5;

s=0:(2*H/(200)):2*H;
t_d=2*H/V_max;

gust_shape=[0 Uds/2*(1-cos(pi*s/H))];
dot_gust_shape=[0 diff(gust_shape)];
dot_dot_gust_shape=[0 diff(dot_gust_shape)];

clear rigid_state
rigid_state.time=[t_start:t_d/(length(gust_shape)-1):t_start+t_d];%wingaero_free_rigid.t_vec;%
rigid_state.signals.values=zeros(22,size(AeroelasticSSM.ASre,2));
rigid_state.signals.values=zeros(length(rigid_state.time),size(AeroelasticSSM.ASre,2));
% control def

%sine doublet
t_start=0.5;
delta_freq=15;
delta_amp=5*pi/180;
n_sines=20;

%control surface inertia calculation
csMass=2;
csCG=sum(squeeze(aircraft_ow.wings(1).wing_segments(6).xyz_te_device),2)/4;
hingeLineInner=aircraft_ow.wings(1).wing_segments(6).xyz_te_device(:,:,1)';
hingeLineOuter=aircraft_ow.wings(1).wing_segments(6).xyz_te_device(:,:,2)';
hingeDirection=(hingeLineOuter-hingeLineInner)/(norm(hingeLineOuter-hingeLineInner));
csInertiaHingelineGlobal=hingeDirection*csMass*norm(abs(cross(hingeLineOuter-hingeLineInner,csCG-hingeLineInner))/abs(hingeLineOuter-hingeLineInner))^2;
%n�chster Strukturknoten
loadPoint=aircraft_structure_ow.node_coords(dsearchn(aircraft_structure_ow.node_coords, hingeLineInner'),:)';
dPH=hingeLineInner-loadPoint;
csInertiaLoadPointGlobal=[0;0;0;csInertiaHingelineGlobal;];


                           
Control_surface_inertia=zeros(612,11);
Control_surface_inertia(aircraft_structure_ow.dof_node_beam(:,2)==dsearchn(aircraft_structure_ow.node_coords, hingeLineInner'),7)=csInertiaLoadPointGlobal;
CsInertiaModal=Aircraft.AeroelasticSSM.modeshapes'*Control_surface_inertia;



%% Sim
sim('test_dynamics',t_sim);

%% read pch
fileID=fopen('C:\Users\rmkseywa\Desktop\Nastran_Workdir\nas_unsteady_moment\checks\unsteady_forcet.pch','r');
nasData=[];
tline = fgets(fileID);
i=0;
nasLine=[];
while ischar(tline)
    if tline(1)~='$'
        if tline(1:6)~='-CONT-'
            if i~=0
                 nasData(i,:)=nasLine;
            end
            i=i+1;
            nasLine=str2double(strsplit(tline));
        else 
            nasLine=[nasLine str2double(strsplit(tline))];
        end
    end
    tline = fgets(fileID);
end
fclose(fileID);

%%
figure;
hold on
grid on

 plot(nasData(:,2), nasData(:,11))
 plot(physical_states.Time, -physical_states.Data(:,611))
 legend('nas','daed')
 
