%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Rigid_Aerodata] = compute_aerotable_rigid(aircraft,ac_state,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
state=ac_state.aerodynamic_state;
deflected=0;
control_surfs=aircraft.control_surfaces;
drag=1;
if nargin==2
    aircraft=aircraft.compute_grid();
    config=0;
else
    deflected=1;
    def=varargin{1};
    aircraft=aircraft.compute_deflected_grid(def);
    config=1;
end

vecConf=1:1:length(ac_state);
vecDeltaR=-30:10:30;
vecDeltaP=-30:10:30;
vecDeltaQ=-30:10:30;
vecCsDefl=-30:5:30;
vecAlpha=-10:2:15;
%vecAlpha=-10:5:15;
vecMa=0.0:0.05:0.2;
%vecMa=0.0:0.2:0.6;
vecBeta=-10:5:10;
%vecBeta=-10:10:10;
veciH=-20:5:20;
conf=0;
Rigid_Aerodata.Static_Coeff.CX=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_Aerodata.Static_Coeff.CY=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_Aerodata.Static_Coeff.CZ=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_Aerodata.Static_Coeff.Cl=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_Aerodata.Static_Coeff.Cm=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_Aerodata.Static_Coeff.Cn=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

CXp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CYp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CZp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Clp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cmp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cnp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

CXq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CYq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CZq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Clq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cmq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cnq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

CXr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CYr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
CZr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Clr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cmr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Cnr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

%counterst
i=1;
j=1;
k=1;
l=1;
tic

if deflected
    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference); 
else
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference); 
end
wingaero=wingaero.f_solve_std();
wingaero=wingaero.f_solve_full();
t_stop=toc

tic
wingaero=wingaero.f_set_state(state);
wingaero=wingaero.set_grid(aircraft.grid,aircraft.panels);
%wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
wingaero=wingaero.f_solve_std();
wingaero=wingaero.f_solve_full();
t_stop2=toc

%time estimation
t_loop=t_stop;
t_loop2=t_stop2;
t_total=0;
%for conf=vecConf
for Ma=vecMa
    for alpha=vecAlpha
        for beta=vecBeta
            t_total=t_total+t_loop;
        end
    end
end
%end

%for conf=vecConf
for Ma=vecMa
    for alpha=vecAlpha
        for beta=vecBeta
            for iH=veciH
                t_total=t_total+t_loop2;
            end
            for DeltaQ=vecDeltaQ
                t_total=t_total+t_loop2;
            end
            for DeltaP=vecDeltaP
                t_total=t_total+t_loop2;
            end
            for DeltaR=vecDeltaR
                t_total=t_total+t_loop2;
            end
            for DeltaR=vecDeltaR
                t_total=t_total+t_loop2;
            end
        end
    end
end
%end

fprintf('Estimated computation time: %f hours \n',t_total/3600);
%vecConf=1;
% static derivatives

Uinf=200;
rho_air=0.7;

for conf_cs=1:length(aircraft.control_surfaces)
    aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{conf_cs},0);
end
if deflected    
    aircraft=aircraft.compute_grid();
    aircraft=aircraft.compute_deflected_grid(def);
    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
else
    aircraft=aircraft.compute_grid();
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
end


coeff_lim=1E-10;
for conf=vecConf
    k=1;
    if length(ac_state)>1
        for conf_cs=1:length(ac_state(conf).aircraft_state.control_deflections)
            aircraft=aircraft.f_set_control_surface(ac_state(conf).aircraft_state.control_surfaces{conf_cs},ac_state(conf).aircraft_state.control_deflections(conf_cs)); %set all control surfaces to zero!!
        end
    end
    if deflected
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(def);
    else
        aircraft=aircraft.compute_grid();
    end

    for Ma=vecMa
        j=1;
        state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);
        aircraft=aircraft.compute_CD_f(state,aircraft.reference.S_ref);
        
        for alpha=vecAlpha
            i=1;
            for beta=vecBeta
                disp(['Ma=' num2str(Ma) ' aoa=' num2str(alpha) ' beta=' num2str(beta)])
                state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);
                if deflected
                    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                else
                    
                    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                end
                %  wingaero=wingaero.f_set_state(state);

                wingaero=wingaero.f_solve_std();
                wingaero=wingaero.f_solve_full();
                if drag==0
                    wingaero.Cdi=0;
                    aircraft.CD_f=0;
                else
                    %wingaero.Cdi=wingaero.compute_trefftz_drag_mex();
                end
                
                C= [wingaero.CX, wingaero.CY, wingaero.Cl, wingaero.CL, wingaero.CM, wingaero.CN];
                Cp=[wingaero.CXp,wingaero.CYp,wingaero.CZp,wingaero.CLp,wingaero.CMp,wingaero.CNp];
                Cq=[wingaero.CXq,wingaero.CYq,wingaero.CZq,wingaero.CLq,wingaero.CMq,wingaero.CNq];
                Cr=[wingaero.CXr,wingaero.CYr,wingaero.CZr,wingaero.CLr,wingaero.CMr,wingaero.CNr];
                % static derivatives
                % static force derivatives
                Rigid_Aerodata.Static_Coeff.CX(l,k,j,i)=-C(1);
                Rigid_Aerodata.Static_Coeff.CY(l,k,j,i)=C(2);
                Rigid_Aerodata.Static_Coeff.CZ(l,k,j,i)=-C(3);
                % static moment derivatives
                Rigid_Aerodata.Static_Coeff.Cl(l,k,j,i)=-C(4);
                Rigid_Aerodata.Static_Coeff.Cm(l,k,j,i)=C(5);
                Rigid_Aerodata.Static_Coeff.Cn(l,k,j,i)=-C(6);
                % damping derivatives
                % p derivatives
                Rigid_Aerodata.Damping_Coeff.CXp(l,k,j,i)=-Cp(1);
                Rigid_Aerodata.Damping_Coeff.CYp(l,k,j,i)=Cp(2);
                Rigid_Aerodata.Damping_Coeff.CZp(l,k,j,i)=-Cp(3);
                Rigid_Aerodata.Damping_Coeff.Clp(l,k,j,i)=Cp(4);
                Rigid_Aerodata.Damping_Coeff.Cmp(l,k,j,i)=Cp(5);
                Rigid_Aerodata.Damping_Coeff.Cnp(l,k,j,i)=Cp(6);
                
                % q derivatives
                Rigid_Aerodata.Damping_Coeff.CXq(l,k,j,i)=-Cq(1);%*cosd(alpha)+Cq(3)*sind(alpha);
                Rigid_Aerodata.Damping_Coeff.CYq(l,k,j,i)=Cq(2);
                Rigid_Aerodata.Damping_Coeff.CZq(l,k,j,i)=-Cq(3);%*cosd(alpha)-Cq(1)*sind(alpha);
                Rigid_Aerodata.Damping_Coeff.Clq(l,k,j,i)=-Cq(4);
                Rigid_Aerodata.Damping_Coeff.Cmq(l,k,j,i)=Cq(5);
                Rigid_Aerodata.Damping_Coeff.Cnq(l,k,j,i)=-Cq(6);
                
                % p derivatives
                Rigid_Aerodata.Damping_Coeff.CXr(l,k,j,i)=-Cr(1);
                Rigid_Aerodata.Damping_Coeff.CYr(l,k,j,i)=Cr(2);
                Rigid_Aerodata.Damping_Coeff.CZr(l,k,j,i)=-Cr(3);
                Rigid_Aerodata.Damping_Coeff.Clr(l,k,j,i)=Cr(4);
                Rigid_Aerodata.Damping_Coeff.Cmr(l,k,j,i)=Cr(5);
                Rigid_Aerodata.Damping_Coeff.Cnr(l,k,j,i)=Cr(6);
                i=i+1;
                %                 C;
                %                 Ma;
                %                 alpha;
            end
            j=j+1;
        end
        k=k+1;
    end
    l=l+1;
end


if length(ac_state)==1
    Rigid_Aerodata.Static_Coeff.CX(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CX(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CY(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CY(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CZ(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CZ(1,:,:,:);
    % static moment derivatives
    Rigid_Aerodata.Static_Coeff.Cl(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cl(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cm(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cm(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cn(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cn(1,:,:,:);
    % damping derivatives
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnp(1,:,:,:);
    
    % q derivatives
    Rigid_Aerodata.Damping_Coeff.CXq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXq(1,:,:,:);%*cosd(alpha)+Cq(3)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.CYq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZq(1,:,:,:);%*cosd(alpha)-Cq(1)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.Clq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnq(1,:,:,:);
    
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnr(1,:,:,:);
end

% control surface derivatives
h=1;
i=1;
j=1;
k=1;
l=1;

CX=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));
CY=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));
CZ=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));
Cl=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));
Cm=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));
Cn=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(vecDeltaQ));

Delta=zeros(1,length(aircraft.control_surfaces));
for conf_cs=1:length(aircraft.control_surfaces)
    aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{conf_cs},0);
end

if deflected
    aircraft=aircraft.compute_grid();
    aircraft=aircraft.compute_deflected_grid(def);
else
    aircraft=aircraft.compute_grid();
end
for conf=vecConf
    k=1;
    if config~=1
        for conf_cs=1:length(config(conf).controls)
            aircraft=aircraft.f_set_control_surface(config(conf).controls{conf_cs},config(conf).deflections*0);
        end
    end
    if deflected
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(def);
    else
        aircraft=aircraft.compute_grid();
    end
    
    for Ma=vecMa
        j=1;
        for alpha=vecAlpha
            i=1;
            for beta=vecBeta
                h=1;
                state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);
                if deflected
                    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                else
                    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                end
                % wingaero=wingaero.f_set_state(state);
                %wingaero=wingaero.set_grid(aircraft.grid,aircraft.panels);
                
                wingaero=wingaero.f_solve_std();
                wingaero=wingaero.f_solve_full();
                if drag==0
                    wingaero.Cdi=0;
                end
                C0= [wingaero.CX, wingaero.CY, wingaero.CZ, wingaero.CL, wingaero.CM, wingaero.CN];
                for cs=1:length(control_surfs)
                   
                    for cs_defl=vecCsDefl
                disp(['Ma=' num2str(Ma) ' aoa=' num2str(alpha) ' beta=' num2str(beta) ' cs:' control_surfs{cs} ' csdefl=' num2str(cs_defl)])
                        aircraft=aircraft.f_set_control_surface(control_surfs{cs},cs_defl);
                        if deflected
                            aircraft=aircraft.compute_grid();
                            aircraft=aircraft.compute_deflected_grid(def);
                            wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                        else
                            aircraft=aircraft.compute_grid();
                            wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
                        end
                        wingaero=wingaero.f_solve_std();
                        if drag==0
                            wingaero.Cdi=0;
                            aircraft.CD_f=0;
                        end
                        %wingaero=wingaero.f_solve_full();
                        C= [wingaero.CX, wingaero.CY, wingaero.CZ, wingaero.CL, wingaero.CM, wingaero.CN]-C0;
                        % static force derivatives
                        CX(l,k,j,i,h)=-C(1);
                        CY(l,k,j,i,h)=C(2);
                        CZ(l,k,j,i,h)=-C(3);
                        % static moment derivatives
                        Cl(l,k,j,i,h)=-C(4);
                        Cm(l,k,j,i,h)=C(5);
                        Cn(l,k,j,i,h)=-C(6);
                        h=h+1;
                        
                    end
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX(l,k,j,i,:)=CX(l,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY(l,k,j,i,:)=CY(l,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ(l,k,j,i,:)=CZ(l,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl(l,k,j,i,:)=Cl(l,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm(l,k,j,i,:)=Cm(l,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn(l,k,j,i,:)=Cn(l,k,j,i,:)']);
                    
                    aircraft=aircraft.f_set_control_surface(control_surfs{cs},0);
                    if deflected
                        aircraft=aircraft.compute_grid();
                        aircraft=aircraft.compute_deflected_grid(def);
                    else
                        aircraft=aircraft.compute_grid();
                    end
                    h=1;
                    control_surfs{cs}
                    Ma
                end
                i=i+1;
            end
            j=j+1;
        end
        k=k+1;
        k
    end
    l=l+1;
    l
end

if config==1
    k=1;
    l=2;
    for Ma=vecMa
        j=1;
        for alpha=vecAlpha
            i=1;
            for beta=vecBeta
                h=1;
                
                for cs=1:length(control_surfs)
                    
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX(1,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY(1,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ(1,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl(1,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm(1,k,j,i,:)']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn(1,k,j,i,:)']);
                    h=1;
                end
                i=i+1;
            end
            j=j+1;
        end
        k=k+1;
        k
    end
end

%% %% ground effect and landing gear dummy
Rigid_Aerodata.Landing_Gear.CX_gear=zeros(length(vecAlpha),1);
Rigid_Aerodata.Landing_Gear.CZ_gear=zeros(length(vecAlpha),1);

Rigid_Aerodata.Ground_Effect.CX_ground=zeros(length(vecConf),length(vecAlpha));
Rigid_Aerodata.Ground_Effect.CZ_ground=zeros(length(vecConf),length(vecAlpha));
Rigid_Aerodata.Ground_Effect.Cm_ground=zeros(length(vecConf),length(vecAlpha));

Xi_fact=zeros(40,1);
l=1;
landing_gear_down=1;
h=0;
for conf=vecConf
    i=1;
    for alpha=vecAlpha
        Rigid_Aerodata.Ground_Effect.CZ(l,i)=0;
        Rigid_Aerodata.Ground_Effect.CX(l,i)=0;
        Rigid_Aerodata.Ground_Effect.Cm(l,i)=0;
        i=i+1;
    end
    l=l+1;
end

if length(ac_state)==1
    i=1;
    for alpha=vecAlpha
        Rigid_Aerodata.Ground_Effect.CZ(2,i)= Rigid_Aerodata.Ground_Effect.CZ(1,i);
        Rigid_Aerodata.Ground_Effect.CX(2,i)=Rigid_Aerodata.Ground_Effect.CX(1,i);
        Rigid_Aerodata.Ground_Effect.Cm(2,i)=Rigid_Aerodata.Ground_Effect.Cm(1,i);
        i=i+1;
    end
end

i=1;
for alpha=vecAlpha
    Rigid_Aerodata.Landing_Gear.CZ(i)=0;
    Rigid_Aerodata.Landing_Gear.CX(i)=0;
    i=i+1;
end

for h=[0:1:40]
    Xi_fact(h+1)=0;
end
Rigid_Aerodata.Ground_Effect.Xi=Xi_fact;
Rigid_Aerodata.Ground_Effect.vecH=0:1:40;
Rigid_Aerodata.Landing_Gear.Xi=Xi_fact;
Rigid_Aerodata.Landing_Gear.vecH=0:1:40;

Rigid_Aerodata.Grid.vecAlpha=vecAlpha;
Rigid_Aerodata.Grid.vecMa=vecMa;
Rigid_Aerodata.Grid.vecBeta=vecBeta;
if config==1
    Rigid_Aerodata.Grid.vecConf=0:1;
else
    Rigid_Aerodata.Grid.vecConf=vecConf-1;
end
%%% TODO!!! automate!!!
for cs=1:length(control_surfs)
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Grid= [' num2str(vecCsDefl) ']']);
end
% vecDeltaSP=vecDeltaR;
mkdir(aircraft.name)
save([aircraft.name '/Rigid_Aerodata'],'Rigid_Aerodata');

end

