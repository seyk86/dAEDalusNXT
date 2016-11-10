%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [AeroelasticSSM]=generate_SSM_from_GAF(Q,k,aircraft_structure,aircraft,n_modes,V_max,rho_air)

%% assuming control modes are at the end, gust modes are not implemented yet! careful!
Q=Q([1:n_modes end-length(aircraft.control_surfaces)+1:end],[1:n_modes end-length(aircraft.control_surfaces)+1:end],:);

if aircraft_structure.modefrequencies(7)/aircraft_structure.modefrequencies(6)>100
    nR=6;%always 6 rigid body modes
    constrained=0;
else
    nR=0;
    constrained=1;
end
nE=n_modes-nR;
nC=length(aircraft.control_surfaces);
nG=0;



modeshapes_rom=real(aircraft_structure.modeshapes)*[eye(n_modes); zeros(size(aircraft_structure.modeshapes,1)-n_modes,n_modes)];

%Mm=modeshapes_rom(1:end,7:n_modes)'*aircraft_structure.Mff(1:end,1:end)*modeshapes_rom(1:end,7:n_modes);
Mm=eye(nE,nE);
if constrained
    Km=modeshapes_rom(1:end,1:n_modes)'*aircraft_structure.Kff(1:end,1:end)*modeshapes_rom(1:end,1:n_modes);
                 
    M_tot_mean=Mm;
    Km=eye(nE,nE).*Km;
    K_tot_mean=Km;
else
    Km=modeshapes_rom(1:end,7:n_modes)'*aircraft_structure.Kff(1:end,1:end)*modeshapes_rom(1:end,7:n_modes);
    [~,~,I_Hat,~,b_Hat_Skew,~,A_Bar]=compute_eqm_matrices_mean([0 0 0]',[0 0 0]',aircraft_structure.node_coords,aircraft_structure.nodal_deflections*0,[0 0 0]',[0 0 0]');
    M_inert=A_Bar*aircraft_structure.Mff*A_Bar';
    Mass=I_Hat'*M_inert*I_Hat;
    Inertia=b_Hat_Skew'*M_inert*b_Hat_Skew;
    zw1=I_Hat'*A_Bar*aircraft_structure.Mff;
    zw2=b_Hat_Skew'*A_Bar*aircraft_structure.Mff;
    zw3=modeshapes_rom(1:end,7:n_modes)'*aircraft_structure.Mff*A_Bar'*I_Hat;
    zw4=modeshapes_rom(1:end,7:n_modes)'*aircraft_structure.Mff*A_Bar'*b_Hat_Skew;

    M_tot_mean=[    Mass                            zeros(3,3)                             zeros(3,length(Mm));
                    zeros(3,3)                      Inertia                                zeros(3,length(Mm));
                    zw3                             zw4                                    Mm];
                
    Km=eye(nE,nE).*Km;
    K_tot_mean=M_tot_mean*0;
    zw5=I_Hat'*A_Bar*aircraft_structure.Kff*modeshapes_rom(1:end,7:n_modes);
    zw6=b_Hat_Skew'*A_Bar*aircraft_structure.Kff*modeshapes_rom(1:end,7:n_modes);
    zw7=Km;
    K_tot_mean(:,7:end)=[zw5;zw6;zw7];
end
 
    



       
mode_idx=[1:n_modes size(Q,1)-nC+1:size(Q,1)];

Q=Q(mode_idx,mode_idx,1:end);
k=k(1:end);

gamma=k;

[A0, A1, A2, Arest] = rogers_state_space_approximation(Q,k,gamma);
[A0ms,A1ms,A2ms,Ems,Dms] = minimum_state_approximation(Q,k,gamma,40);

nL=length(gamma);
AeroelasticSSM.gamma=gamma;
AeroelasticSSM.nR=nR;
AeroelasticSSM.nE=nE;
AeroelasticSSM.nL=nL;
AeroelasticSSM.nC=nC;
AeroelasticSSM.nG=nG;

AeroelasticSSM.A0ms=A0ms;
AeroelasticSSM.A1ms=A1ms;
AeroelasticSSM.A2ms=A2ms;
AeroelasticSSM.Dms=Dms;
AeroelasticSSM.Ems=Ems;

AeroelasticSSM.A0=A0;
AeroelasticSSM.A1=A1;
AeroelasticSSM.A2=A2;
AeroelasticSSM.Arest=Arest;

AeroelasticSSM.Ms=M_tot_mean(1:nR+nE,1:nR+nE);
AeroelasticSSM.Ks=K_tot_mean(1:nR+nE,1:nR+nE);

AeroelasticSSM.bb=0.5*aircraft.reference.c_ref;
AeroelasticSSM.S_ref=aircraft.reference.S_ref;

AeroelasticSSM.modeshapes=aircraft_structure.modeshapes(:,nR+1:nR+nE);
AeroelasticSSM.acc_dof_idx=[1 2 3 4 5 6];


for V=1:0.5:V_max
%    [AS,ASre,ASee,ASce,ASge,ASrr,AScc,ASgg,Mae,Cae,Kae,Alag,Ks,Ms,ASpre]=generate_roger_SSM_ext(A0,A1,A2,Arest,M_tot_mean(1:nR+nE,1:nR+nE),K_tot_mean(1:nR+nE,1:nR+nE),aircraft,gamma,V,rho_air,nC,nG);
%    [msAS,msASrr,msASre,msASce,msASge,msASee,Kae,Cae,Mae]=generate_minimum_state_SSM(A0ms,A1ms,A2ms,Ems,Dms,M_tot_mean(1:nR+nE,1:nR+nE),K_tot_mean(1:nR+nE,1:nR+nE),aircraft,gamma,V,rho_air,nC,nG);
    [AeroelasticSSM]=generate_minimum_state_elasticSSM(A0ms,A1ms,A2ms,Ems,Dms,M_tot_mean(nR+1:nR+nE,nR+1:nR+nE),K_tot_mean(nR+1:nR+nE,nR+1:nR+nE),aircraft,gamma,V,rho_air,AeroelasticSSM);
    [AeroelasticSSM]=generate_roger_elasticSSM(A0,A1,A2,Arest,M_tot_mean(nR+1:nR+nE,nR+1:nR+nE),K_tot_mean(nR+1:nR+nE,nR+1:nR+nE),aircraft,gamma,V,rho_air,AeroelasticSSM);
    %[msAS,msASre,msASce,msASge,msASee,msKae,msCae,msMae,msA0rs,msASrs,msAScrs,msASgrs,msASrr]=minimum_state_SSM_coupled_dynamics(A0ms,A1ms,A2ms,Ems,Dms,[1 2 3 4 5 6],mode_idx(7:end-nG-nC),mode_idx(end-nG-nC+1:end-nG),mode_idx(end-nG+1:end),gamma,length(gamma),aircraft,aircraft_state,Mse,Kse,V_start,rho,0,0);
    eigenvAA=eig(AeroelasticSSM.ASee);
     eigenvmsAA=eig(AeroelasticSSM.msASee);
    
    %[EVec,eigenvAA]=eig(AS);
    %sys=ss(AS,ones(length(AS),1),ones(1,length(AS)),0)
%     plot(real(eigenvAA), imag(eigenvAA),'o','MarkerFaceColor',[V/V_max 0 1-V/V_max],'MarkerEdgeColor',[V/V_max 0 1-V/V_max])
    hold on
%     eigenvAA=eig(msAn);
%     %sys=ss(AS,ones(length(AS),1),ones(1,length(AS)),0);
     plot(real(eigenvmsAA), imag(eigenvmsAA),'o','MarkerFaceColor',[V/V_max*0.5 1-V/V_max 0],'MarkerEdgeColor',[V/V_max*0.5  1-V/V_max 0])
     if any(real(eigenvmsAA)>0)
         disp(['re>0, V ',num2str(V)]);
     end
end

AeroelasticSSM.selGAC=eye(nE,length(AeroelasticSSM.msASee));
AeroelasticSSM.selGACdot=[zeros(nE,nE) eye(nE,length(AeroelasticSSM.msASee)-nE)];

%mkdir(aircraft.name)
%save([aircraft.name '/AESSM'],'AeroelasticSSM');
