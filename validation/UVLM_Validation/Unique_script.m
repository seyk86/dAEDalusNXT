%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
tic
clear
% close all

% generate required paths for program execution
addpath(genpath('../../aerodynamics'));
addpath(genpath('../../geometry'));
addpath(genpath('../../aircraft'));
addpath(genpath('../../stdlib'));
addpath(genpath('../../input'));

UVLM_settings=class_UVLM_computation_settings();
UVLM_settings.debug=0;
UVLM_settings.wakelength_factor=0.9;
UVLM_settings.n_osc=8;
%minimum steps per period
UVLM_settings.spp=8;
aero_solver_settings=class_aero_solver_settings;

%% define aircraft geometry
aircraft=class_aircraft('AR40wing_mod.xml',1);
% aircraft.plot
% aircraft.weights=weights;
%     aircraft.grid_settings.x_max_grid_size=.05; %~40 panels for chord 2
%     aircraft.grid_settings.x_max_grid_size=.1; %~20 panels for chord 2
aircraft.grid_settings.x_max_grid_size=.125; %~16 panels for chord 2
%   aircraft.grid_settings.x_max_grid_size=.2; %~10 panels for chord 2
% aircraft.grid_settings.x_max_grid_size=.7;
%     aircraft.grid_settings.y_max_grid_size=30; % Span 100
aircraft.grid_settings.y_max_grid_size=40; % Span 100
% aircraft.grid_settings.x_max_grid_size=0.18;
% aircraft.grid_settings.y_max_grid_size=0.25;
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();

Uinf=50;
alpha=0;
beta=0;
Ma=0.0;

%     rho_air=0.397;
rho_air=1.225;
state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);
% wingaero.plot_grid;

% wingaero=wingaero.compute_influence_coeff_matrix();
wingaero=class_UVLM_solver(aircraft.name,aircraft.grid,aircraft.is_te,aircraft.panels,state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);

a=-0.4;
b=1;
c=1;

%% Amplitude

ah=1*pi/180;

hh=0.2;

%% Loop for Heave, Pitch (dAEDalus, Theodorsen, Garrick)
for ii=1:size(hh,2)
    
    %Reduced frequency range
    k_UVLM=[0.01 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0]; % Reduced frequency range
%     k_UVLM=[1 1.5 2];
%     k_UVLM=[0.01 0.25 0.5 1 1.5 2];
%     k_UVLM=[2];
    
    % Preallocation
    Q_UVLM=zeros(2,2,size(k_UVLM,2));
    Q=zeros(3,3,size(k_UVLM,2));
    Q_theodorsen=Q_UVLM;
    
    for kk=1:size(k_UVLM,2)
        % Solve heave
        wingaero = wingaero.solve_unsteady_heave(hh(ii), k_UVLM(kk));
        Q_UVLM(1,1,kk) = wingaero.Cl_complex(1)+1i*wingaero.Cl_complex(2);
        Q_UVLM(1,2,kk) = wingaero.CM_complex(1)+1i*wingaero.CM_complex(2);
        
        
        % Solve Garrick function
        generate_garrick_data_mod;
        
        
        % Plots for heave motion
        figure_handle = figure('Name',['Garrick Heave ' num2str(k_UVLM(kk)) ],'NumberTitle','off');
        hold on;
        plot(wingaero.alpha_s(end/2:end)*180/pi,wingaero.Cdi2(end/2:end)/2,'black')
        hold on
        plot(alpha_eff(end/2:end)*180/pi,Cd_heave(end/2:end),'b.')
        grid on
        xlabel('\alpha_{eff} [deg]')
        ylabel('C_d []');
        legend('dAEDalusNXT UVLM','Garrick');
        filename=['Garrick Heave ' num2str(k_UVLM(kk)) '.fig'];
        export_fig(filename)
        
        
        %Solve pitch
        wingaero=wingaero.solve_unsteady_pitch(ah*180/pi, k_UVLM(kk));
        Q_UVLM(2,1,kk)= wingaero.Cl_complex(1)+1i*wingaero.Cl_complex(2);
        Q_UVLM(2,2,kk)= wingaero.CM_complex(1)+1i*wingaero.CM_complex(2);
        
        
        % Solve Garrick function
        generate_garrick_data_mod;
        
        
        % Plots for pitch motion
        figure_handle = figure('Name',['Garrick Pitch ' num2str(k_UVLM(kk)) ],'NumberTitle','off');
        hold on;
        plot(wingaero.alpha_s(end/2:end)*180/pi,wingaero.Cdi(end/2:end)/2,'black')
        hold on
        plot(alpha_pitch(end/2:end)*180/pi,Cd_pitch(end/2:end),'b.')
        xlabel('\alpha [deg]')
        ylabel('C_d [ ]');
        grid on
        legend('dAEDalusNXT UVLM','Garrick');
        filename=['Garrick Pitch ' num2str(k_UVLM(kk)) '.fig'];
        export_fig(filename)
        
        
        % Theodorsen function
        [~,~,~,~,Q_theodorsen(:,:,kk)]=Theodorsen_Crepaldi(rho_air,a,b,c,Uinf,hh(ii),ah,k_UVLM(kk));
        
        
        % Theodorsen parameters from Bonin implementation
        [Q(:,:,kk)]=generate_ts_data_mod_function(a,b,c,k_UVLM(kk));
        Q(:,:,kk)=Q(:,:,kk)';
    end
    
end

Q_Daedalus=Q_UVLM;

%% Plot Real and Imaginary Part from Bonin implementation

figure;
subplot(2,2,1) % Heave, CZ
plot(k_UVLM,-real(squeeze(Q_Daedalus(1,1,:))/(hh/1)),'-r+')
hold on
plot(k_UVLM,imag(squeeze(-Q_Daedalus(1,1,:))/(hh/1)),'-b+')
xlabel('k')
ylabel('CZ')
grid on
plot(k_UVLM,real(squeeze(Q(1,1,:))),'-rx')
hold on
plot(k_UVLM,-imag(squeeze(Q(1,1,:))),'-bx')
legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
title('Heave')


subplot(2,2,3) % Heave, CM
plot(k_UVLM,real(squeeze(Q_Daedalus(1,2,:)/(hh))),'-r+')
hold on
plot(k_UVLM,imag(squeeze(Q_Daedalus(1,2,:)/(hh))),'-b+')
xlabel('k')
ylabel('CM')
grid on
plot(k_UVLM,real(squeeze(Q(1,2,:))),'-rx')
hold on
plot(k_UVLM,-imag(squeeze(Q(1,2,:))),'-bx')
legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
title('Heave')

subplot(2,2,2) % Pitch, CZ
plot(k_UVLM,real(squeeze(Q_Daedalus(2,1,:))/(ah)),'-r+')
hold on
plot(k_UVLM,imag(squeeze(Q_Daedalus(2,1,:))/(ah)),'-b+')
xlabel('k')
ylabel('CZ')
grid on
plot(k_UVLM,real(squeeze(Q(2,1,:))),'-rx')
hold on
plot(k_UVLM,-imag(squeeze(Q(2,1,:))),'-bx')
legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
title('Pitch')


subplot(2,2,4) % Pitch, CM
plot(k_UVLM,-real(squeeze(Q_Daedalus(2,2,:))/(ah)),'-r+')
hold on
plot(k_UVLM,imag(squeeze(-Q_Daedalus(2,2,:))/(ah)),'-b+')
xlabel('k')
ylabel('CM')
grid on
plot(k_UVLM,real(squeeze(Q(2,2,:))),'-rx')
hold on
plot(k_UVLM,-imag(squeeze(Q(2,2,:))),'-bx')
legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
title('Pitch')

%% Plot Magnitude and Phase from Bonin implementation

figure;
subplot(2,2,1) % Heave, CZ
plot(k_UVLM,abs(squeeze(Q_Daedalus(1,1,:))/(hh/1)),'-r+')
hold on
plot(k_UVLM,angle(squeeze(-Q_Daedalus(1,1,:))/(hh/1)),'-b+')
xlabel('k')
ylabel('CZ')
grid on
plot(k_UVLM,abs(squeeze(Q(1,1,:))),'-rx')
hold on
plot(k_UVLM,-angle(squeeze(Q(1,1,:))),'-bx')
legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
title('Heave')


subplot(2,2,3) % Heave, CM
plot(k_UVLM,abs(squeeze(Q_Daedalus(1,2,:)/(hh))),'-r+')
hold on
plot(k_UVLM,angle(squeeze(Q_Daedalus(1,2,:)/(hh))),'-b+')
xlabel('k')
ylabel('CM')
grid on
plot(k_UVLM,abs(squeeze(Q(1,2,:))),'-rx')
hold on
plot(k_UVLM,-angle(squeeze(Q(1,2,:))),'-bx')
legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
title('Heave')

subplot(2,2,2) % Pitch, CZ
plot(k_UVLM,abs(squeeze(Q_Daedalus(2,1,:))/(ah)),'-r+')
hold on
plot(k_UVLM,angle(squeeze(Q_Daedalus(2,1,:))/(ah)),'-b+')
xlabel('k')
ylabel('CZ')
grid on
plot(k_UVLM,abs(squeeze(Q(2,1,:))),'-rx')
hold on
plot(k_UVLM,-angle(squeeze(Q(2,1,:))),'-bx')
legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
title('Pitch')


subplot(2,2,4) % Pitch, CM
plot(k_UVLM,abs(squeeze(Q_Daedalus(2,2,:))/(ah)),'-r+')
hold on
plot(k_UVLM,angle(squeeze(-Q_Daedalus(2,2,:))/(ah)),'-b+')
xlabel('k')
ylabel('CM')
grid on
plot(k_UVLM,abs(squeeze(Q(2,2,:))),'-rx')
hold on
plot(k_UVLM,-angle(squeeze(Q(2,2,:))),'-bx')
legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
title('Pitch')

%% Error in Magnitude and Phase

%Relative Error
% Heave, CZ
Error_abs_cz_heave=abs((abs(Q_Daedalus(1,1,:)/hh)-abs(Q(1,1,:)))./abs(Q(1,1,:)));
Error_angle_cz_heave=abs((angle(-Q_Daedalus(1,1,:)/hh)+angle(Q(1,1,:)))./angle(Q(1,1,:)));

% Heave, CM
Error_abs_cm_heave=abs((abs(Q_Daedalus(1,2,:)/hh)-abs(Q(1,2,:)))./abs(Q(1,2,:)));
Error_angle_cm_heave=abs((angle(Q_Daedalus(1,2,:)/hh)+angle(Q(1,2,:)))./angle(Q(1,2,:)));

% Pitch, CZ
Error_abs_cz_pitch=abs((abs(Q_Daedalus(2,1,:)/ah)-abs(Q(2,1,:)))./abs(Q(2,1,:)));
Error_angle_cz_pitch=abs((angle(Q_Daedalus(2,1,:)/ah)+angle(Q(2,1,:)))./angle(Q(2,1,:)));

% Pitch, CM
Error_abs_cm_pitch=abs((abs(Q_Daedalus(2,2,:)/ah)-abs(Q(2,2,:)))./abs(Q(2,2,:)));
Error_angle_cm_pitch=abs((angle(-Q_Daedalus(2,2,:)/ah)+angle(Q(2,2,:)))./angle(Q(2,2,:)));

% % Absolute error
% % Heave, CZ
% Error_abs_cz_heave=abs((abs(Q_Daedalus(1,1,:)/hh)-abs(Q(1,1,:))));
% Error_angle_cz_heave=abs((angle(-Q_Daedalus(1,1,:)/hh)+angle(Q(1,1,:))));
% 
% % Heave, CM
% Error_abs_cm_heave=abs((abs(Q_Daedalus(1,2,:)/hh)-abs(Q(1,2,:))));
% Error_angle_cm_heave=abs((angle(Q_Daedalus(1,2,:)/hh)+angle(Q(1,2,:))));
% 
% % Pitch, CZ
% Error_abs_cz_pitch=abs((abs(Q_Daedalus(2,1,:)/ah)-abs(Q(2,1,:))));
% Error_angle_cz_pitch=abs((angle(Q_Daedalus(2,1,:)/ah)+angle(Q(2,1,:))));
% 
% % Pitch, CM
% Error_abs_cm_pitch=abs((abs(Q_Daedalus(2,2,:)/ah)-abs(Q(2,2,:))));
% Error_angle_cm_pitch=abs((angle(-Q_Daedalus(2,2,:)/ah)+angle(Q(2,2,:))));

figure(200)
hold on
subplot(2,2,1) 
hold on
plot(k_UVLM,squeeze(Error_abs_cz_heave(:,:,:)),'-b+')
xlabel('k')
ylabel('CZ Magnitude error')
title('Heave')
subplot(2,2,3) 
hold on
plot(k_UVLM,squeeze(Error_abs_cm_heave(:,:,:)),'-b+')
xlabel('k')
ylabel('CM Magnitude error')
title('Heave')
subplot(2,2,2) 
hold on
plot(k_UVLM,squeeze(Error_abs_cz_pitch(:,:,:)),'-b+')
xlabel('k')
ylabel('CZ Magnitude error')
title('Pitch')
subplot(2,2,4)
hold on
plot(k_UVLM,squeeze(Error_abs_cm_pitch(:,:,:)),'-b+')
xlabel('k')
ylabel('CM Magnitude error')
title('Pitch')

figure(201)
hold on
subplot(2,2,1) 
hold on
plot(k_UVLM,squeeze(Error_angle_cz_heave(:,:,:)),'-b+')
xlabel('k')
ylabel('CZ Phase error')
title('Heave')
subplot(2,2,3) 
hold on
plot(k_UVLM,squeeze(Error_angle_cm_heave(:,:,:)),'-b+')
xlabel('k')
ylabel('CM Phase error')
title('Heave')
subplot(2,2,2) 
hold on
plot(k_UVLM,squeeze(Error_angle_cz_pitch(:,:,:)),'-b+')
xlabel('k')
ylabel('CZ Phase error')
title('Pitch')
subplot(2,2,4)
hold on
plot(k_UVLM,squeeze(Error_angle_cm_pitch(:,:,:)),'-b+')
xlabel('k')
ylabel('CM Phase error')
title('Pitch')

figure(200)
subplot(2,2,1);legend('Type1','Type2')
subplot(2,2,2);legend('Type1','Type2')
subplot(2,2,3);legend('Type1','Type2')
subplot(2,2,4);legend('Type1','Type2')
figure(201)
subplot(2,2,1);legend('Type1','Type2')
subplot(2,2,2);legend('Type1','Type2')
subplot(2,2,3);legend('Type1','Type2')
subplot(2,2,4);legend('Type1','Type2')


%% Error in Real and Imaginary
% 
% % Relative Error
% % Heave, CZ
% Error_real_cz_heave=abs((-real(Q_Daedalus(1,1,:)/hh)-real(Q(1,1,:)))./real(Q(1,1,:)));
% Error_imag_cz_heave=abs((imag(-Q_Daedalus(1,1,:)/hh)+imag(Q(1,1,:)))./imag(Q(1,1,:)));
% 
% % Heave, CM
% Error_real_cm_heave=abs((real(Q_Daedalus(1,2,:)/hh)-real(Q(1,2,:)))./real(Q(1,2,:)));
% Error_imag_cm_heave=abs((imag(Q_Daedalus(1,2,:)/hh)+imag(Q(1,2,:)))./imag(Q(1,2,:)));
% 
% % Pitch, CZ
% Error_real_cz_pitch=abs((real(Q_Daedalus(2,1,:)/ah)-real(Q(2,1,:)))./real(Q(2,1,:)));
% Error_imag_cz_pitch=abs((imag(Q_Daedalus(2,1,:)/ah)+imag(Q(2,1,:)))./imag(Q(2,1,:)));
% 
% % Pitch, CM
% Error_real_cm_pitch=abs((-real(Q_Daedalus(2,2,:)/ah)-real(Q(2,2,:)))./real(Q(2,2,:)));
% Error_imag_cm_pitch=abs((imag(-Q_Daedalus(2,2,:)/ah)+imag(Q(2,2,:)))./imag(Q(2,2,:)));
% 
% % %Absolute Error
% % % Heave, CZ
% % Error_real_cz_heave=abs((-real(Q_Daedalus(1,1,:)/hh)-real(Q(1,1,:))));
% % Error_imag_cz_heave=abs((imag(-Q_Daedalus(1,1,:)/hh)+imag(Q(1,1,:))));
% % 
% % % Heave, CM
% % Error_real_cm_heave=abs((real(Q_Daedalus(1,2,:)/hh)-real(Q(1,2,:))));
% % Error_imag_cm_heave=abs((imag(Q_Daedalus(1,2,:)/hh)+imag(Q(1,2,:))));
% % 
% % % Pitch, CZ
% % Error_real_cz_pitch=abs((real(Q_Daedalus(2,1,:)/ah)-real(Q(2,1,:))));
% % Error_imag_cz_pitch=abs((imag(Q_Daedalus(2,1,:)/ah)+imag(Q(2,1,:))));
% % 
% % % Pitch, CM
% % Error_real_cm_pitch=abs((-real(Q_Daedalus(2,2,:)/ah)-real(Q(2,2,:))));
% % Error_imag_cm_pitch=abs((imag(-Q_Daedalus(2,2,:)/ah)+imag(Q(2,2,:))));
% 
% figure(202)
% hold on
% subplot(2,2,1) 
% hold on
% plot(k_UVLM,squeeze(Error_real_cz_heave(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CZ Real error')
% title('Heave')
% subplot(2,2,3) 
% hold on
% plot(k_UVLM,squeeze(Error_real_cm_heave(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CM Real error')
% title('Heave')
% subplot(2,2,2) 
% hold on
% plot(k_UVLM,squeeze(Error_real_cz_pitch(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CZ Real error')
% title('Pitch')
% subplot(2,2,4)
% hold on
% plot(k_UVLM,squeeze(Error_real_cm_pitch(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CM Real error')
% title('Pitch')
% 
% figure(203)
% hold on
% subplot(2,2,1) 
% hold on
% plot(k_UVLM,squeeze(Error_imag_cz_heave(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CZ Imaginary error')
% title('Heave')
% subplot(2,2,3) 
% hold on
% plot(k_UVLM,squeeze(Error_imag_cm_heave(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CM Imaginary error')
% title('Heave')
% subplot(2,2,2) 
% hold on
% plot(k_UVLM,squeeze(Error_imag_cz_pitch(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CZ Imaginary error')
% title('Pitch')
% subplot(2,2,4)
% hold on
% plot(k_UVLM,squeeze(Error_imag_cm_pitch(:,:,:)),'-bx')
% xlabel('k')
% ylabel('CM Imaginary error')
% title('Pitch')

%% Plot Real and Imaginary values (from Theodorsen_Crepaldi)

% chord=2*b;
% figure;
% subplot(2,2,1) % Heave, CZ
% plot(k_UVLM,-real(squeeze(Q_Daedalus(1,1,:))/(hh/1)),'-r+')
% hold on
% plot(k_UVLM,-imag(squeeze(Q_Daedalus(1,1,:))/(hh/1)),'-b+')
% xlabel('k')
% ylabel('CZ')
% grid on
% plot(k_UVLM,real(squeeze(Q_theodorsen(1,1,:))),'-rx')
% hold on
% plot(k_UVLM,imag(squeeze(Q_theodorsen(1,1,:))),'-bx')
% legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
% title('Heave')
%
% subplot(2,2,3) % Heave, CM
% plot(k_UVLM,real(squeeze(Q_Daedalus(1,2,:)/(hh))),'-r+')
% hold on
% plot(k_UVLM,imag(squeeze(Q_Daedalus(1,2,:)/(hh))),'-b+')
% xlabel('k')
% ylabel('CM')
% grid on
% plot(k_UVLM,real(squeeze(Q_theodorsen(1,2,:))),'-rx')
% hold on
% plot(k_UVLM,imag(squeeze(Q_theodorsen(1,2,:))),'-bx')
% legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
% title('Heave')
%
% subplot(2,2,2) % Pitch, CZ
% plot(k_UVLM,real(squeeze(Q_Daedalus(2,1,:))/(ah)),'-r+')
% hold on
% plot(k_UVLM,imag(squeeze(Q_Daedalus(2,1,:))/(ah)),'-b+')
% xlabel('k')
% ylabel('CZ')
% grid on
% plot(k_UVLM,real(squeeze(Q_theodorsen(2,1,:))),'-rx')
% hold on
% plot(k_UVLM,imag(squeeze(Q_theodorsen(2,1,:))),'-bx')
% legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
% title('Pitch')
%
% subplot(2,2,4) % Pitch, CM
% plot(k_UVLM,-real(squeeze(Q_Daedalus(2,2,:))/(ah)),'-r+')
% hold on
% plot(k_UVLM,-imag(squeeze(Q_Daedalus(2,2,:))/(ah)),'-b+')
% xlabel('k')
% ylabel('CM')
% grid on
% plot(k_UVLM,real(squeeze(Q_theodorsen(2,2,:))),'-rx')
% hold on
% plot(k_UVLM,imag(squeeze(Q_theodorsen(2,2,:))),'-bx')
% legend('Re dAEDalus','Imag dAEDalus','Re Theodorsen','Imag Theodorsen')
% title('Pitch')
%
%
% %% Plot Magnitude and Phase
%
% figure;
% subplot(2,2,1) % Heave, CZ
% plot(k_UVLM,abs(squeeze(Q_Daedalus(1,1,:))/(hh/1)),'-r+')
% hold on
% plot(k_UVLM,angle(squeeze(-Q_Daedalus(1,1,:))/(hh/1)),'-b+')
% xlabel('k')
% ylabel('CZ')
% grid on
% plot(k_UVLM,abs(squeeze(Q_theodorsen(1,1,:))),'-rx')
% hold on
% plot(k_UVLM,angle(squeeze(Q_theodorsen(1,1,:))),'-bx')
% legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
% title('Heave')
%
%
% subplot(2,2,3) % Heave, CM
% plot(k_UVLM,abs(squeeze(Q_Daedalus(1,2,:)/(hh))),'-r+')
% hold on
% plot(k_UVLM,angle(squeeze(Q_Daedalus(1,2,:)/(hh))),'-b+')
% xlabel('k')
% ylabel('CM')
% grid on
% plot(k_UVLM,abs(squeeze(Q_theodorsen(1,2,:))),'-rx')
% hold on
% plot(k_UVLM,angle(squeeze(Q_theodorsen(1,2,:))),'-bx')
% legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
% title('Heave')
%
% subplot(2,2,2) % Pitch, CZ
% plot(k_UVLM,abs(squeeze(Q_Daedalus(2,1,:))/(ah)),'-r+')
% hold on
% plot(k_UVLM,angle(squeeze(Q_Daedalus(2,1,:))/(ah)),'-b+')
% xlabel('k')
% ylabel('CZ')
% grid on
% plot(k_UVLM,abs(squeeze(Q_theodorsen(2,1,:))),'-rx')
% hold on
% plot(k_UVLM,angle(squeeze(Q_theodorsen(2,1,:))),'-bx')
% legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
% title('Pitch')
%
%
% subplot(2,2,4) % Pitch, CM
% plot(k_UVLM,abs(squeeze(Q_Daedalus(2,2,:))/(ah)),'-r+')
% hold on
% plot(k_UVLM,angle(squeeze(-Q_Daedalus(2,2,:))/(ah)),'-b+')
% xlabel('k')
% ylabel('CM')
% grid on
% plot(k_UVLM,abs(squeeze(Q_theodorsen(2,2,:))),'-rx')
% hold on
% plot(k_UVLM,angle(squeeze(Q_theodorsen(2,2,:))),'-bx')
% legend('Magnit dAEDalus','Phase dAEDalus','Magnit Theodorsen','Phase Theodorsen')
% title('Pitch')

toc


