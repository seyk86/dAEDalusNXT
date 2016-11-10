%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
clear all
%close all
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
UVLM_settings.spp=16;
aero_solver_settings=class_aero_solver_settings;

%% define aircraft geometry
%aircraft=class_aircraft('AGARD445.6Wing.xml',1);

%aircraft=class_aircraft('optimwing.xml',1)

aircraft=class_aircraft('AR40wing.xml',1);
% aircraft.plot
% aircraft.weights=weights;
aircraft.grid_settings.x_max_grid_size=0.13;
aircraft.grid_settings.y_max_grid_size=40;
% aircraft.grid_settings.x_max_grid_size=0.18;
% aircraft.grid_settings.y_max_grid_size=0.25;
aircraft.grid_settings.wake=2;
aircraft=aircraft.compute_grid();
  
Uinf=10;
alpha=0;
beta=0;
Ma=0.0;

rho_air=0.397;
state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);

% wingaero.plot_grwid;
 tic
% wingaero=wingaero.compute_influence_coeff_matrix();

wingaero=class_UVLM_solver(aircraft.name,aircraft.grid,aircraft.is_te,aircraft.panels,state,aircraft.grid_wake,aircraft.panels_wake,aircraft.reference,UVLM_settings);

k_UVLM=[0.1 1];
 hh=[0.2 0.02];
 ah=1*pi/180;
 aa=-0.5;
 
 Q_UVLM=zeros(3,3,length( k_UVLM));
 for kk=1:length( k_UVLM)
     wingaero=wingaero.solve_unsteady_heave(hh(kk), k_UVLM(kk));
     generate_garrick_data;
     figure_handle = figure('Name',['garrick_heave' num2str(k_UVLM(kk)) ],'NumberTitle','off');
     hold on;
     set(gcf,'PaperUnits','centimeters')
     xSize = 7.5; ySize = 7.5;
     xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
     set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
     set(gcf,'Position',[200 50 xSize*50 ySize*50])
     pos = get(gcf, 'DefaultAxesPosition');
     set(gcf, 'DefaultAxesPosition',[pos(1)-0.06, pos(2), pos(3)-0.06, pos(4)]);
     plot(wingaero.alpha_s(end/2:end)*180/pi,wingaero.Cdi2(end/2:end)/2,'black')
     hold on
     plot(alpha_eff(end/2:end)*180/pi,Cd_heave(end/2:end),'b.')
     grid on
     xlim([-1.2 1.2])
     xlabel('\alpha_{eff} [deg]')
     ylabel('C_d []');
    % ylim([-19E-4 2E-4])
     legend('dAEDalusNXT UVLM','Garrick');
     filename=['garrick_heave' num2str(k_UVLM(kk)) '.fig'];
     export_fig(filename)
     %close(figure_handle);
    
     wingaero=wingaero.solve_unsteady_pitch(ah*180/pi, k_UVLM(kk));
     generate_garrick_data;
     figure_handle = figure('Name',['garrick_pitch' num2str(k_UVLM(kk)) ],'NumberTitle','off');
     hold on;
     set(gcf,'PaperUnits','centimeters')
     xSize = 7.5; ySize = 7.5;
     xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
     set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
     set(gcf,'Position',[200 50 xSize*50 ySize*50])
     pos = get(gcf, 'DefaultAxesPosition');
     set(gcf, 'DefaultAxesPosition',[pos(1)-0.06, pos(2), pos(3)-0.06, pos(4)]);
     plot(wingaero.alpha_s(end/2:end)*180/pi,wingaero.Cdi(end/2:end)/2,'black')
     hold on
     plot(alpha_pitch(end/2:end)*180/pi,Cd_pitch(end/2:end),'b.')
     xlabel('\alpha [deg]')
     ylabel('C_d [ ]');
     grid on
     xlim([-1 1])
  %   ylim([-19E-4 2E-4])
     legend('dAEDalusNXT UVLM','Garrick');
     filename=['garrick_pitch' num2str(k_UVLM(kk)) '.fig'];
     export_fig(filename)
    % close(figure_handle);
 
 end
