%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [output_args] = compare_free_flight_models(Models)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


legendtxt= ['legend(h,''' Models{1}.Name ''''];
for i=2:length(Models)    
legendtxt=[legendtxt ',''' Models{i}.Name ''''];
end
legendtxt=[legendtxt ')'];


hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Xe(:,1),Models{i}.Xe(:,3))
end
hold off
eval(legendtxt)

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Euler(:,1));
end
eval(legendtxt)
ylabel('\Phi [rad/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Euler(:,2));
end
eval(legendtxt)
ylabel('\Theta [rad/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Euler(:,3));
end
eval(legendtxt)
ylabel('\Psi [rad/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.pqr(:,1));
end
eval(legendtxt)
ylabel('p [rad/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.pqr(:,2));
end
eval(legendtxt)
ylabel('q [rad/s]')
xlabel('t [s]')


figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.pqr(:,3));
end
eval(legendtxt)
ylabel('r [rad/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Alpha*180/pi);
end
eval(legendtxt)
ylabel('Alpha [�]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Beta*180/pi);
end
eval(legendtxt)

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.V);
end
eval(legendtxt)
ylabel('V [m/s]')
xlabel('t [s]')

% 

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.qinf(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('qinf[]')
xlabel('t [s]')
% 

% figure
% hold on
% for i=1:length(Models)    
% h(i)=plot(Models{i}.Time,Models{i}.CX(1:length(Models{i}.Time)));
% end
% eval(legendtxt)
% ylabel('C_X []')
% xlabel('t [s]')
% 
% 
figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.CZ(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('C_Z []')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.CY(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('C_Y []')
xlabel('t [s]')


figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.CM(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('C_M []')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.CL(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('C_L []')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.CN(1:length(Models{i}.Time)));
end
eval(legendtxt)
ylabel('C_N []')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Vb(1:length(Models{i}.Time),1));
end
eval(legendtxt)
ylabel('V_b [m/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Vb(1:length(Models{i}.Time),3));
end
eval(legendtxt)
ylabel('V_b [m/s]')
xlabel('t [s]')


figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Vb(1:length(Models{i}.Time),1));
end
eval(legendtxt)
ylabel('V_e [m/s]')
xlabel('t [s]')

figure
hold on
for i=1:length(Models)    
h(i)=plot(Models{i}.Time,Models{i}.Vb(1:length(Models{i}.Time),3));
end
eval(legendtxt)
ylabel('V_e [m/s]')
xlabel('t [s]')


% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.Vb,'r')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.Vb(1:length(MeanAxisModal.Time),:),'c')
% plot(MeanAxisPartitioned.Time,MeanAxisPartitioned.Vb(1:length(MeanAxisPartitioned.Time),:),'k')
% 
% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.Ve,'r')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.Ve(1:length(MeanAxisModal.Time),:),'c')
% plot(MeanAxisPartitioned.Time,MeanAxisPartitioned.Ve(1:length(MeanAxisPartitioned.Time),:),'k')
% 
% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.V,'r')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.V(1:length(MeanAxisModal.Time)),'c')
% plot(FixedAxis.Time,FixedAxis.V(1:length(FixedAxis.Time)),'g')
% plot(MeanAxisID.Time,MeanAxisID.V(1:length(MeanAxisID.Time)),'y')
% plot(MeanAxisPartitioned.Time,MeanAxisPartitioned.V(1:length(MeanAxisPartitioned.Time)),'k')
% 
% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.elastic_forces(:,5)+MeanAxisHybrid.rigid_forces(:,5),'g')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.CM(2:end)*1/2*1.225.*MeanAxisModal.V(1:length(MeanAxisModal.Time))'.^2*aircraft.reference.c_ref*aircraft.reference.S_ref,'b')
% 
% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.elastic_forces(:,3)+MeanAxisHybrid.rigid_forces(:,3),'g')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.CZ(2:end)*1/2*1.225.*MeanAxisModal.V(1:length(MeanAxisModal.Time))'.^2*aircraft.reference.S_ref,'b')
% 
% figure
% plot(MeanAxisHybrid.Time,MeanAxisHybrid.rigid_forces(:,3),'g')
% hold on
% plot(MeanAxisModal.Time,MeanAxisModal.CZ(2:end)*1/2*1.225.*MeanAxisModal.V(1:length(MeanAxisModal.Time))'.^2*aircraft.reference.S_ref,'b')
% plot(MeanAxisPartitioned.Time,MeanAxisPartitioned.CZ(1:end)*1/2*1.225.*MeanAxisPartitioned.V(1:length(MeanAxisPartitioned.Time))'.^2*aircraft.reference.S_ref,'k')
% 

end

