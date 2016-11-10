%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function MeanAxisHybridSimulink=simulate_free_flying_hybrid_simulink(modelname,aircraft,Init_Flight_State,Aircraft,AeroelasticSSM,ctrl_input,t_vec,maneuver_name)


mkdir(['results/' aircraft.name],maneuver_name);
%% generate required paths for program execution
% initialize flight state variables
MeanAxisHybridSimulink.Xe=zeros(length(t_vec),3);
MeanAxisHybridSimulink.Euler=zeros(length(t_vec),3);
MeanAxisHybridSimulink.Vb=zeros(length(t_vec),3);
MeanAxisHybridSimulink.Ve=zeros(length(t_vec),3);
MeanAxisHybridSimulink.pqr=zeros(length(t_vec),3);
MeanAxisHybridSimulink.pqr_dot=zeros(length(t_vec),3);
MeanAxisHybridSimulink.Alpha=zeros(length(t_vec),1);
MeanAxisHybridSimulink.Beta=zeros(length(t_vec),1);
MeanAxisHybridSimulink.V=zeros(length(t_vec),1);
MeanAxisHybridSimulink.qinf=zeros(length(t_vec),1);
MeanAxisHybridSimulink.Ma=zeros(length(t_vec),1);


[Environment]=f_environment_init();
Environment.ACTIVATE_GUST=0;

control_surface_states_0=Init_Flight_State.AircraftState.ControlDeflections;
control_surface_states=zeros(length(aircraft.control_deflections),1);

ctrl_state=zeros(1,length(control_surface_states)*3);
ctrl_state_0=ctrl_state;
for n_ctrl=1:length(control_surface_states)
   ctrl_state_0(n_ctrl)=Init_Flight_State.AircraftState.ControlDeflections(n_ctrl);
end 

SimMode=1;
control_in.time=ctrl_input.time;
control_in.signals.values=zeros(length(control_in.time),ctrl_input.signals.dimensions+4+7);
control_in.signals.dimensions=ctrl_input.signals.dimensions+4+7;

for n_ctrl=1:length(control_surface_states)
    control_in.signals.values(:,2+n_ctrl)=ctrl_input.signals.values(:,n_ctrl)+Init_Flight_State.AircraftState.ControlDeflections(n_ctrl);
end


control_in.signals.values(:,1)=Init_Flight_State.AircraftState.DeltaThr;
control_in.signals.values(:,2)=Init_Flight_State.AircraftState.DeltaThr;


MeanAxisHybridSimulink.Alpha(1)=Init_Flight_State.Alpha;
%ASE=AeroelasticSSM.msASee^(-1)*(-AeroelasticSSM.msASre(:,:)*[0 0 0 0 Init_Flight_State.Alpha  0   zeros(1,6)  zeros(1,6)]'-AeroelasticSSM.msASce*(ctrl_state_0')*pi/180);

Aircraft.AeroelasticSSM=AeroelasticSSM;
%Init_Flight_State.x_stat_ms=AeroelasticSSM.msASee^-1*1/2*Init_Flight_State.rho*Init_Flight_State.VCAS^2*aircraft.reference.S_ref*[zeros(AeroelasticSSM.nE,1);AeroelasticSSM.Q_0(7:6+AeroelasticSSM.nE,1);zeros(size(AeroelasticSSM.msASre(:,:),1)-2*AeroelasticSSM.nE,1)];

Aircraft.ModelName=[aircraft.name '_sim'];
%Aircraft.TrimSurfaces=aircraft.control_surfaces;
Aircraft.TrimSurfaces=aircraft.control.trim_surfaces;

save([aircraft.name '/' aircraft.name '_SIMDATA'],'Aircraft','Aircraft');
%Init_Flight_State.AeroelasticState=ASE;

assignin('base','Init_Flight_State',Init_Flight_State);
assignin('base','Environment',Environment);
assignin('base','control_in',control_in);
assignin('base','SimMode',SimMode);
assignin('base','Aircraft',Aircraft);
% found no other solution.. T_FLAP is needed in integrated simulink
% m-function but the struct AIRCRAFT is too big to be processed in the
% m-file (crash)..
assignin('base','T_FLAP',Aircraft.T_FLAP)

time=sim(modelname,max(t_vec));
MeanAxisHybridSimulink.Time=time;
MeanAxisHybridSimulink.Xe=XeOut.Data;
MeanAxisHybridSimulink.Euler=EulerOut.Data;
MeanAxisHybridSimulink.Vb=VbOut.Data;
MeanAxisHybridSimulink.Ve=VeOut.Data;
MeanAxisHybridSimulink.pqr=pqrOut.Data;

MeanAxisHybridSimulink.Gamma=gammaOut.Data;
%MeanAxisHybridSimulink.pqr_dot=pqr_dotOut.Data;
MeanAxisHybridSimulink.Alpha=AlphaOut.Data;
MeanAxisHybridSimulink.Beta=-BetaOut.Data;
MeanAxisHybridSimulink.V=VOut.Data;
%MeanAxisHybridSimulink.qinf=qinfOut.Data;
%MeanAxisHybridSimulink.Ma=MaOut.Data;
MeanAxisHybridSimulink.CX=-aero_coeff.Data(:,1);
MeanAxisHybridSimulink.CZ=-aero_coeff.Data(:,3);
MeanAxisHybridSimulink.CM=aero_coeff.Data(:,5);
MeanAxisHybridSimulink.CY=aero_coeff.Data(:,2);
MeanAxisHybridSimulink.CL=-aero_coeff.Data(:,4);
MeanAxisHybridSimulink.CN=-aero_coeff.Data(:,6);
%MeanAxisHybridSimulink.CM=momentout.Data(:,2)./(qinfOut.Data*aircraft.reference.S_ref*aircraft.reference.c_ref);
%MeanAxisHybridSimulink.CZtot=forceout.Data(:,3)./(qinfOut.Data*aircraft.reference.c_ref*aircraft.reference.S_ref);
%MeanAxisHybridSimulink.x_AE=generalized_states.Data;
%MeanAxisHybridSimulink.Aeroelastic_Force_Out=Aeroelastic_Force_Out;
%MeanAxisHybridSimulink.Gravity_Force_Out=gravityforce;
save(['results/' aircraft.name '/' maneuver_name '/MeanAxisHybridSimulink'],'MeanAxisHybridSimulink','MeanAxisHybridSimulink');
% aeroblock_in.time=MeanAxisPartitionedQS.Time;
% aeroblock_in.signals.values=zeros(length(MeanAxisPartitionedQS.Time),8);
% aeroblock_in.signals.dimensions=8;
% aeroblock_in.signals.values(:,1)=MeanAxisPartitionedQS.pqr(:,1);
% aeroblock_in.signals.values(:,2)=MeanAxisPartitionedQS.pqr(:,2);
% aeroblock_in.signals.values(:,3)=MeanAxisPartitionedQS.pqr(:,3);
% aeroblock_in.signals.values(:,4)=MeanAxisPartitionedQS.Beta;
% aeroblock_in.signals.values(:,5)=MeanAxisPartitionedQS.Alpha;
% aeroblock_in.signals.values(:,7)=MeanAxisPartitionedQS.V;
% aeroblock_in.signals.values(:,8)=1.225;
