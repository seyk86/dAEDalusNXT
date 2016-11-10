%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function x=Lsfit_osc(omega,time,y)

% COSINE least squares fitting for the input function

%% lsfit parameters
options = optimset('TolX',1E-20,'TolFun',1E-20,'Display','off');
amp_lim=1e-8;
LB=[0 -pi -Inf];
UB=[Inf pi Inf];

%% Initial
Init_Amplitude=(abs(max(y(length(time)/4:end)))+abs(min(y(length(time)/4:end))))/2;
Init_Phase=-1.5;
Init_Offset=0.5*max(y(length(time)/4:end))+0.5*min(y(length(time)/4:end));

%% Fitting
x = lsqcurvefit(@(x,time) x(1)*cos(omega*time+x(2))+x(3),[Init_Amplitude Init_Phase Init_Offset],time(length(time)/2:end),y(length(time)/2:end),LB,UB,options);

if x(1)<amp_lim
    x(1)=0;
    x(2)=0;
end

end
