%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Uinf] = a2bf(V_A,alpha,beta,Ma_corr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 Uinf=V_A*[cosd(alpha*Ma_corr)*cosd(beta*Ma_corr)  sind(beta*Ma_corr)*cosd(alpha*Ma_corr) sind(alpha*Ma_corr)];

end

