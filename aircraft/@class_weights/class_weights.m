%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_weights
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        MZFW=0;
        MTOW=0;
        MLW=0;
        W=0;
        OWE=0;
        Fuel=0;
        
        
        WingSystemsEstimate=0;
        WingSkinEstimate=0; 
        WingInitialGuess=0;
        WingTotal=0;
        
        
        FuselagePAX;
        Fuselage
        FuselageNonStructuralEstimate=0;          % Includes payload, crew equipment and accommodation accessories
        FuselageSystemsEstimate=0;
        
        PassengerWeight=80; %weight of person in kg
        BaggageWeight=30; % weight of baggage in kg
        SeatWeight=10; %weight of passenger aset in kg

    end
    
    methods
    end
    
end

