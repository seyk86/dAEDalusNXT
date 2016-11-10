%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%------------------------------------------------------------------------
%
% Calc_Starting_point_in_WGS84.m
%
% Calculates the WGS84 coordinates of a desired starting point
% defined by a given distance from an airport. This function is mainly used
% for finding landing starting points, since often the desired distance from
% an airport is known, but the corresponding WGS coordinates are unknown.
% However, FlightSim requests a starting point defined by WGS coordinates
%
% INPUT:
% --obligatory
% - x_rw: desired distance of starting point from the runway reference
%         point in (negative) runway direction in meters
% - y_rw: desired lateral offset of the starting point from the runway
%         reference point in meters. A positive value indicates a starting
%         position on the right-hand side of the runway (in flying
%         direction).
% ---optioncal (but in this order!): If omitted, the munich airport will be
%                                   used as a default
% - rw_lat: geodetic latitude of the runway reference point in degrees
% - rw_long: geodetic longitude of the runway reference point in degrees
% - rw_heading: Orientation of the runway relative to the NED- System. The 
%               angle is given in rad.
%
% OUTPUT:
% - lat: geodetic latitude of the starting point in degrees
% - long: geodetic longitude of the starting point in degrees
%
%-------------------------------------------------------------------------
%
% Author: Jens Dodenh�ft                          
%
%-------------------------------------------------------------------------
%
% CHANGE_LOG:
% - 03-Mar-2014: Initial Version
%
%%------------------------------------------------------------------------
function [lat,lon] = Calc_Starting_point_in_WGS84(x_rw,y_rw,varargin)

%Get Airport location and its runway orientation:
if(length(varargin)==3)
   %The function was called with airport location and orientation. Read
   %these inputs:    
   airport_lat=varargin{1};
   airport_long=varargin{2};
   rw_heading=varargin{3};
else
    %Set Munich airport as default:
    airport_lat= 48.36315;  %old:48.36295;
    airport_long=11.7726; %old:11.77;
    rw_heading=83.5*(pi/180);
end

%For test reasons:
% fprintf('Inputs in varargin(%d):\n',length(varargin));
% %fprintf('   %d\n', varargin{1});
% %fprintf('   %d\n', varargin{2});    
% %fprintf('   %d\n', varargin{3});
% fprintf('Airport values:');
% fprintf('lat:   %d\n', airport_lat);
% fprintf('long:   %d\n', airport_long);
% fprintf('runway heading:   %d\n', rw_heading);


%Build desired position in runway coordinate system:
position_runway_system=[-x_rw;y_rw;0];
%Transformation: Runway System --> Flat Earth System on Runway
T=[cos(rw_heading) sin(rw_heading) 0;-sin(rw_heading) cos(rw_heading) 0;0 0 1]';
%Keep in mind that rotation matrices are orthogonal matrices. Thus, the
%inverse of the transformation M_RW_O (that is the needed M_O_RW) can be
%determined by transposing the matrix. This is applied here.

position_flat_earth=T*position_runway_system;
%Calculate latitude and longitude:
WGS_coordinates=flat2lla((position_flat_earth)',[airport_lat airport_long],0,0,'WGS84');
lat=WGS_coordinates(1);
lon=WGS_coordinates(2);
end

